/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2018 André Müller (muellan@uni-mainz.de)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *****************************************************************************/

#include <sstream>
#include <unordered_map>
#include <algorithm> // sort
#include <set>
#include <unordered_set>

#include "dna_encoding.h"
#include "printing.h"
#include "sequence_io.h"
#include "sequence_view.h"
#include "alignment.h"
#include "candidates.h"
#include "querying.h"
#include "matches_per_target.h"

#include "classification.h"


namespace mc {

using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::flush;


/*************************************************************************//**
 *
 * @brief makes non-owning view to (sub)sequence
 *
 *****************************************************************************/
template<class Sequence>
inline auto
make_view_from_window_range(const Sequence& s, const window_range& range,
                            int size, int stride)
    -> decltype(make_view(s.begin(),s.end()))
{
    auto end = s.begin() + (stride * range.end) + size;
    if(end > s.end()) end = s.end();

    return make_view(s.begin() + (stride * range.beg), end);
}



/*************************************************************************//**
 *
 * @brief performs a semi-global alignment
 *
 *****************************************************************************/
template<class Subject>
alignment<default_alignment_scheme::score_type,typename sequence::value_type>
make_semi_global_alignment(const sequence_query& query,
                           const Subject& subject)
{
    std::size_t score  = 0;
    std::size_t scorer = 0;

    const auto scheme = default_alignment_scheme{};

    //compute alignment
    auto align = align_semi_global(query.seq1, subject, scheme);
    score = align.score;
    //reverse complement
    auto query1r = make_reverse_complement(query.seq1);
    auto alignr = align_semi_global(query1r, subject, scheme);
    scorer = alignr.score;

    //align paired read as well
    if(!query.seq2.empty()) {
        score += align_semi_global_score(query.seq2, subject, scheme);
        auto query2r = make_reverse_complement(query.seq2);
        scorer += align_semi_global_score(query2r, subject, scheme);
    }

    return (score > scorer) ? align : alignr;
}


/*************************************************************************//**
 *
 * @brief returns query taxon (ground truth for precision tests)
 *
 *****************************************************************************/
const taxon*
ground_truth(const database& db, const string& header)
{
    //try to extract query id and find the corresponding target in database
    const taxon* tax = nullptr;
    tax = db.taxon_with_name(extract_ncbi_accession_version_number(header));
    if(tax) return db.next_ranked_ancestor(tax);

    tax = db.taxon_with_name(extract_ncbi_accession_number(header));
    if(tax) return db.next_ranked_ancestor(tax);

    //try to extract id from header
    tax = db.taxon_with_id(extract_taxon_id(header));
    if(tax) return db.next_ranked_ancestor(tax);

    //try to find entire header as sequence identifier
    tax = db.taxon_with_name(header);
    if(tax) return db.next_ranked_ancestor(tax);

    return nullptr;
}



/*************************************************************************//**
 *
 * @brief removes hits of a specific taxon (on a rank) from database matches;
 *        can be very slow!
 *
 *****************************************************************************/
void remove_hits_on_rank(const database& db,
                         const taxon& tax, taxon_rank rank,
                         matches_per_location& hits)
{
    const taxon* excl = db.ancestor(tax,rank);

    matches_per_location maskedHits;
    for(const auto& hit : hits) {
        auto t = db.ancestor(hit.loc.tax, rank);
        if(t != excl) {
            maskedHits.push_back(hit);
        }
    }

    hits.swap(maskedHits);
}



/*************************************************************************//**
 *
 * @brief prepares ground truth based evaluation;
 *        might remove hits from 'allhits' if clade exclusion is turned on
 *
 *****************************************************************************/
void prepare_evaluation(const database& db,
                        const evaluation_options& opt,
                        sequence_query& query,
                        matches_per_location& allhits)
{
    if(opt.precision ||
       (opt.determineGroundTruth) ||
       (opt.excludeRank != taxon_rank::none) )
    {
        query.groundTruth = ground_truth(db, query.header);
    }

    //clade exclusion
    if(opt.excludeRank != taxon_rank::none && query.groundTruth) {
        remove_hits_on_rank(db, *query.groundTruth, opt.excludeRank, allhits);
    }
}



/*************************************************************************//**
 *
 * @brief lowest common ancestral taxon of several candidate taxa
 *
 *****************************************************************************/
const taxon*
lowest_common_ancestor(const database& db,
                       const classification_options& opt,
                       const classification_candidates& cand)
{
    using hit_count = match_candidate::count_type;

    if(cand.empty() || !cand[0].tax) return nullptr;

    if(cand.size() == 1) {
        return (cand[0].tax->rank() <= opt.highestRank) ? cand[0].tax : nullptr;
    }

    if(cand.size() == 2) {
        // (cand[1].hits > cand[0].hits - opt.hitsMin)
        const taxon* tax = db.ranked_lca(cand[0].tax, cand[1].tax);

        //classify if rank is below or at the highest rank of interest
        return (tax && tax->rank() <= opt.highestRank) ? tax : nullptr;
    }

    // begin lca with first candidate
    const taxon* lca_taxon = cand[0].tax;
    hit_count threshold = (cand[0].hits - opt.hitsMin) * opt.hitsDiffFraction;

    for(auto i = cand.begin()+1; i != cand.end(); ++i) {
        // include all candidates with hits above threshold
        if(i->hits > threshold) {
            // include candidate in lca
            lca_taxon = db.ranked_lca(lca_taxon, i->tax);
            // exit early if lca rank already too high
            if(!lca_taxon || lca_taxon->rank() > opt.highestRank)
                return nullptr;
        } else {
            break;
        }
    }
    return lca_taxon;
}



/*************************************************************************//**
 *
 * @brief classification cadidates + derived best classification
 *
 *****************************************************************************/
struct classification
{
    classification(classification_candidates cand):
        candidates{std::move(cand)}, best{nullptr}
    {}

    // classification(classification&& c):
    //     candidates{std::move(c.candidates)}, best{c.best}
    // {}

    classification_candidates candidates;
    const taxon* best;
};



/*************************************************************************//**
 *
 * @brief generate classification candidates
 *
 *****************************************************************************/
classification_candidates
make_classification_candidates(const database& db,
                               const classification_options& opt,
                               const sequence_query& query,
                               const matches_per_location& allhits)
{
    candidate_generation_rules rules;

    rules.maxWindowsInRange = window_id( 2 + (
        std::max(query.seq1.size() + query.seq2.size(), opt.insertSizeMax) /
        db.target_window_stride() ));

    rules.mergeBelow    = opt.lowestRank;
    rules.maxCandidates = opt.maxNumCandidatesPerQuery;

    return classification_candidates{db, allhits, rules};

}



/*************************************************************************//**
 *
 * @brief  classify using top matches/candidates
 *
 *****************************************************************************/
const taxon*
classify(const database& db, const classification_options& opt,
         const classification_candidates& cand)
{
    if(cand.empty()) return nullptr;

    if(cand.size() == 1) {
        return (cand[0].hits >= opt.hitsMin) ? cand[0].tax : nullptr;
    }

    //sum of top-2 hits < threshold => considered not classifiable
    if((cand[0].hits + cand[1].hits) < opt.hitsMin) return nullptr;

    //either top 2 are the same sequences with at least 'hitsMin' many hits
    //(checked before) or hit difference between these top 2 is above threshold
    if( (cand[0].tax == cand[1].tax)
        || (cand[0].hits - cand[1].hits >= opt.hitsMin) )
    {
        //return top candidate
        return cand[0].tax;
    }

    return lowest_common_ancestor(db, opt, cand);
}



/*************************************************************************//**
 *
 * @brief classify using all database matches
 *
 *****************************************************************************/
classification
classify(const database& db,
         const classification_options& opt,
         const sequence_query& query,
         const matches_per_location& allhits)
{
    classification cls { make_classification_candidates(db, opt, query, allhits) };

    cls.best = classify(db, opt, cls.candidates);

    return cls;
}



/*************************************************************************//**
 *
 * @brief classify using only matches in tgtMatches
 *
 *****************************************************************************/
void
update_classification(const database& db,
                      const classification_options& opt,
                      classification& cls,
                      const matches_per_target& tgtMatches)
{
    // size_t candidateCount = 0;
    // for(auto it = cls.candidates.begin(); it != cls.candidates.end(); ++it) {
    //     if(tgtMatches.find(it->tax) != tgtMatches.end())
    //         candidateCount += 1;
    // }
    // if(candidateCount == 0)
    //     //no candidates would remain -> keep old classification
    //     return;
    // if(candidateCount == cls.candidates.size())
    //     //nothing changed -> keep old classification
    //     return;

    for(auto it = cls.candidates.begin(); it != cls.candidates.end();) {
        if(tgtMatches.find(it->tax) == tgtMatches.end())
            it = cls.candidates.erase(it);
        else
            ++it;
    }
    cls.best = classify(db, opt, cls.candidates);
}


/*************************************************************************//**
 *
 * @brief  choose classification closest to ground truth from top candidates or lca
 *
 *****************************************************************************/
const taxon*
choose_best_classification(const database& db,
                           const classification_options& opt,
                           const taxon* groundTruth,
                           const classification_candidates& cand)
{
    if(!groundTruth)
        return nullptr;

    const taxon* best = nullptr;
    if(cand.size() >= 2) {
        best = db.ranked_lca(cand[0].tax, cand[1].tax);
    }
    auto lca = db.ranked_lca(best, groundTruth);
    auto lowestCorrectRank = lca ? lca->rank() : taxon_rank::none;
    
    for(const auto& c : cand) {
        const auto lca = db.ranked_lca(c.tax, groundTruth);
        if(lca && lca->rank() < lowestCorrectRank) {
            best = c.tax;
            lowestCorrectRank = lca->rank();
        }
    }

    // if(best && groundTruth && best->rank() != groundTruth->rank()) {
    lca = db.ranked_lca(best, groundTruth);
    if(lca && lca->rank() > opt.highestRank) {
        best = nullptr;
        lowestCorrectRank = taxon_rank::none;
    }

    return best;
}


/*************************************************************************//**
 *
 * @brief classify closest to ground truth using all database matches
 *
 *****************************************************************************/
classification
classify_best(const database& db,
              const classification_options& opt,
              const sequence_query& query,
              const matches_per_location& allhits)
{
    classification cls { make_classification_candidates(db, opt, query, allhits) };

    cls.best = choose_best_classification(db, opt, query.groundTruth, cls.candidates);

    return cls;
}


/*************************************************************************//**
 *
 * @brief add difference between result and truth to statistics
 *
 *****************************************************************************/
void update_coverage_statistics(const database& db,
                                const sequence_query& query,
                                const classification& cls,
                                classification_statistics& stats)
{
    if(!query.groundTruth) return;
    //check if taxa are covered in DB
    for(const taxon* tax : db.ranks(query.groundTruth)) {
        if(tax) {
            auto r = tax->rank();
            if(db.covers(*tax)) {
                if(!cls.best || r < cls.best->rank()) { //unclassified on rank
                    stats.count_coverage_false_neg(r);
                } else { //classified on rank
                    stats.count_coverage_true_pos(r);
                }
            }
            else {
                if(!cls.best || r < cls.best->rank()) { //unclassified on rank
                    stats.count_coverage_true_neg(r);
                } else { //classified on rank
                    stats.count_coverage_false_pos(r);
                }
            }
        }
    }
}



/*************************************************************************//**
 *
 * @brief evaluate classification of one query
 *
 *****************************************************************************/
void evaluate_classification(
    const database& db,
    const evaluation_options& opt,
    const sequence_query& query,
    const classification& cls,
    classification_statistics& statistics)
{
    if(opt.precision) {
        auto lca = db.ranked_lca(cls.best, query.groundTruth);
        auto lowestCorrectRank = lca ? lca->rank() : taxon_rank::none;

        statistics.assign_known_correct(
            cls.best ? cls.best->rank() : taxon_rank::none,
            query.groundTruth ? query.groundTruth->rank() : taxon_rank::none,
            lowestCorrectRank);

        //check if taxa of assigned target are covered
        if(opt.taxonCoverage) {
            update_coverage_statistics(db, query, cls, statistics);
        }

    } else {
        statistics.assign(cls.best ? cls.best->rank() : taxon_rank::none);
    }
}


/*************************************************************************//**
 *
 * @brief evaluate candidates of one query
 *
 *****************************************************************************/
void evaluate_candidates(
    const database& db,
    const evaluation_options& opt,
    const classification_options& opt_classify,
    const sequence_query& query,
    const classification& cls,
    classification_statistics& statistics)
{
    if(opt.precision) {
        // const taxon* best = cls.best;
        const taxon* best = nullptr;
        if(cls.candidates.size() >= 2) {
            best = db.ranked_lca(cls.candidates[0].tax, cls.candidates[1].tax);
        }
        auto lca = db.ranked_lca(best, query.groundTruth);
        auto lowestCorrectRank = lca ? lca->rank() : taxon_rank::none;
        
        for(const auto& cand : cls.candidates) {
            const auto lca = db.ranked_lca(cand.tax, query.groundTruth);
            if(lca && lca->rank() < lowestCorrectRank) {
                best = cand.tax;
                lowestCorrectRank = lca->rank();
            }
        }

        // if(best && query.groundTruth && best->rank() != query.groundTruth->rank()) {
        lca = db.ranked_lca(best, query.groundTruth);
        if(lca && lca->rank() > opt_classify.highestRank) {
            best = nullptr;
            lowestCorrectRank = taxon_rank::none;
        }

        statistics.assign_known_correct(
            best ? best->rank() : taxon_rank::none,
            query.groundTruth ? query.groundTruth->rank() : taxon_rank::none,
            lowestCorrectRank);

        //check if taxa of assigned target are covered
        // if(opt.taxonCoverage) {
            // update_coverage_statistics(db, query, cls, statistics);
        // }

    } else {
        statistics.assign(cls.best ? cls.best->rank() : taxon_rank::none);
    }
}


/*************************************************************************//**
 *
 * @brief compute alignment of top hits and optionally show it
 *
 *****************************************************************************/
void show_alignment(std::ostream& os,
                    const database& db,
                    const classification_output_options& opt,
                    const sequence_query& query,
                    const classification_candidates& tophits)
{
    //try to align to top target
    const taxon* tgtTax = tophits[0].tax;
    if(tgtTax && tgtTax->rank() == taxon_rank::Sequence) {
        const auto& src = tgtTax->source();
        try {
            //load candidate file and forward to sequence
            auto reader = make_sequence_reader(src.filename);
            reader->skip(src.index-1);

            if(reader->has_next()) {
                auto tgtSequ = reader->next().data;
                auto subject = make_view_from_window_range(
                               tgtSequ, tophits[0].pos,
                               db.target_window_size(),
                               db.target_window_stride());

                auto align = make_semi_global_alignment(query, subject);

                //print alignment to top candidate
                const auto w = db.target_window_stride();
                os  << '\n'
                    << opt.format.comment << "  score  " << align.score
                    << "  aligned to "
                    << src.filename << " #" << src.index
                    << " in range [" << (w * tophits[0].pos.beg)
                    << ',' << (w * tophits[0].pos.end + w) << "]\n"
                    << opt.format.comment << "  query  " << align.query << '\n'
                    << opt.format.comment << "  target " << align.subject;
            }
        }
        catch(std::exception& e) {
            if(opt.showErrors) cerr << e.what() << '\n';
        }
    }
}



/*************************************************************************//**
 *
 * @brief print header line for mapping table
 *
 *****************************************************************************/
void show_query_mapping_header(std::ostream& os,
                               const classification_output_options& opt)
{
    if(opt.mapViewMode == map_view_mode::none) return;

    const auto& colsep = opt.format.column;

    os << opt.format.comment << "TABLE_LAYOUT: ";

    if(opt.showQueryIds) os << "query_id" << colsep;

    os << "query_header" << colsep;

    if(opt.showGroundTruth) {
        show_taxon_header(os, opt, "truth_");
        os << colsep;
    }

    if(opt.showAllHits) os << "all_hits" << colsep;
    if(opt.showTopHits) os << "top_hits" << colsep;
    if(opt.showLocations) os << "candidate_locations" << colsep;

    show_taxon_header(os, opt);

    os << '\n';
}



/*************************************************************************//**
 *
 * @brief shows one query mapping line
 *        [query id], query_header, classification [, [top|all]hits list]
 *
 *****************************************************************************/
void show_query_mapping(
    std::ostream& os,
    const database& db,
    const classification_output_options& opt,
    const sequence_query& query,
    const classification& cls,
    const matches_per_location& allhits)
{
    if(opt.mapViewMode == map_view_mode::none ||
        (opt.mapViewMode == map_view_mode::mapped_only && !cls.best))
    {
        return;
    }

    const auto& colsep = opt.format.column;

    if(opt.showQueryIds) os << query.id << colsep;

    //print query header (first contiguous string only)
    auto l = query.header.find(' ');
    if(l != string::npos) {
        auto oit = std::ostream_iterator<char>{os, ""};
        std::copy(query.header.begin(), query.header.begin() + l, oit);
    }
    else {
        os << query.header;
    }
    os << colsep;

    if(opt.showGroundTruth) {
        show_taxon(os, db, opt, query.groundTruth);
        os << colsep;
    }

    if(opt.showAllHits) {
        show_matches(os, db, allhits, opt.lowestRank);
        os << colsep;
    }
    if(opt.showTopHits) {
        show_matches(os, db, cls.candidates, opt.lowestRank);
        os << colsep;
    }
    if(opt.showLocations) {
        show_candidate_ranges(os, db, cls.candidates);
        os << colsep;
    }

    show_taxon(os, db, opt, cls.best);

    if(opt.showAlignment && cls.best) {
        show_alignment(os, db, opt, query, cls.candidates);
    }

    os << '\n';
}



/*************************************************************************//**
 *
 * @brief todo
 *
 *****************************************************************************/
void show_read_features(
    std::ostream& os,
    const database& db,
    const classification_output_options& opt,
    const sequence_query& query,
    const classification& cls,
    const matches_per_location& allhits)
{
    // if(opt.mapViewMode == map_view_mode::none ||
    //     (opt.mapViewMode == map_view_mode::mapped_only && !cls.best))
    // {
    //     return;
    // }

    // const auto& colsep = opt.format.column;

    // if(opt.showQueryIds) os << query.id << colsep;

    // //print query header (first contiguous string only)
    // auto l = query.header.find(' ');
    // if(l != string::npos) {
    //     auto oit = std::ostream_iterator<char>{os, ""};
    //     std::copy(query.header.begin(), query.header.begin() + l, oit);
    // }
    // else {
    //     os << query.header;
    // }
    // os << colsep;

    // if(opt.showGroundTruth) {
    //     show_taxon(os, db, opt, query.groundTruth);
    //     os << colsep;
    // }

    // if(opt.showAllHits) {
    //     show_matches(os, db, allhits, opt.lowestRank);
    //     os << colsep;
    // }
    // if(opt.showTopHits) {
    //     show_matches(os, db, cls.candidates, opt.lowestRank);
    //     os << colsep;
    // }
    // if(opt.showLocations) {
    //     show_candidate_ranges(os, db, cls.candidates);
    //     os << colsep;
    // }

    // show_taxon(os, db, opt, cls.best);

    // if(opt.showAlignment && cls.best) {
    //     show_alignment(os, db, opt, query, cls.candidates);
    // }

    // os << '\n';
}



/*************************************************************************//**
 *
 * @brief filter out targets which have a coverage percentage below a percentile of all
 *        coverage percentages
 *
 *****************************************************************************/
void filter_targets_by_grouped_coverage(
    const database& db,
    matches_per_target& tgtMatches,
    const float percentile)
{
    using coverage_percentage = std::pair<const taxon*, float>;

    std::unordered_map<const taxon*, window_id> group_numwindows(tgtMatches.size());
    std::unordered_map<const taxon*, window_id> group_sizes(tgtMatches.size());
    std::unordered_map<const taxon*, vector<const taxon*> > group_children(tgtMatches.size());

    //sum number of coverred windows and sizes for each parent
    for(const auto& mapping : tgtMatches) {
        const taxon* target = mapping.first;
        const window_id targetSize = target->source().windows;
        std::set<window_id> hitWindows;
        for(const auto& candidate : mapping.second) {
            for(const auto& windowMatch : candidate.matches) {
                hitWindows.emplace(windowMatch.win);
            }
        }
        const taxon* parent = db.parent(target);
        group_numwindows[parent] += hitWindows.size();
        group_sizes[parent] += targetSize;
        group_children[parent].emplace_back(target);
    }

    //calculate coverage percentages
    vector<coverage_percentage> coveragePercentages;
    coveragePercentages.reserve(tgtMatches.size());

    float coveragePercentagesSum = 0;
    for(const auto& group_size : group_sizes) {
        const float covP = float(group_numwindows[group_size.first]) / group_size.second;
        coveragePercentagesSum += covP;
        coveragePercentages.emplace_back(group_size.first, covP);
    }

    //sort by coverage descending
    std::sort(coveragePercentages.begin(), coveragePercentages.end(),
        [](coverage_percentage& a, coverage_percentage& b){
            return a.second < b.second;
    });

    //filter out targets
    float coveragePercentagesPartSum = 0;
    for(auto it = coveragePercentages.begin(); it != coveragePercentages.end(); ++it) {
        coveragePercentagesPartSum += it->second;
        if(coveragePercentagesPartSum > percentile / 100 * coveragePercentagesSum) {
            // coveragePercentages.resize(std::distance(it, coveragePercentages.rend()));
            break;
        }
        for(const auto& child : group_children[it->first])
            tgtMatches.erase(child);
    }
}


/*************************************************************************//**
 *
 * @brief filter out targets which have a larger gap than maxGapSize
 *
 *****************************************************************************/
void filter_targets_by_gap(
    matches_per_target& tgtMatches,
    const window_id maxGapSize)
{
    for(auto it = tgtMatches.begin(); it != tgtMatches.end();) {
        const taxon* target = it->first;
        std::set<window_id> hitWindows;
        for(const auto& candidate : it->second) {
            for(const auto& windowMatch : candidate.matches) {
                hitWindows.emplace(windowMatch.win);
            }
        }
        const window_id targetSize = target->source().windows;

        window_id longestGapSize = 0;
        window_id prevWin = 0;
        for(const auto& win : hitWindows) {
            const window_id gapSize = win - prevWin;
            if(gapSize > longestGapSize)
                longestGapSize = gapSize;
            prevWin = win;
        }
        const window_id gapSize = (targetSize-1) - prevWin;
        if(gapSize > longestGapSize)
            longestGapSize = gapSize;

        if(longestGapSize > maxGapSize)
            it = tgtMatches.erase(it);
        else
            ++it;
    }
}


/*************************************************************************//**
 *
 * @brief filter out targets which have a coverage percentage below a percentile of all
 *        coverage percentages                 
 *
 *****************************************************************************/
void filter_targets_by_coverage(
    matches_per_target& tgtMatches,
    const float percentile)
{
    using coverage_percentage = std::pair<const taxon*, float>;

    vector<coverage_percentage> coveragePercentages;
    coveragePercentages.reserve(tgtMatches.size());

    float coveragePercentagesSum = 0;

    //calculate coverage percentages
    for(const auto& mapping : tgtMatches) {
        const taxon* target = mapping.first;
        const window_id targetSize = target->source().windows;
        std::unordered_set<window_id> hitWindows;
        for(const auto& candidate : mapping.second) {
            for(const auto& windowMatch : candidate.matches) {
                hitWindows.emplace(windowMatch.win);
            }
        }
        const float covP = float(hitWindows.size()) / targetSize;
        coveragePercentagesSum += covP;
        coveragePercentages.emplace_back(target, covP);
    }

    //sort by coverage descending
    std::sort(coveragePercentages.begin(), coveragePercentages.end(),
        [](coverage_percentage& a, coverage_percentage& b){
            return a.second < b.second;
    });

    //filter out targets
    float coveragePercentagesPartSum = 0;
    for(auto it = coveragePercentages.begin(); it != coveragePercentages.end(); ++it) {
        coveragePercentagesPartSum += it->second;
        if(coveragePercentagesPartSum > percentile / 100 * coveragePercentagesSum) {
            // coveragePercentages.resize(std::distance(it, coveragePercentages.rend()));
            break;
        }
        tgtMatches.erase(it->first);
    }
}

/*************************************************************************//**
 *
 * @brief filter out targets
 *        - which have no uniquely mapped reads
 *        - which have a coverage percentage below a percentile of all
 *          coverage percentages                 
 *
 *****************************************************************************/
void filter_targets_by_coverage(
    matches_per_target& uniqueTgtMatches, 
    matches_per_target& tgtMatches,
    const float percentile)
{
    if(!uniqueTgtMatches.empty()) {
        //filter unique matches by coverage
        filter_targets_by_coverage(uniqueTgtMatches, percentile);

        //discard targets without uniques
        for(auto it = tgtMatches.begin(); it != tgtMatches.end();) {
            if(uniqueTgtMatches.find(it->first) == uniqueTgtMatches.end())
                it = tgtMatches.erase(it);
            else
                ++it;
        }
    }

    //filter all matches by coverage
    filter_targets_by_coverage(tgtMatches, percentile);
}

/*************************************************************************//**
 *
 * @brief per-batch buffer for output and (target -> hits) lists
 *
 *****************************************************************************/
struct query_mapping
{
    sequence_query query;
    classification cls;
    // matches_per_location allhits;
};

/*************************************************************************//**
 *
 * @brief per-batch buffer for output and (target -> hits) lists
 *
 *****************************************************************************/
struct mappings_buffer
{
    std::ostringstream out;
    std::unordered_map<query_id, query_mapping> queryMappings;
    matches_per_target hitsPerTarget;
    matches_per_target uniqueHitsPerTarget;
};

/*************************************************************************//**
 *
 * @brief default classification scheme with
 *        additional target->hits list generation and output
 *        try to map each read to a taxon with the lowest possible rank
 *
 *****************************************************************************/
void map_queries_to_targets_default(
    const vector<string>& infiles,
    const database& db, const query_options& opt,
    classification_results& results)
{
    //global target -> query_id/win:hits... list
    matches_per_target uniqueTgtMatches;
    matches_per_target tgtMatches;
    std::unordered_map<query_id, query_mapping> queryMappings;

    //makes an empty batch buffer
    //(each batch might be processed by a different thread)
    const auto makeBatchBuffer = [] { return mappings_buffer(); };

    //updates buffer with the database answer of a single query
    const auto processQuery = [&] (mappings_buffer& buf,
        sequence_query&& query, matches_per_location& allhits)
    {
        if(query.empty()) return;

        prepare_evaluation(db, opt.evaluate, query, allhits);

        // auto cls = classify(db, opt.classify, query, allhits);
        auto cls = classify_best(db, opt.classify, query, allhits);

        if(opt.output.showHitsPerTargetList) {
            //insert best candidate if classified as unique sequence
            // if(cls.best && cls.best->rank() == taxon_rank::Sequence)
            //     buf.uniqueHitsPerTarget.insert(query.id, allhits, cls.candidates,
            //                                    cls.candidates[0].hits - 1);
            //insert all candidates relevant for classification
            // if(cls.best && cls.best->rank() != taxon_rank::none) // this only works if all targets are ranked 
            if((cls.candidates.size() == 1 && cls.candidates[0].hits >= opt.classify.hitsMin) ||
               (cls.candidates.size() >= 2 && (cls.candidates[0].hits + cls.candidates[1].hits) >= opt.classify.hitsMin))
                buf.hitsPerTarget.insert(query.id, allhits, cls.candidates,
                                         // (cls.candidates[0].hits - opt.classify.hitsMin));
                                         // (cls.candidates[0].hits - opt.classify.hitsMin) * opt.classify.hitsDiffFraction);
                                         0);
        }
        else {
            // evaluate_classification(db, opt.evaluate, query, cls, results.statistics);
            evaluate_candidates(db, opt.evaluate, opt.classify, query, cls, results.statistics);

            show_query_mapping(buf.out, db, opt.output, query, cls, allhits);
        }

        sequence().swap(query.seq1);
        // query.seq1.swap(sequence);
        if(query.seq2.size())
            sequence().swap(query.seq2);
            // query.seq2.swap(sequence{});

        buf.queryMappings.emplace(query.id, query_mapping{std::move(query), std::move(cls)});
    };


    const auto finalizeBatch = [&] (mappings_buffer&& buf) {
        if(opt.output.showHitsPerTargetList) {
            queryMappings.insert(std::make_move_iterator(buf.queryMappings.begin()),
                                 std::make_move_iterator(buf.queryMappings.end()) );
            //merge batch (target->hits) lists into global one
            uniqueTgtMatches.merge(std::move(buf.uniqueHitsPerTarget));
            tgtMatches.merge(std::move(buf.hitsPerTarget));
        }
        //write output buffer to output stream when batch is finished
        results.mapout << buf.out.str();
    };


    const auto showStatus = [&] (const std::string& msg, float progress) {
        if(opt.output.mapViewMode != map_view_mode::none) {
            results.mapout << opt.output.format.comment << msg << '\n';
        }
        if(progress >= 0.0f) {
            show_progress_indicator(results.status, progress);
        }
    };

    //run (parallel) database queries according to processing options
    query_database(infiles, db, opt.process,
                   makeBatchBuffer, processQuery, finalizeBatch,
                   showStatus);

    //target -> hits list?
    if(tgtMatches.empty()) return;

    //filter uniqueTgtMatches, delete targets from uniqueTgtMatches and tgtMatches
    //filter tgtMatches, delete targets from tgtMatches
    // filter_targets_by_coverage(uniqueTgtMatches, tgtMatches);
    // if(!uniqueTgtMatches.empty()) {
    //     //filter unique matches by coverage
    //     // filter_targets_by_coverage(uniqueTgtMatches, percentile);

    //     //discard targets without uniques
    //     for(auto it = tgtMatches.begin(); it != tgtMatches.end();) {
    //         if(uniqueTgtMatches.find(it->first) == uniqueTgtMatches.end())
    //             it = tgtMatches.erase(it);
    //         else
    //             ++it;
    //     }
    // }
    //filter all matches by coverage
    if(opt.classify.covPercentile > 0)
        filter_targets_by_coverage(tgtMatches, opt.classify.covPercentile);
    // filter_targets_by_grouped_coverage(db, tgtMatches, opt.classify.covPercentile);

    if(opt.classify.maxGapSize < std::numeric_limits<window_id>::max())
        filter_targets_by_gap(tgtMatches, opt.classify.maxGapSize);

    //serial
    // for(auto& mapping : queryMappings) {
    //     // classify using only targets left in tgtMatches
    //     update_classification(db, opt.classify, mapping.second.cls, tgtMatches);

    //     evaluate_classification(db, opt.evaluate, mapping.second.query, mapping.second.cls, results.statistics);

    //     show_query_mapping(results.mapout, db, opt.output, mapping.second.query, mapping.second.cls, matches_per_location{});
    // }

    //parallel
    parallel_queue queue(opt.process.numThreads);
    const auto load = 32 * queue.concurrency();
    std::mutex mtx;

    auto it = queryMappings.begin();
    while(it != queryMappings.end()) {
        if(queue.unsafe_waiting() < load) {
            auto batchEnd = size_t(std::distance(it, queryMappings.end())) > opt.process.batchSize ?
                            std::next(it, opt.process.batchSize) : queryMappings.end();
            queue.enqueue([&, it, batchEnd] {
                std::ostringstream bufout;
                for(auto i = it; i != batchEnd; ++i) {
                    // classify using only targets left in tgtMatches
                    update_classification(db, opt.classify, i->second.cls, tgtMatches);

                    evaluate_classification(db, opt.evaluate, i->second.query, i->second.cls, results.statistics);

                    show_query_mapping(bufout, db, opt.output, i->second.query, i->second.cls, matches_per_location{});
                }
                std::lock_guard<std::mutex> lock(mtx);
                results.mapout << bufout.str();
            }); //enqueue
            it = batchEnd;
        }
    }
    //wait for all enqueued tasks to finish
    queue.wait(); 


    // results.auxout
    //     << opt.output.format.comment
    //     << "--- list of hits for each reference sequence ---\n"
    //     << opt.output.format.comment
    //     << "window start position within sequence = window_index * window_stride(="
    //     << db.query_window_stride() << ")\n";

    // tgtMatches.sort_match_lists();
    // show_matches_per_targets(results.auxout, db, tgtMatches, opt.output);
    // show_num_matches_per_targets(results.auxout, db, tgtMatches, opt.output);
    show_features_of_targets(results.auxout, db, tgtMatches, opt.output);
}



/*************************************************************************//**
 *
 * @brief default classification scheme & output
 *        try to map each read to a taxon with the lowest possible rank
 *
 *****************************************************************************/
void map_queries_to_targets(const vector<string>& infiles,
                            const database& db, const query_options& opt,
                            classification_results& results)
{
    if(opt.output.mapViewMode != map_view_mode::none) {
        show_query_mapping_header(results.mapout, opt.output);
    }

    map_queries_to_targets_default(infiles, db, opt, results);
}

} // namespace mc

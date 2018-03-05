/*****************************************************************************
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

#include "dna_encoding.h"
#include "printing.h"
#include "sequence_io.h"
#include "sequence_view.h"
#include "alignment.h"
#include "candidates.h"
#include "query_options.h"
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


/*****************************************************************************
 *
 * @brief makes non-owning view to (sub)sequence
 *
 *****************************************************************************/
template<class Sequence, class Index>
inline auto
make_view_from_window_range(
    const Sequence& s, const window_range<Index>& range, int size, int stride)
    -> decltype(make_view(s.begin(),s.end()))
{
    auto end = s.begin() + (stride * range.end) + size;
    if(end > s.end()) end = s.end();

    return make_view(s.begin() + (stride * range.beg), end);
}



/*****************************************************************************
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


/*****************************************************************************
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



/*****************************************************************************
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



/*****************************************************************************
 *
 * @brief
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



/*****************************************************************************
 *
 * @brief lowest common taxon of several classification candidate
 *        sequences / taxa
 *
 * @param trustedMajority  fraction of total feature hits that a
 *                         candidate (taxon) must have in order to be
 *                         chosen as classification result
 *
 * @param lowestRank       rank to start investigation on
 * @param highestRank
 *
 *****************************************************************************/
const taxon*
lowest_common_taxon(
    const database& db,
    const classification_candidates& cand,
    float trustedMajority,
    taxon_rank lowestRank = taxon_rank::Sequence,
    taxon_rank highestRank = taxon_rank::Domain)
{
    if(!cand[0].tax) return nullptr;
    if(cand.size() < 2) {
        if(cand[0].tax->rank() <= highestRank) return cand[0].tax;
    }
    else if(cand.size() < 3) {
        const taxon* tax = db.ranked_lca(cand[0].tax, cand[1].tax);

        //classify if rank is below or at the highest rank of interest
        if(tax && tax->rank() <= highestRank) return tax;
    }
    else {
        if(lowestRank == taxon_rank::Sequence) ++lowestRank;

        std::unordered_map<const taxon*,int> scores;
        scores.rehash(2*cand.size());

        int totalscore = 0;
        for(int i = 0, n = cand.size(); i < n; ++i) {
            totalscore += cand[i].hits;
        }
        float threshold = totalscore * trustedMajority;

        for(auto rank = lowestRank; rank <= highestRank; ++rank) {
            //hash-count taxon id occurrences on rank 'r'
            for(int i = 0, n = cand.size(); i < n; ++i) {
                //use target id instead of taxon if at sequence level
                const taxon* tax = db.ancestor(cand[i].tax, rank);
                auto score = cand[i].hits;
                auto it = scores.find(tax);
                if(it != scores.end()) {
                    it->second += score;
                } else {
                    scores.insert(it, {tax, score});
                }
            }
            //determine taxon with most votes
            const taxon* toptax = nullptr;
            int topscore = 0;
            for(const auto& x : scores) {
                if(x.second > topscore) {
                    toptax   = x.first;
                    topscore = x.second;
                }
            }
            //if enough candidates (weighted by their hits)
            //agree on a taxon => classify as such
            if(toptax && topscore >= threshold) {
                return toptax;
            }
            scores.clear();
        }
    }
    //candidates couldn't agree on a taxon on any relevant taxonomic rank
    return nullptr;
}



/*************************************************************************//**
 *
 * @brief
 *
 *****************************************************************************/
struct classification
{
    classification(classification_candidates cand):
        candidates{std::move(cand)}, best{nullptr}
    {}

    classification_candidates candidates;
    const taxon* best;
};



/*****************************************************************************
 *
 * @brief  classify using top candidates
 *
 *****************************************************************************/
const taxon*
classify(const database& db, const classification_options& opt,
         const classification_candidates& cand)
{
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

    return lowest_common_taxon(db, cand, opt.hitsDiffFraction,
                               opt.lowestRank, opt.highestRank);
}



/*************************************************************************//**
 *
 * @brief
 *
 *****************************************************************************/
classification
classify(const database& db,
         const classification_options& opt,
         const sequence_query& query,
         const matches_per_location& allhits)
{
    auto numWindows = window_id( 2 + (
        std::max(query.seq1.size() + query.seq2.size(), opt.insertSizeMax) /
        db.target_window_stride() ));

    classification cls {
        classification_candidates{db, allhits, numWindows, opt.lowestRank} };

    cls.best = classify(db, opt, cls.candidates);

    return cls;
}



/*****************************************************************************
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



/*****************************************************************************
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
            reader->skip(src.index);

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
                    << opt.comment << "  score  " << align.score
                    << "  aligned to "
                    << src.filename << " #" << src.index
                    << " in range [" << (w * tophits[0].pos.beg)
                    << ',' << (w * tophits[0].pos.end + w) << "]\n"
                    << opt.comment << "  query  " << align.query << '\n'
                    << opt.comment << "  target " << align.subject;
            }
        }
        catch(std::exception& e) {
            if(opt.showErrors) cerr << e.what() << '\n';
        }
    }
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

    if(opt.showQueryIds) {
        os << query.id << opt.separator;
    }

    //print query header (first contiguous string only)
    auto l = query.header.find(' ');
    if(l != string::npos) {
        auto oit = std::ostream_iterator<char>{os, ""};
        std::copy(query.header.begin(), query.header.begin() + l, oit);
    }
    else {
        os << query.header;
    }
    os << opt.separator;

    if(opt.showGroundTruth) {
        show_taxon(os, db, opt, query.groundTruth);
        os << opt.separator;
    }

    if(opt.showAllHits) {
        show_matches(os, db, allhits, opt.lowestRank);
        os << opt.separator;
    }
    if(opt.showTopHits) {
        show_matches(os, db, cls.candidates, opt.lowestRank);
        os << opt.separator;
    }
    if(opt.showLocations) {
        show_candidate_ranges(os, db, cls.candidates);
        os << opt.separator;
    }

    show_taxon(os, db, opt, cls.best);

    if(opt.showAlignment && cls.best) {
        show_alignment(os, db, opt, query, cls.candidates);
    }

    os << '\n';
}



/*************************************************************************//**
 *
 * @brief default classification scheme & output
 *        try to map each read to a taxon with the lowest possible rank
 *
 *****************************************************************************/
void map_queries_to_targets_default(
    const vector<string>& infiles,
    const database& db, const query_options& opt,
    classification_results& results)
{
    using obuf_t = std::ostringstream;

    const auto makeBatchBuffer = [] { return obuf_t(); };


    const auto processQuery = [&] (
        obuf_t& buf, sequence_query&& query, matches_per_location&& allhits)
    {
        if(query.empty()) return;

        prepare_evaluation(db, opt.evaluate, query, allhits);

        auto cls = classify(db, opt.classify, query, allhits);

        evaluate_classification(db, opt.evaluate, query, cls, results.statistics);

        show_query_mapping(buf, db, opt.output, query, cls, allhits);
    };


    const auto finalizeBatch = [&] (const obuf_t& buf) {
        //write buffer to output stream when batch is finished
        results.out << buf.str();
    };


    const auto showStatus = [&] (const std::string& msg, float progress) {
        if(opt.output.mapViewMode != map_view_mode::none) {
            results.out << opt.output.comment << msg << '\n';
        }
        if(progress >= 0.0f) {
            show_progress_indicator(results.status, progress);
        }
    };


    //run (parallel) database queries according to processing options
    query_database(infiles, db, opt.process,
                   makeBatchBuffer, processQuery, finalizeBatch,
                   showStatus);
}



/*************************************************************************//**
 *
 * @brief per-batch buffer for output and (target -> hits) lists
 *
 *****************************************************************************/
struct mappings_buffer
{
    std::ostringstream out;
    matches_per_target hitsPerTarget;
};

/*************************************************************************//**
 *
 * @brief default classification scheme with
 *        additional target->hits list generation and output
 *        try to map each read to a taxon with the lowest possible rank
 *
 *****************************************************************************/
void map_queries_to_targets_and_vice_versa(
    const vector<string>& infiles,
    const database& db, const query_options& opt,
    classification_results& results)
{
    //global target -> query_id/win:hits... list
    matches_per_target tgtMatches;

    //batch buffer (each batch might be processed by a different thread)
    const auto makeBatchBuffer = [] { return mappings_buffer(); };


    const auto processQuery = [&] (mappings_buffer& buf,
        sequence_query&& query, matches_per_location&& allhits)
    {
        if(query.empty()) return;

        prepare_evaluation(db, opt.evaluate, query, allhits);

        auto cls = classify(db, opt.classify, query, allhits);

        buf.hitsPerTarget.insert(query.id, allhits, cls.candidates);

        evaluate_classification(db, opt.evaluate, query, cls, results.statistics);

        show_query_mapping(buf.out, db, opt.output, query, cls, allhits);
    };


    const auto finalizeBatch = [&] (mappings_buffer&& buf) {
        //merge batch (target->hits) lists into global one
        tgtMatches.merge(std::move(buf.hitsPerTarget));
        //write output buffer to output stream when batch is finished
        results.out << buf.out.str();
    };


    const auto showStatus = [&] (const std::string& msg, float progress) {
        if(opt.output.mapViewMode != map_view_mode::none) {
            results.out << opt.output.comment << msg << '\n';
        }
        if(progress >= 0.0f) {
            show_progress_indicator(results.status, progress);
        }
    };


    //run (parallel) database queries according to processing options
    query_database(infiles, db, opt.process,
                   makeBatchBuffer, processQuery, finalizeBatch,
                   showStatus);

    if(tgtMatches.empty()) return;

    std::ostream* tgtMatchesOut = &std::cout;
    std::ofstream tgtMatchesFile;
    if(!opt.output.hitsPerTargetOutfile.empty()) {
        tgtMatchesFile.open(opt.output.hitsPerTargetOutfile);
        if(tgtMatchesFile.good()) {
            results.out << opt.output.comment
                << " (genome -> matches) list will be written to "
                << opt.output.hitsPerTargetOutfile << '\n';
            tgtMatchesOut = &tgtMatchesFile;
        }
    }

    if(opt.output.mapViewMode != map_view_mode::none && !tgtMatchesFile->good()) {
        results.out << opt.output.comment
            << "------------------------------------------------------------\n";
    }

    *tgtMatchesOut  << opt.output.comment
        << " note: window start position within genome = window_index * window_stride(="
        << db.query_window_stride() << ")\n";

    tgtMatches.sort_match_lists();
    show_matches_per_targets(*tgtMatchesOut, tgtMatches, opt.output);
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
    if(opt.output.showHitsPerTargetList) {
        map_queries_to_targets_and_vice_versa(infiles, db, opt, results);
    }
    //default (and the only functionality available prior to version 0.20)
    else {
        map_queries_to_targets_default(infiles, db, opt, results);
    }
}

} // namespace mc

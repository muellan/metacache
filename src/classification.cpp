
#include "dna_encoding.h"
#include "print_info.h"
#include "print_results.h"
#include "sequence_io.h"
#include "sequence_view.h"
#include "alignment.h"
#include "candidates.h"
#include "query_options.h"

#include "classification.h"


namespace mc {

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
template<class Query, class Subject, class AlignmentScheme>
alignment<typename AlignmentScheme::score_type,typename Query::value_type>
make_alignment(const Query& query1, const Query& query2,
               const Subject& subject,
               const AlignmentScheme& scheme)
{
    std::size_t score  = 0;
    std::size_t scorer = 0;

    //compute alignment
    auto align = align_semi_global(query1, subject, scheme);
    score = align.score;
    //reverse complement
    auto query1r = make_reverse_complement(query1);
    auto alignr = align_semi_global(query1r, subject, scheme);
    scorer = alignr.score;

    //align paired read as well
    if(!query2.empty()) {
        score += align_semi_global_score(query2, subject, scheme);
        auto query2r = make_reverse_complement(query2);
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
    if(cand.size() < 3) {
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



/*****************************************************************************
 *
 * @brief  analyze hit statistics from database
 * @param  query1  read
 * @param  query2  paired read
 *
 *****************************************************************************/
const taxon*
classification(const database& db, const query_options& opt,
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



/*****************************************************************************
 *
 * @brief prints a classification to stream 'os'
 *
 *****************************************************************************/
void show_classification(std::ostream& os,
                         const database& db, const query_options& opt,
                         const taxon* tax)
{
    if(!tax || tax->rank() > opt.highestRank) {
        if(opt.separateTaxaInfo) {
            auto rmax = opt.showLineage ? opt.highestRank : opt.lowestRank;
            for(auto r = opt.lowestRank; r <= rmax; ++r) {
                switch(opt.showTaxaAs) {
                    default:
                    case taxon_print_mode::rank_name:
                        os << "--" << opt.outSeparator << "--";
                        break;
                    case taxon_print_mode::rank_id:
                        os << "--" << opt.outSeparator << taxonomy::none_id();
                        break;
                    case taxon_print_mode::rank_name_id:
                        os << "--" << opt.outSeparator << "--"
                           << opt.outSeparator << taxonomy::none_id();
                        break;
                }
                if(r < rmax) os << opt.outSeparator;
            }
        }
        else {
            os << "--";
        }
    }
    else {
        auto rmin = opt.lowestRank < tax->rank() ? tax->rank() : opt.lowestRank;
        auto rmax = opt.showLineage ? opt.highestRank : rmin;

        if(opt.separateTaxaInfo) {
            show_ranks(os, db.ranks(tax), opt.showTaxaAs, rmin, rmax, opt.outSeparator);
        }
        else {
            show_ranks(os, db.ranks(tax), opt.showTaxaAs, rmin, rmax, "");
        }
    }
}



/*****************************************************************************
 *
 * @brief add difference between result and truth to statistics
 *
 *****************************************************************************/
void update_coverage_statistics(const database& db,
                                const taxon* result, const taxon* truth,
                                classification_statistics& stats)
{
    if(!truth) return;
    //check if taxa are covered in DB
    for(const taxon* tax : db.ranks(truth)) {
        if(tax) {
            auto r = tax->rank();
            if(db.covers(*tax)) {
                if(!result || r < result->rank()) { //unclassified on rank
                    stats.count_coverage_false_neg(r);
                } else { //classified on rank
                    stats.count_coverage_true_pos(r);
                }
            }
            else {
                if(!result || r < result->rank()) { //unclassified on rank
                    stats.count_coverage_true_neg(r);
                } else { //classified on rank
                    stats.count_coverage_false_pos(r);
                }
            }
        }
    }
}



/*****************************************************************************
 *
 * @brief compute alignment of top hits and optionally show it
 *
 *****************************************************************************/
void use_alignment(const database& db, const query_options& opt,
                   const sequence& query1, const sequence& query2,
                   const classification_candidates& tophits,
                   classification_statistics& stats,
                   std::ostream& os)
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
                auto subject = make_view_from_window_range(tgtSequ,
                                                           tophits[0].pos,
                                                           db.target_window_size(),
                                                           db.target_window_stride());

                auto align = make_alignment(query1, query2, subject,
                                            needleman_wunsch_scheme{});

                stats.register_alignment_score(align.score);

                //print alignment to top candidate
                if(align.score > 0 && opt.showAlignment) {
                    const auto w = db.target_window_stride();

                    os << '\n'
                        << opt.comment << "  score  " << align.score
                        << "  aligned to "
                        << src.filename << " #" << src.index
                        << " in range [" << (w * tophits[0].pos.beg)
                        << ',' << (w * tophits[0].pos.end) << "]\n"
                        << opt.comment << "  query  " << align.query << '\n'
                        << opt.comment << "  target " << align.subject;
                }
            }
        } catch(std::exception& e) {
            if(opt.showErrors) cerr << e.what() << '\n';
        }
    }
}



/*****************************************************************************
 *
 * @brief process database hit results: classify & print mapping
 *
 *****************************************************************************/
void process_database_answer(
     const database& db, const query_options& opt,
     const string& header,
     const sequence& query1, const sequence& query2,
     matches_per_location&& hits,
     classification_statistics& stats, std::ostream& os)
{
    if(header.empty()) return;

    //preparation -------------------------------
    const taxon* groundTruth = nullptr;
    if(opt.testPrecision ||
       (opt.mapViewMode != map_view_mode::none && opt.showGroundTruth) ||
       (opt.excludedRank != taxon_rank::none) )
    {
        groundTruth = ground_truth(db, header);
    }
    //clade exclusion
    if(opt.excludedRank != taxon_rank::none && groundTruth) {
        remove_hits_on_rank(db, *groundTruth, opt.excludedRank, hits);
    }

    //classify ----------------------------------
    auto numWindows = window_id( 2 + (
        std::max(query1.size() + query2.size(), opt.insertSizeMax) /
        db.target_window_stride() ));

    classification_candidates tophits {db, hits, numWindows, opt.lowestRank};
    const taxon* cls = classification(db, opt, tophits );

    //evaluate classification -------------------
    if(opt.testPrecision) {
        auto lca = db.ranked_lca(cls, groundTruth);
        auto lowestCorrectRank = lca ? lca->rank() : taxon_rank::none;

        stats.assign_known_correct(
            cls ? cls->rank() : taxon_rank::none,
            groundTruth ? groundTruth->rank() : taxon_rank::none,
            lowestCorrectRank);

        //check if taxa of assigned target are covered
        if(opt.testCoverage) {
            update_coverage_statistics(db, cls, groundTruth, stats);
        }

    } else {
        stats.assign(cls ? cls->rank() : taxon_rank::none);
    }

    //show mapping ------------------------------
    bool showMapping = opt.mapViewMode == map_view_mode::all ||
        (opt.mapViewMode == map_view_mode::mapped_only && cls);

    if(showMapping) {
        //print query header (first contiguous string only)
        auto l = header.find(' ');
        if(l != string::npos) {
            auto oit = std::ostream_iterator<char>{os, ""};
            std::copy(header.begin(), header.begin() + l, oit);
        }
        else {
            os << header;
        }
        os << opt.outSeparator;

        if(opt.showGroundTruth) {
            if(groundTruth) {
                show_ranks(os, db.ranks(groundTruth),
                    opt.showTaxaAs, groundTruth->rank(),
                    opt.showLineage ? opt.highestRank : groundTruth->rank(),
                    opt.outSeparator);
            } else {
                os << "n/a";
            }
            os << opt.outSeparator;
        }

        //print results
        if(opt.showAllHits) {
            show_matches(os, db, hits, opt.lowestRank);
            os << opt.outSeparator;
        }
        if(opt.showTopHits) {
            show_matches(os, db, tophits, opt.lowestRank);
            os << opt.outSeparator;
        }
        if(opt.showLocations) {
            show_candidate_ranges(os, db, tophits);
            os << opt.outSeparator;
        }
        show_classification(os, db, opt, cls);
    }

    //optional alignment ------------------------------
    if(opt.testAlignment && cls) {
        use_alignment(db, opt, query1, query2, tophits, stats, os);
    }

    if(showMapping) os << '\n';
}


} // namespace mc

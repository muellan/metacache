/*****************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2017 André Müller (muellan@uni-mainz.de)
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

#include <functional>
#include <limits>

#include "args_handling.h"
#include "timer.h"
#include "filesys_utility.h"
#include "cmdline_utility.h"
#include "sequence_io.h"
#include "sequence_view.h"
#include "alignment.h"
#include "parallel_task_queue.h"
#include "candidates.h"
#include "print_info.h"
#include "print_results.h"


namespace mc {

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::flush;


/*****************************************************************************
 * @brief pairing of queries
 *****************************************************************************/
enum class pairing_mode : unsigned char {
    none, files, sequences
};


/*****************************************************************************
 * @brief alignment mode
 *****************************************************************************/
enum class align_mode : unsigned char {
    none, semi_global
};


/*****************************************************************************
 * @brief how to show mapping
 *****************************************************************************/
enum class map_view_mode : unsigned char {
    none, mapped_only, all
};



/*****************************************************************************
 *
 * @brief query options
 *
 *****************************************************************************/
struct query_options
{
    //-----------------------------------------------------
    //output options & formatting
    //-----------------------------------------------------
    //prefix for each non-mapping line
    string comment = "# ";
    //separates individual mapping fields
    string outSeparator = "\t|\t";
    //show database properties
    bool showDBproperties = false;
    //make a separate output file for each input file
    bool splitOutput = false;
    //show top candidate sequences and their associated k-mer hash hit count
    bool showTopHits = false;
    //show all k-mer-hash hits in database for each given read
    bool showAllHits = false;
    //how to show classification (read mappings), if 'none', only summary will be shown
    map_view_mode mapViewMode = map_view_mode::all;
    //show candidate position(s) in reference sequence(s)
    bool showLocations = false;
    //show known taxon (or complete lineage if 'showLineage' on)
    bool showGroundTruth = false;
    //show all ranks that a sequence could be classified on
    bool showLineage = false;
    //what to show of a taxon
    taxon_print_mode showTaxaAs = taxon_print_mode::rank_name;
    bool separateTaxaInfo = false;
    bool showAlignment = false;
    //show error messages?
    bool showErrors = true;

    //-----------------------------------------------------
    //classification options
    //-----------------------------------------------------
    pairing_mode pairing = pairing_mode::none;
    //ranks/taxa to classify on and to show
    taxon_rank lowestRank  = taxon_rank::Sequence;
    taxon_rank highestRank = taxon_rank::Domain;
    //ground truth rank to exclude (for clade exclusion test)
    taxon_rank excludedRank = taxon_rank::none;
    //
    std::uint16_t hitsMin  = 0;  //< 1 : deduced from database parameters
    float hitsDiffFraction = 0.9;
    //maximum range in sequence that read (pair) is expected to be in
    std::size_t insertSizeMax = 0;
    //limit number of reads per file (for quick tests)
    std::uint_least64_t queryLimit =
        std::numeric_limits<std::uint_least64_t>::max();

    //-----------------------------------------------------
    //analysis options
    //-----------------------------------------------------
    //test precision (ground truth must be available)
    bool testPrecision = false;
    bool testCoverage = false;
    bool testAlignment = false;
    //additional file with query -> ground truth mapping
    string sequ2taxonPreFile;

    //-----------------------------------------------------
    //query sampling scheme
    //-----------------------------------------------------
    int sketchlen = -1;  //< 0 : use value from database
    int winlen = -1;     //< 0 : use value from database
    int winstride = -1;  //< 0 : use value from database

    //-----------------------------------------------------
    // tuning parameters
    //-----------------------------------------------------
    int maxLocationsPerFeature = -1; //< 0 : use value from database
    int numThreads = std::thread::hardware_concurrency();
    int batchSize = 0; // automatic
    bool removeOverpopulatedFeatures = false;


    //-----------------------------------------------------
    //filenames
    //-----------------------------------------------------
    string dbfile;
    std::vector<string> infiles;
    string outfile;
};




/*****************************************************************************
 *
 * @brief command line args -> query options
 *
 *****************************************************************************/
query_options
get_query_options(const args_parser& args,
                  const query_options& defaults = query_options{})
{
    query_options opt;

    //files
    opt.dbfile = database_name(args);
    opt.infiles = sequence_filenames(args);

    opt.showDBproperties =  defaults.showDBproperties || args.contains("verbose");

    //pairing
    if(args.contains({"paired_files", "paired-files",
                      "pair_files", "pair-files", "pairfiles"}))
    {
        if(opt.infiles.size() > 1) {
            opt.pairing = pairing_mode::files;
            std::sort(opt.infiles.begin(), opt.infiles.end());
        }
    }
    else if(args.contains({"paired_sequences", "paired-sequences",
                           "pair_sequences", "pair-sequences",
                           "pair_sequ", "pair-sequ",
                           "pair_seq", "pair-seq",
                           "pairsequ", "pairseq", "paired"}) )
    {
        opt.pairing = pairing_mode::sequences;
    }

    //analysis options
    opt.testCoverage = defaults.testCoverage || args.contains("coverage");

    opt.testPrecision = defaults.testPrecision || opt.testCoverage ||
                        args.contains("precision");


    //classification ranks
    opt.lowestRank  = defaults.lowestRank;
    auto lowestRank = args.get<string>("lowest", "");
    if(!lowestRank.empty()) {
        auto r = taxonomy::rank_from_name(lowestRank);
        if(r < taxonomy::rank::root) opt.lowestRank = r;
    }

    opt.highestRank = defaults.highestRank;
    auto highestRank = args.get<string>("highest", "");
    if(!highestRank.empty()) {
        auto r = taxonomy::rank_from_name(highestRank);
        if(r < taxon_rank::root) opt.highestRank = r;
    }

    if(opt.lowestRank  > opt.highestRank) opt.lowestRank  = opt.highestRank;
    if(opt.highestRank < opt.lowestRank ) opt.highestRank = opt.lowestRank;

    opt.excludedRank = defaults.excludedRank;
    auto excludedRank = args.get<string>("exclude", "");
    if(!excludedRank.empty()) {
        auto r = taxonomy::rank_from_name(excludedRank);
        if(r < taxon_rank::root) opt.excludedRank = r;
    }

    //query sketching
    opt.sketchlen = args.get<int>("sketchlen", defaults.sketchlen);
    opt.winlen    = args.get<int>("winlen", defaults.winlen);
    opt.winstride = args.get<int>("winstride",
                                    opt.winlen > 0 ? opt.winlen : defaults.winstride);

    //query classification
    opt.hitsMin  = args.get<std::uint16_t>("hitmin", defaults.hitsMin);
    opt.hitsDiffFraction = args.get<float>("hitdiff", defaults.hitsDiffFraction);
    //interprest numbers > 1 as percentage
    if(opt.hitsDiffFraction > 1) opt.hitsDiffFraction *= 0.01;

    //alignment
    opt.showAlignment = defaults.testAlignment ||
                        args.contains({"showalign", "show-align", "show_align"});

    opt.testAlignment = opt.showAlignment ||
                          args.contains({"align", "alignment",
                                         "testalign", "test-align"});

    //output formatting
    opt.showLineage = defaults.showLineage || args.contains("lineage");

    opt.outfile = args.get<string>("out", defaults.outfile);

    opt.splitOutput = defaults.splitOutput ||
                      args.contains({"splitout","split-out"});

    if(opt.outfile.empty()) {
        opt.outfile = args.get<string>({"splitout","split-out"}, "");
    }

    opt.outSeparator = args.get<string>("separator", defaults.outSeparator);

    opt.showLocations = defaults.showLocations || args.contains("locations");

    opt.showTopHits = defaults.showTopHits ||
                      args.contains({"tophits", "top-hits"});

    opt.showAllHits = defaults.showAllHits ||
                      args.contains({"allhits", "all-hits"});


    opt.separateTaxaInfo = args.contains({"separate-cols",
                                          "separatecols", "separate_cols"});

    if(args.contains({"taxidsonly","taxids-only","taxids_only",
                      "taxidonly", "taxid-only", "taxid_only"}))
    {
        opt.showTaxaAs = taxon_print_mode::rank_id;
    }
    else if(args.contains({"taxids", "taxid"})) {
        opt.showTaxaAs = taxon_print_mode::rank_name_id;
    }
    else {
        opt.showTaxaAs = defaults.showTaxaAs;
    }

    opt.mapViewMode = defaults.mapViewMode;
    if(args.contains({"nomap","no-map","noshowmap","nomapping","nomappings"})) {
        opt.mapViewMode = map_view_mode::none;
    }
    else if(args.contains({"mapped-only", "mapped_only", "mappedonly"})) {
        opt.mapViewMode = map_view_mode::mapped_only;
    }
    else if(args.contains({"showallmap","show-all-map","show-all-mappings"})) {
        opt.mapViewMode = map_view_mode::all;
    }

    //showing hits changes the mapping mode!
    if(opt.mapViewMode == map_view_mode::none && opt.showTopHits) {
        opt.mapViewMode = map_view_mode::mapped_only;
    }
    else if(opt.showAllHits) opt.mapViewMode = map_view_mode::all;

    opt.showGroundTruth = defaults.showGroundTruth ||
                          args.contains({"ground-truth", "ground_truth",
                                         "groundtruth"});

    opt.insertSizeMax = args.get<std::size_t>({"insertsize", "insert-size"},
                                                defaults.insertSizeMax);

    opt.queryLimit = args.get<std::uint_least64_t>({"query-limit"},
                                                   defaults.queryLimit);

    //database tuning parameters
    opt.maxLocationsPerFeature = args.get<int>({"max_locations_per_feature",
                                                "max-locations-per-feature"},
                                                defaults.maxLocationsPerFeature);

    opt.removeOverpopulatedFeatures = defaults.removeOverpopulatedFeatures ||
        args.contains({"remove-overpopulated-features",
                       "remove_overpopulated_features" });

    opt.numThreads = args.get<int>("threads", defaults.numThreads);
    opt.batchSize = args.get<int>("batch-size", defaults.batchSize);

    opt.showErrors = defaults.showErrors && !args.contains({"-noerr","-noerrors"});

    return opt;
}



/*****************************************************************************
 *
 *
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
 *
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



/*****************************************************************************
 *
 * @brief classifies single sequences from a single source
 *
 *****************************************************************************/
void classify(parallel_queue& queue,
    const database& db, const query_options& opt,
    sequence_reader& reader, std::ostream& os, classification_statistics& stats)
{
    const auto load = 32 * queue.concurrency();
    const auto batchSize = opt.batchSize > 0 ? opt.batchSize
                                             : 4096 * queue.concurrency();
    std::mutex mtx;
    std::atomic<std::uint64_t> queryLimit(opt.queryLimit);

    while(reader.has_next() && queryLimit > 0) {
        if(queue.unsafe_waiting() < load) {
            queue.enqueue([&]{
                std::ostringstream obuf;

                for(std::size_t i = 0;
                    i < batchSize && queryLimit > 0 && reader.has_next();
                    ++i, --queryLimit)
                {
                    auto query = reader.next();
                    if(!query.header.empty()) {
                        auto matches = db.matches(query.data);

                        process_database_answer(db, opt,
                            query.header, query.data, sequence{},
                            std::move(matches), stats, obuf);
                    }
                }
                //flush output buffer
                std::lock_guard<std::mutex> lock(mtx);
                os << obuf.str();
            }); //enqueue
        }
    }
    //wait for all enqueued tasks to finish
    queue.wait();
}



/*****************************************************************************
 *
 * @brief classifies pairs of sequences from 2 sources
 *
 *****************************************************************************/
void classify_pairs(parallel_queue& queue,
                    const database& db, const query_options& opt,
                    sequence_reader& reader1, sequence_reader& reader2,
                    std::ostream& os, classification_statistics& stats)
{
    const auto load = 32 * queue.concurrency();
    const auto batchSize = opt.batchSize > 0 ? opt.batchSize
                                             : 4096 * queue.concurrency();
    std::mutex mtx1;
    std::mutex mtx2;
    std::atomic<std::uint64_t> queryLimit(opt.queryLimit);

    while(reader1.has_next() && reader2.has_next() && queryLimit > 0) {
        if(queue.unsafe_waiting() < load) {
            queue.enqueue([&]{
                std::ostringstream obuf;

                for(std::size_t i = 0; i < batchSize && queryLimit > 0;
                    ++i, --queryLimit)
                {
                    std::unique_lock<std::mutex> lock1(mtx1);
                    if(!reader1.has_next() || !reader2.has_next()) break;
                    //make sure that both queries are read in sequence
                    auto query1 = reader1.next();
                    auto query2 = reader2.next();
                    lock1.unlock();

                    if(!query1.header.empty()) {
                        auto matches = db.matches(query1.data);
                        db.accumulate_matches(query2.data, matches);

                        process_database_answer(db, opt,
                            query1.header, query1.data, query2.data,
                            std::move(matches), stats, obuf);
                    }
                }
                //flush output buffer
                std::lock_guard<std::mutex> lock(mtx2);
                os << obuf.str();
            }); //enqueue
        }
    }
    //wait for all enqueued tasks to finish
    queue.wait();
}



/*****************************************************************************
 *
 * @brief classify sequences per input file
 *
 *****************************************************************************/
void classify_per_file(parallel_queue& queue,
                       const database& db, const query_options& opt,
                       const std::vector<string>& infilenames,
                       std::ostream& os, classification_statistics& stats)
{
    for(size_t i = 0; i < infilenames.size(); ++i) {
        const auto& fname = infilenames[i];
        if(opt.mapViewMode != map_view_mode::none) {
            os << opt.comment << fname << '\n';
        }
        else if(!opt.outfile.empty() && infilenames.size() > 1) {
            show_progress_indicator(i / float(infilenames.size()));
        }
        try {
            auto reader = make_sequence_reader(fname);
            if(opt.pairing == pairing_mode::sequences) {
                classify_pairs(queue, db, opt, *reader, *reader, os, stats);
            } else {
                classify(queue, db, opt, *reader, os, stats);
            }
        }
        catch(std::exception& e) {
            if(opt.mapViewMode != map_view_mode::none) {
                os << "FAIL: " << e.what() << endl;
            }
        }
    }
}



/*****************************************************************************
 *
 * @brief classifies paired-end sequences in pairs of files
 *
 *****************************************************************************/
void classify_on_file_pairs(parallel_queue& queue,
                            const database& db, const query_options& opt,
                            const std::vector<string>& infilenames,
                            std::ostream& os, classification_statistics& stats)
{
    //pair up reads from two consecutive files in the list
    for(size_t i = 0; i < infilenames.size(); i += 2) {
        const auto& fname1 = infilenames[i];
        const auto& fname2 = infilenames[i+1];
        if(opt.mapViewMode != map_view_mode::none) {
            os << opt.comment << fname1 << " + " << fname2 << '\n';
        }
        else if(!opt.outfile.empty() && infilenames.size() > 1) {
            show_progress_indicator(i / float(infilenames.size()));
        }
        try {
            auto reader1 = make_sequence_reader(fname1);
            auto reader2 = make_sequence_reader(fname2);
            classify_pairs(queue, db, opt, *reader1, *reader2, os, stats);
        }
        catch(std::exception& e) {
            if(opt.mapViewMode != map_view_mode::none) {
                os << "FAIL: " << e.what() << endl;
            }
        }
    }
}



/*****************************************************************************
 *
 * @brief prints classification parameters
 *
 *****************************************************************************/
void show_classification_parameters(const query_options& opt, std::ostream& os)
{
    const auto& comment = opt.comment;
    if(opt.mapViewMode != map_view_mode::none) {
        os << comment << "Reporting per-read mappings (non-mapping lines start with '"
           << comment << "').\n";

        os << comment
           << "Output will be constrained to ranks from '"
           << taxonomy::rank_name(opt.lowestRank) << "' to '"
           << taxonomy::rank_name(opt.highestRank) << "'.\n";

        if(opt.showLineage) {
            os << comment
               << "The complete lineage will be reported "
               << "starting with the lowest match.\n";
        }
        else {
            os << comment
               << "Only the lowest matching rank will be reported.\n";
        }
    }
    else {
        os << comment << "Per-Read mappings will not be shown.\n";
    }
    if(opt.excludedRank != taxon_rank::none) {
        os << comment << "Clade Exclusion on Rank: "
           << taxonomy::rank_name(opt.excludedRank);
    }
    if(opt.pairing == pairing_mode::files) {
        os << comment << "File based paired-end mode:\n"
           << comment << "  Reads from two consecutive files will be interleaved.\n"
           << comment << "  Max insert size considered " << opt.insertSizeMax << ".\n";
    }
    else if(opt.pairing == pairing_mode::sequences) {
        os << comment << "Per file paired-end mode:\n"
           << comment << "  Reads from two consecutive sequences in each file will be paired up.\n"
           << comment << "  Max insert size considered " << opt.insertSizeMax << ".\n";
    }

    if(opt.testAlignment) {
        os << comment << "Query sequences will be aligned to best candidate target => SLOW!\n";
    }

    os << comment << "Using " << opt.numThreads << " threads\n";

}



/*****************************************************************************
 *
 * @brief prints classification parameters & results
 *        decides how to handle paired sequences
 *
 *****************************************************************************/
void classify_sequences(const database& db, const query_options& opt,
                        const std::vector<string>& infilenames,
                        std::ostream& os)
{
    show_classification_parameters(opt, os);

    if(opt.outfile.empty()) os.flush();

    parallel_queue queue(opt.numThreads);

    timer time;
    time.start();

    classification_statistics stats;
    if(opt.pairing == pairing_mode::files) {
        classify_on_file_pairs(queue, db, opt, infilenames, os, stats);
    }
    else {
        classify_per_file(queue, db, opt, infilenames, os, stats);
    }
    if(!opt.outfile.empty() || opt.mapViewMode == map_view_mode::none) {
        clear_current_line();
    }

    time.stop();

    //show results
    const auto numQueries = (opt.pairing == pairing_mode::none)
                            ? stats.total() : 2 * stats.total();

    const auto speed = numQueries / time.minutes();
    const auto& comment = opt.comment;
    os << comment << "queries: " << numQueries << '\n'
       << comment << "time:    " << time.milliseconds() << " ms\n"
       << comment << "speed:   " << speed << " queries/min\n";

    if(stats.total() > 0) {
        show_classification_statistics(os, stats, comment);
    } else {
        os << comment << "No valid query sequences found." << endl;
    }

    if(opt.outfile.empty()) os.flush();
}



/*****************************************************************************
 *
 * @brief runs classification on several input files and writes the results
 *        to an output file or stdout
 *
 *****************************************************************************/
void process_input_files(const database& db, const query_options& opt,
                         const std::vector<string>& infilenames,
                         const string& outfilename)
{

    if(outfilename.empty()) {
        classify_sequences(db, opt, infilenames, cout);
    }
    else {
        std::ofstream os {outfilename, std::ios::out};

        if(os.good()) {
            cout << "Output will be redirected to file: " << outfilename << endl;

            classify_sequences(db, opt, infilenames, os);
        }
        else {
            throw file_write_error{
                "Could not write to output file " + outfilename};
        }
    }
}



/*****************************************************************************
 *
 * @brief runs classification on all input files
 *
 *****************************************************************************/
void process_input_files(const database& db, const query_options& opt)
{
    if(opt.infiles.empty()) {
        cerr << "No input filenames provided!" << endl;
        return;
    }
    else {
        bool anyreadable = false;
        for(const auto& f : opt.infiles) {
            if(file_readable(f)) { anyreadable = true; break; }
        }
        if(!anyreadable) {
            throw std::runtime_error{
                        "None of the query sequence files could be opened"};
        }
    }

    //process files / file pairs separately
    if(opt.splitOutput) {
        string outfile;
        //process each input file pair separately
        if(opt.pairing == pairing_mode::files && opt.infiles.size() > 1) {
            for(std::size_t i = 0; i < opt.infiles.size(); i += 2) {
                const auto& f1 = opt.infiles[i];
                const auto& f2 = opt.infiles[i+1];
                if(!opt.outfile.empty()) {
                    outfile = opt.outfile + "_" + extract_filename(f1)
                                            + "_" + extract_filename(f2)
                                            + ".txt";
                }
                process_input_files(db, opt,
                                    std::vector<string>{f1,f2}, outfile);
            }
        }
        //process each input file separately
        else {
            for(const auto& f : opt.infiles) {
                if(!opt.outfile.empty()) {
                    outfile = opt.outfile + "_" + extract_filename(f) + ".txt";
                }
                process_input_files(db, opt,
                                    std::vector<string>{f}, outfile);
            }
        }
    }
    //process all input files at once
    else {
        process_input_files(db, opt, opt.infiles, opt.outfile);
    }
}



/*****************************************************************************
 *
 * @brief modify database content and sketching scheme according to
 *        command line options
 *
 *****************************************************************************/
void configure_database_according_to_query_options(
    database& db, const query_options& opt)
{
    if(opt.removeOverpopulatedFeatures) {
        auto old = db.feature_count();

        auto maxlpf = opt.maxLocationsPerFeature - 1;
        if(maxlpf < 0 || maxlpf >= database::max_supported_locations_per_feature())
            maxlpf = database::max_supported_locations_per_feature() - 1;

        maxlpf = std::min(maxlpf, db.max_locations_per_feature() - 1);
        if(maxlpf > 0) { //always keep buckets with size 1
            cout << "\nRemoving features with more than "
                 << maxlpf << " locations... " << std::flush;

            auto rem = db.remove_features_with_more_locations_than(maxlpf);

            cout << rem << " of " << old << " removed." << endl;
        }
        //in case new max is less than the database setting
        db.max_locations_per_feature(opt.maxLocationsPerFeature);
    }
    else if(opt.maxLocationsPerFeature > 1) {
        db.max_locations_per_feature(opt.maxLocationsPerFeature);
        cout << "max locations per feature set to "
             << opt.maxLocationsPerFeature << endl;
    }

    //use a different sketching scheme for querying?
    if(opt.winlen > 0)    db.query_window_size(opt.winlen);
    if(opt.winstride > 0) db.query_window_stride(opt.winstride);

    if(opt.sketchlen > 0) {
        auto s = db.query_sketcher();
        s.sketch_size(opt.sketchlen);
        db.query_sketcher(std::move(s));
    }
}



/*****************************************************************************
 *
 * @brief sets up some query options according to database parameters
 *        or command line options
 *
 *****************************************************************************/
void configure_query_options_according_to_database(
    query_options& opt, const database& db)
{
    //deduce hit threshold from database?
    if(opt.hitsMin < 1) {
        auto sks = db.target_sketcher().sketch_size();
        if(sks >= 6) {
            opt.hitsMin = static_cast<int>(sks / 3.0);
        } else if (sks >= 4) {
            opt.hitsMin = 2;
        } else {
            opt.hitsMin = 1;
        }
    }
}



/*****************************************************************************
 *
 * @brief primitive REPL mode for repeated querying using the same database
 *
 *****************************************************************************/
void run_interactive_query_mode(const database& db, const query_options& initOpt)
{
    while(true) {
        cout << "$> " << std::flush;

        string input;
        std::getline(std::cin, input);
        if(input.empty() || input.find(":q") == 0) {
            cout << "Terminate." << endl;
            return;
        }
        else if(input.find("#") == 0)
        {
            //comment line, do nothing
        }
        else {
            //tokenize input into whitespace-separated words and build args list
            std::vector<string> args {"query", initOpt.dbfile};
            std::istringstream iss(input);
            while(iss >> input) { args.push_back(input); }

            //read command line options (use initial ones as defaults)
            auto opt = get_query_options(args_parser{std::move(args)}, initOpt);
            configure_query_options_according_to_database(opt, db);

            try {
                process_input_files(db, opt);
            }
            catch(std::exception& e) {
                if(initOpt.showErrors) cerr << e.what() << endl;
            }
        }
    }
}



/*****************************************************************************
 *
 * @brief    run query reads against pre-built database
 *           entry point for query mode
 *
 * @details  note that precision (positive predictive value) and
 *           clade exclusion testing is much slower and only intended
 *           for classification performance evaluation
 *
 *****************************************************************************/
void main_mode_query(const args_parser& args)
{
    auto opt = get_query_options(args);

    auto db = make_database<database>(opt.dbfile);

    configure_database_according_to_query_options(db, opt);

    if(opt.showDBproperties) {
        print_static_properties(db);
        print_content_properties(db);
    }

    if(!opt.infiles.empty()) {
        cerr << "Classifying query sequences." << endl;

        configure_query_options_according_to_database(opt, db);
        process_input_files(db, opt);
    }
    else {
        cout << "No input files provided.\n"
            "Running in interactive mode:\n"
            " - Enter input file name(s) and command line options and press return.\n"
            " - The initially given command line options will be used as defaults.\n"
            " - All command line options that would modify the database are ignored.\n"
            " - Each line will be processed separately.\n"
            " - Lines starting with '#' will be ignored.\n"
            " - Enter an empty line or press Ctrl-D to quit MetaCache.\n"
            << endl;

        run_interactive_query_mode(db, opt);
    }
}


} // namespace mc

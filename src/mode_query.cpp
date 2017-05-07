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
    std::string comment = "# ";
    //separates individual mapping fields
    std::string outSeparator = "\t|\t";
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
    taxon_print_mode showTaxaAs = taxon_print_mode::name_only;
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
    float hitsDiff = 0.5;
    //maximum range in sequence that read (pair) is expected to be in
    std::size_t insertSizeMax = 0;

    //-----------------------------------------------------
    //analysis options
    //-----------------------------------------------------
    //test precision (ground truth must be available)
    bool testPrecision = false;
    bool testCoverage = false;
    bool testAlignment = false;
    //additional file with query -> ground truth mapping
    std::string sequ2taxonPreFile;

    //-----------------------------------------------------
    //query sampling scheme
    //-----------------------------------------------------
    int sketchlen = -1;  //< 0 : use value from database
    int winlen = -1;     //< 0 : use value from database
    int winstride = -1;  //< 0 : use value from database

    //-----------------------------------------------------
    // tuning parameters
    //-----------------------------------------------------
    float maxLoadFactor = -1;        //< 0 : use value from database
    int maxTargetsPerSketchVal = -1; //< 0 : use value from database
    int numThreads = 1;

    //-----------------------------------------------------
    //filenames
    //-----------------------------------------------------
    std::string dbfile;
    std::vector<std::string> infiles;
    std::string outfile;
};




/*****************************************************************************
 *
 * @brief command line args -> query options
 *
 *****************************************************************************/
query_options
get_query_options(const args_parser& args)
{
    const query_options defaults;

    query_options param;

    //files
    param.dbfile = database_name(args);

    param.infiles = sequence_filenames(args);
    if(param.infiles.empty()) {
        throw std::runtime_error{"No query sequence filenames provided"};
    }
    else {
        bool anyreadable = false;
        for(const auto& f : param.infiles) {
            if(file_readable(f)) { anyreadable = true; break; }
        }
        if(!anyreadable) {
            throw std::runtime_error{
                        "None of the query sequence files could be opened"};
        }
    }

    param.showDBproperties = args.contains("verbose");

    //pairing
    if(args.contains({"paired_files", "paired-files",
                      "pair_files", "pair-files", "pairfiles"}))
    {
        if(param.infiles.size() > 1) {
            param.pairing = pairing_mode::files;
            std::sort(param.infiles.begin(), param.infiles.end());
        }
    }
    else if(args.contains({"paired_sequences", "paired-sequences",
                           "pair_sequences", "pair-sequences",
                           "pair_sequ", "pair-sequ",
                           "pair_seq", "pair-seq",
                           "pairsequ", "pairseq", "paired"}) )
    {
        param.pairing = pairing_mode::sequences;
    }

    //analysis options
    param.testCoverage = args.contains("coverage");
    param.testPrecision = param.testCoverage || args.contains("precision");

    //classification ranks
    auto lowestRank = args.get<std::string>("lowest", "");
    if(!lowestRank.empty()) {
        auto r = taxonomy::rank_from_name(lowestRank);
        if(r < taxonomy::rank::root) {
            param.lowestRank = r;
        }
    }

    auto highestRank = args.get<std::string>("highest", "");
    if(!highestRank.empty()) {
        auto r = taxonomy::rank_from_name(highestRank);
        if(r < taxon_rank::root) param.highestRank = r;
    }

    if(param.lowestRank > param.highestRank)
        param.lowestRank = param.highestRank;

    if(param.highestRank < param.lowestRank )
        param.highestRank = param.lowestRank;

    auto excludedRank = args.get<std::string>("exclude", "");
    if(!excludedRank.empty()) {
        auto r = taxonomy::rank_from_name(excludedRank);
        if(r < taxon_rank::root) param.excludedRank = r;
    }

    //query sketching
    param.sketchlen = args.get<int>("sketchlen", defaults.sketchlen);
    param.winlen    = args.get<int>("winlen", defaults.winlen);
    param.winstride = args.get<int>("winstride",
                                    param.winlen > 0 ? param.winlen : defaults.winstride);

    //query classification
    param.hitsMin  = args.get<std::uint16_t>("hitmin", defaults.hitsMin);
    param.hitsDiff = args.get<float>("hitdiff", defaults.hitsDiff);

    //alignment
    param.showAlignment = args.contains("showalign");

    param.testAlignment = param.showAlignment ||
                          args.contains({"align", "alignment",
                                         "testalign", "test-align"});

    //output formatting
    param.showLineage = args.contains("lineage");

    param.outfile = args.get<std::string>("out", "");

    if(param.outfile.empty()) {
        param.outfile = args.get<std::string>({"splitout","split-out"}, "");
        param.splitOutput = true;
    }
    else {
        param.splitOutput = args.contains({"splitout","split-out"});
    }

    param.outSeparator = args.get<std::string>("separator", "\t|\t");

    param.showLocations = args.contains("locations");

    param.showTopHits = args.contains({"tophits", "top-hits"});
    param.showAllHits = args.contains({"allhits", "all-hits"});

    if(args.contains({"taxidsonly","taxids-only","taxids_only",
                      "taxidonly", "taxid-only", "taxid_only"}))
    {
        param.showTaxaAs = taxon_print_mode::id_only;
    }
    else if(args.contains("taxids") || args.contains("taxid")) {
        param.showTaxaAs = taxon_print_mode::id_name;
    }
    else {
        param.showTaxaAs = taxon_print_mode::name_only;
    }

    if(args.contains({"nomap","no-map","noshowmap","nomapping","nomappings"})) {
        param.mapViewMode = map_view_mode::none;
    }
    else if(args.contains({"mapped-only", "mapped_only", "mappedonly"})) {
        param.mapViewMode = map_view_mode::mapped_only;
    }

    //showing hits changes the mapping mode!
    if(param.mapViewMode == map_view_mode::none && param.showTopHits) {
        param.mapViewMode = map_view_mode::mapped_only;
    }
    else if(param.showAllHits) param.mapViewMode = map_view_mode::all;

    param.showGroundTruth = args.contains({"ground-truth", "ground_truth",
                                           "groundtruth"});

    param.insertSizeMax = args.get<std::size_t>({"insertsize", "insert-size"},
                                                defaults.insertSizeMax);

    //database tuning parameters
    param.maxLoadFactor = args.get<float>({"max-load-fac", "max_load_fac"},
                                          defaults.maxLoadFactor);

    param.maxTargetsPerSketchVal = args.get<int>({"max_locations_per_feature"
                                                  "max-locations-per-feature"},
                                                defaults.maxTargetsPerSketchVal);

    param.numThreads = args.get<int>("threads",
                                     std::thread::hardware_concurrency());

    param.showErrors = !args.contains({"-noerr","-noerrors"});

    return param;
}



/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
template<class Sequence, class Index>
inline auto
make_view_from_window_range(
    const Sequence& s, const index_range<Index>& winRange, int stride = 1)
    -> decltype(make_view(s.begin(),s.end()))
{
    auto end = s.begin() + (stride * winRange.end);
    if(end > s.end()) end = s.end();

    return make_view(s.begin() + (stride * winRange.beg), end);
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
ground_truth(const database& db, const std::string& header)
{
    //try to extract query id and find the corresponding target in database
    const taxon* tax = nullptr;
    tax = db.taxon_with_name(extract_ncbi_accession_version_number(header));
    if(tax) return db.next_ranked(tax);

    tax = db.taxon_with_name(extract_ncbi_accession_number(header));
    if(tax) return db.next_ranked(tax);

    //try to extract id from header
    tax = db.taxon_with_id(extract_taxon_id(header));
    if(tax) return db.next_ranked(tax);

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
    const taxon* excl = db.ranks(tax)[int(rank)];

    auto maskedHits = hits;
    for(const auto& hit : hits) {
        const taxon* tgt = db.taxon_of_target(hit.first.tgt);
        if(tgt) {
            auto t = db.ranks(*tgt)[int(rank)];
            if(t != excl) maskedHits.insert(hit);
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
    taxon_rank lowestRank = taxon_rank::subSpecies,
    taxon_rank highestRank = taxon_rank::Domain)
{
    if(cand.count() < 3) {
        const taxon* tax = db.ranked_lca(db.taxon_of_target(cand.target(0)),
                                         db.taxon_of_target(cand.target(1)) );

        //classify if rank is below or at the highest rank of interest
        if(tax && tax->rank() <= highestRank) return tax;
    }
    else {
        if(lowestRank == taxon_rank::Sequence) {
            lowestRank = taxonomy::next_main_rank(lowestRank);
        }

        std::unordered_map<const taxon*,int> scores;
        scores.rehash(2*cand.count());

        for(auto r = lowestRank; r <= highestRank; r = taxonomy::next_main_rank(r)) {
            //hash-count taxon id occurrences on rank 'r'
            int totalscore = 0;
            for(int i = 0, n = cand.count(); i < n; ++i) {
                //use target id instead of taxon if at sequence level
                const taxon* tax = db.ranks(db.taxon_of_target(cand.target(i)))[int(r)];
                auto score = cand.hits(i);
                totalscore += score;
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
            if(toptax && topscore >= (totalscore * trustedMajority)) {
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
sequence_classification(
    const database& db, const query_options& param,
    const classification_candidates& cand)
{
    //sum of top-2 hits < threshold => considered not classifiable
    if((cand.hits(0) + cand.hits(1)) < param.hitsMin) return nullptr;

    //either top 2 are the same sequences with at least 'hitsMin' many hits
    //(checked before) or hit difference between these top 2 is above threshhold
    if( (cand.target(0) == cand.target(1))
        || (cand.hits(0) - cand.hits(1)) >= param.hitsMin)
    {
        //return top candidate
        return db.taxon_of_target(cand.target(0));
    }

    return lowest_common_taxon(db, cand, param.hitsDiff,
                               param.lowestRank, param.highestRank);
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
        os << "--";
    }
    else {
        auto rmin = opt.lowestRank < tax->rank() ? tax->rank() : opt.lowestRank;
        auto rmax = opt.showLineage ? opt.highestRank : rmin;

        show_ranks(os, db.ranks(tax), opt.showTaxaAs, rmin, rmax);
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
 * @brief process database hit results: classify & print mapping
 *
 *****************************************************************************/
void process_database_answer(
     const database& db, const query_options& opt,
     const std::string& header,
     const sequence& query1, const sequence& query2,
     matches_per_location&& hits, std::ostream& os, classification_statistics& stats)
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

    classification_candidates tophits {hits, numWindows};
    const taxon* cls = sequence_classification(db, opt, tophits );

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
        if(l != std::string::npos) {
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
                    opt.showLineage ? opt.highestRank : groundTruth->rank());
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
        //check which of the top sequences has a better alignment
        const taxon* tgtTax = db.taxon_of_target(tophits.target(0));
        if(tgtTax) {
            const auto& src = tgtTax->source();
            try {
                //load candidate file and forward to sequence
                auto reader = make_sequence_reader(src.filename);
                reader->skip(src.index);

                if(reader->has_next()) {
                    auto tgtSequ = reader->next().data;
                    auto subject = make_view_from_window_range(tgtSequ,
                                        tophits.window(0),
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
                           << " in range [" << (w * tophits.window(0).beg)
                           << ',' << (w * tophits.window(0).end) << "]\n"
                           << opt.comment << "  query  " << align.query << '\n'
                           << opt.comment << "  target " << align.subject;
                    }
                }
            } catch(std::exception& e) {
                if(opt.showErrors) std::cerr << e.what() << '\n';
            }
        }
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
    std::mutex mtx;

    const auto load = 32 * queue.concurrency();
    const auto batchSize = 4096 * queue.concurrency();

//    std::vector<std::ostringstream> outbufs(load);

    while(reader.has_next()) {
        if(queue.unsafe_waiting() < load) {
            queue.enqueue([&]{
                std::ostringstream obuf;

                for(std::size_t i = 0; i < batchSize; ++i) {
                    auto query = reader.next();
                    if(!query.header.empty()) {
                        auto matches = db.matches(query.data);

                        process_database_answer(db, opt,
                            query.header, query.data, sequence{},
                            std::move(matches), obuf, stats);
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
                    const database& db, const query_options& param,
                    sequence_reader& reader1, sequence_reader& reader2,
                    std::ostream& os, classification_statistics& stats)
{
    const auto load = 32 * queue.concurrency();
    const auto batchSize = 4096 * queue.concurrency();

    std::mutex mtx1;
    std::mutex mtx2;

    while(reader1.has_next() && reader2.has_next()) {
        if(queue.unsafe_waiting() < load) {
            queue.enqueue([&]{
                std::ostringstream obuf;

                for(std::size_t i = 0; i < batchSize; ++i) {
                    std::unique_lock<std::mutex> lock1(mtx1);
                    //make sure that both queries are read in sequence
                    auto query1 = reader1.next();
                    auto query2 = reader2.next();
                    lock1.unlock();

                    if(!query1.header.empty()) {
                        auto matches = db.matches(query1.data);
                        db.accumulate_matches(query2.data, matches);

                        process_database_answer(db, param,
                            query1.header, query1.data, query2.data,
                            std::move(matches), obuf, stats);
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
                       const database& db, const query_options& param,
                       const std::vector<std::string>& infilenames,
                       std::ostream& os, classification_statistics& stats)
{
    for(size_t i = 0; i < infilenames.size(); ++i) {
        const auto& fname = infilenames[i];
        if(param.mapViewMode != map_view_mode::none) {
            os << param.comment << fname << '\n';
        }
        else if(!param.outfile.empty() && infilenames.size() > 1) {
            show_progress_indicator(i / float(infilenames.size()));
        }
        try {
            auto reader = make_sequence_reader(fname);
            if(param.pairing == pairing_mode::sequences) {
                classify_pairs(queue, db, param, *reader, *reader, os, stats);
            } else {
                classify(queue, db, param, *reader, os, stats);
            }
        }
        catch(std::exception& e) {
            if(param.mapViewMode != map_view_mode::none) {
                os << "FAIL: " << e.what() << std::endl;
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
                            const database& db, const query_options& param,
                            const std::vector<std::string>& infilenames,
                            std::ostream& os, classification_statistics& stats)
{
    //pair up reads from two consecutive files in the list
    for(size_t i = 0; i < infilenames.size(); i += 2) {
        const auto& fname1 = infilenames[i];
        const auto& fname2 = infilenames[i+1];
        if(param.mapViewMode != map_view_mode::none) {
            os << param.comment << fname1 << " + " << fname2 << '\n';
        }
        else if(!param.outfile.empty() && infilenames.size() > 1) {
            show_progress_indicator(i / float(infilenames.size()));
        }
        try {
            auto reader1 = make_sequence_reader(fname1);
            auto reader2 = make_sequence_reader(fname2);
            classify_pairs(queue, db, param, *reader1, *reader2, os, stats);
        }
        catch(std::exception& e) {
            if(param.mapViewMode != map_view_mode::none) {
                os << "FAIL: " << e.what() << std::endl;
            }
        }
    }
}



/*****************************************************************************
 *
 * @brief prints classification parameters & results
 *        decides how to handle paired sequences
 *
 *****************************************************************************/
void classify_sequences(const database& db, const query_options& param,
                        const std::vector<std::string>& infilenames,
                        std::ostream& os)
{
    const auto& comment = param.comment;

    if(param.mapViewMode != map_view_mode::none) {
        os << comment << "Reporting per-read mappings (non-mapping lines start with '"
           << comment << "').\n";

        os << comment
           << "Output will be constrained to ranks from '"
           << taxonomy::rank_name(param.lowestRank) << "' to '"
           << taxonomy::rank_name(param.highestRank) << "'.\n";

        if(param.showLineage) {
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
    if(param.excludedRank != taxon_rank::none) {
        os << comment << "Clade Exclusion on Rank: "
           << taxonomy::rank_name(param.excludedRank);
    }
    if(param.pairing == pairing_mode::files) {
        os << comment << "File based paired-end mode:\n"
           << comment << "  Reads from two consecutive files will be interleaved.\n"
           << comment << "  Max insert size considered " << param.insertSizeMax << ".\n";
    }
    else if(param.pairing == pairing_mode::sequences) {
        os << comment << "Per file paired-end mode:\n"
           << comment << "  Reads from two consecutive sequences in each file will be paired up.\n"
           << comment << "  Max insert size considered " << param.insertSizeMax << ".\n";
    }

    if(param.testAlignment) {
        os << comment << "Query sequences will be aligned to best candidate target => SLOW!\n";
    }

    os << comment << "Using " << param.numThreads << " threads\n";
    if(param.outfile.empty()) os.flush();

    parallel_queue queue(param.numThreads);

    timer time;
    time.start();

    classification_statistics stats;
    if(param.pairing == pairing_mode::files) {
        classify_on_file_pairs(queue, db, param, infilenames, os, stats);
    }
    else {
        classify_per_file(queue, db, param, infilenames, os, stats);
    }
    if(!param.outfile.empty() || param.mapViewMode == map_view_mode::none) {
        clear_current_line();
    }

    time.stop();

    //show results
    const auto numQueries = (param.pairing == pairing_mode::none)
                            ? stats.total() : 2 * stats.total();

    const auto speed = numQueries / time.minutes();
    os << comment << "queries: " << numQueries << '\n'
       << comment << "time:    " << time.milliseconds() << " ms\n"
       << comment << "speed:   " << speed << " queries/min\n";

    if(stats.total() > 0) {
        show_classification_statistics(os, stats, comment);
    } else {
        os << comment << "No valid query sequences found." << std::endl;
    }

    if(param.outfile.empty()) os.flush();
}



/*****************************************************************************
 *
 * @brief descides where to put the classification results
 *
 *****************************************************************************/
void process_input_files(const database& db, const query_options& param,
                         const std::vector<std::string>& infilenames,
                         const std::string& outfilename)
{

    if(outfilename.empty()) {
        classify_sequences(db, param, infilenames, std::cout);
    }
    else {
        std::ofstream os {outfilename, std::ios::out};

        if(os.good()) {
            std::cout << "Output will be redirected to file: "
                      << outfilename << std::endl;

            classify_sequences(db, param, infilenames, os);
        }
        else {
            throw file_write_error{
                "Could not write to output file " + outfilename};
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
    std::cout << "Classifying query sequences." << std::endl;

    auto param = get_query_options(args);

    auto db = make_database<database>(param.dbfile);

    //configure database
    if(param.maxLoadFactor > 0) {
        db.max_load_factor(param.maxLoadFactor);
    }
    if(param.maxTargetsPerSketchVal > 1) {
        db.remove_features_with_more_locations_than(param.maxTargetsPerSketchVal);
    }

    //deduced query parameters
    if(param.hitsMin < 1) {
        auto sks = db.target_sketcher().sketch_size();
        if(sks >= 6) {
            param.hitsMin = static_cast<int>(sks / 3.0);
        } else if (sks >= 4) {
            param.hitsMin = 2;
        } else {
            param.hitsMin = 1;
        }
    }

    //use a different sketching scheme for querying?
    if(param.winlen > 0)    db.query_window_size(param.winlen);
    if(param.winstride > 0) db.query_window_stride(param.winstride);

    if(param.sketchlen > 0) {
        auto s = db.query_sketcher();
        s.sketch_size(param.sketchlen);
        db.query_sketcher(std::move(s));
    }

    if(param.showDBproperties) {
        print_properties(db);
        std::cout << '\n';
    }

    //process files / file pairs separately
    if(param.splitOutput) {
        std::string outfile;
        //process each input file pair separately
        if(param.pairing == pairing_mode::files && param.infiles.size() > 1) {
            for(std::size_t i = 0; i < param.infiles.size(); i += 2) {
                const auto& f1 = param.infiles[i];
                const auto& f2 = param.infiles[i+1];
                if(!param.outfile.empty()) {
                    outfile = param.outfile + "_" + extract_filename(f1)
                                            + "_" + extract_filename(f2)
                                            + ".txt";
                }
                process_input_files(db, param,
                                    std::vector<std::string>{f1,f2}, outfile);
            }
        }
        //process each input file separately
        else {
            for(const auto& f : param.infiles) {
                if(!param.outfile.empty()) {
                    outfile = param.outfile + "_" + extract_filename(f) + ".txt";
                }
                process_input_files(db, param,
                                    std::vector<std::string>{f}, outfile);
            }
        }
    }
    //process all input files at once
    else {
        process_input_files(db, param, param.infiles, param.outfile);
    }
}


} // namespace mc

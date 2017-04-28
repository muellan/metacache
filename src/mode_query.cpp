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

#include "modes_common.h"


namespace mc {

/*****************************************************************************
 * @brief type aliases
 *****************************************************************************/
using parallel_queue = parallel_function_queue;


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
    //prefix for each non-mapping line
    std::string comment = "# ";
    //separates individual mapping fields
    std::string outSeparator = "\t|\t";

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
    bool weightHitsWithWindows = false;

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

    param.weightHitsWithWindows = args.contains("hitweight");

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

    return param;
}



/*****************************************************************************
 *
 * @brief  print target.window hit statistics from database
 *
 *****************************************************************************/
template<int n>
void show_candidate_ranges(std::ostream& os,
    const database& db,
    const matches_in_contiguous_window_range_top<n,target_id>& cand)
{
    const auto w = db.target_window_stride();

    for(int i = 0; i < n; ++i) {
        os << '[' << (w * cand.window(i).beg)
           << ',' << (w * cand.window(i).end) << "] ";
    }
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
 * @brief aligns one query (pair) against the top 2 candidate windows
 *
 * @param  query1  read
 * @param  query2  paired read
 *
 * @return index of best candidate or -1 if ambiguous or unsatisfactory
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
 * @pre cand.count() >= 2
 *
 *****************************************************************************/
template<int maxn>
const taxon*
lowest_common_taxon(
    const database& db,
    const matches_in_contiguous_window_range_top<maxn,target_id>& cand,
    float trustedMajority,
    taxon_rank lowestRank = taxon_rank::subSpecies,
    taxon_rank highestRank = taxon_rank::Domain)
{
    if(maxn < 3 || cand.count() < 3) {
        auto tax = &db.ranked_lca_of_targets(cand.target_id(0),
                                             cand.target_id(1));

        //classify if rank is below or at the highest rank of interest
        if(tax->rank <= highestRank) return tax;
    }
    else {
        if(lowestRank == taxon_rank::Sequence)
            lowestRank = taxonomy::next_main_rank(lowestRank);

//        std::cout << "vote (" << trustedMajority << "): ";

        std::unordered_map<taxon_id,int> scores;
        scores.rehash(2*cand.count());

        for(auto r = lowestRank; r <= highestRank; r = taxonomy::next_main_rank(r)) {
            //hash-count taxon id occurrences on rank 'r'
            int totalscore = 0;
            for(int i = 0, n = cand.count(); i < n; ++i) {
                //use target id instead of taxon if at sequence level
                auto taxid = db.ranks_of_target(cand.target_id(i))[int(r)];
                auto score = cand.hits(i);
                totalscore += score;
                auto it = scores.find(taxid);
                if(it != scores.end()) {
                    it->second += score;
                } else {
                    scores.insert(it, {taxid, score});
                }
            }

//            std::cout << "\n    " << taxonomy::rank_name(r) << " ";

            //determine taxon id with most votes
            taxon_id toptid = 0;
            int topscore = 0;
            for(const auto& x : scores) {

//                std::cout << x.first << ":" << x.second << ", ";

                if(x.second > topscore) {
                    toptid = x.first;
                    topscore = x.second;
                }
            }

            //if enough candidates (weighted by their hits)
            //agree on a taxon => classify as such
            if(toptid > 0 && topscore >= (totalscore * trustedMajority)) {

//                std::cout << "  => classified " << taxonomy::rank_name(r)
//                          << " " << toptid << '\n';

                return &db.taxon_with_id(toptid);
            }

            scores.clear();
        }

//        std::cout << '\n';
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
template<int n>
classification
sequence_classification(
    const database& db, const query_options& param,
    const matches_in_contiguous_window_range_top<n,target_id>& cand)
{
    auto wc = param.weightHitsWithWindows
                ? (cand.covered_windows() > 1) ? cand.covered_windows()-1 : 1
                : 1;

    //sum of top-2 hits < threshold => considered not classifiable
    if((cand.hits(0) + cand.hits(1)) < wc*param.hitsMin) {
        return classification{};
    }

    //either top 2 are the same sequences with at least 'hitsMin' many hits
    //(checked before) or hit difference between these top 2 is above threshhold
    if( (cand.target_id(0) == cand.target_id(1))
        || (cand.hits(0) - cand.hits(1)) >= wc*param.hitsMin)
    {
        //return top candidate
        auto tid = cand.target_id(0);
        return classification{tid, &db.taxon_of_target(tid)};
    }

    return classification{
        lowest_common_taxon(db, cand, param.hitsDiff,
                               param.lowestRank, param.highestRank) };
}



/*****************************************************************************
 *
 * @brief prints a classification to stream 'os'
 *
 *****************************************************************************/
void show_classification(std::ostream& os,
                         const database& db, const query_options& param,
                         const classification& cls)
{
    if(cls.sequence_level()) {
        auto rmax = param.showLineage ? param.highestRank : param.lowestRank;

        show_ranks_of_target(os, db, cls.target(),
                             param.showTaxaAs, param.lowestRank, rmax);
    }
    else if(cls.has_taxon()) {
        if(cls.rank() > param.highestRank) {
            os << "--";
        }
        else {
            auto rmin = param.lowestRank < cls.rank() ? cls.rank() : param.lowestRank;
            auto rmax = param.showLineage ? param.highestRank : rmin;

            show_ranks(os, db, db.ranks(cls.tax()),
                       param.showTaxaAs, rmin, rmax);
        }
    }
    else {
        os << "--";
    }
}



/*****************************************************************************
 *
 * @brief removes hits of a taxon on a certain rank; can be very slow!
 *
 *****************************************************************************/
void remove_hits_on_rank(const database& db,
                         taxon_rank rank, taxon_id taxid,
                         match_result& res)
{
    auto maskedRes = res;
    for(const auto& orig : res) {
        auto r = db.ranks_of_target(orig.first.tgt)[int(rank)];
        if(r != taxid) {
            maskedRes.insert(orig);
        }
    }
    using std::swap;
    swap(maskedRes,res);
}



/*****************************************************************************
 *
 * @brief process database hit results: classify & print mapping
 *
 *****************************************************************************/
void process_database_answer(
     const database& db, const query_options& param,
     const std::string& header,
     const sequence& query1, const sequence& query2,
     match_result&& hits, std::ostream& os, classification_statistics& stats)
{
    if(header.empty()) return;

    //preparation -------------------------------
    classification groundTruth;
    if(param.testPrecision ||
       (param.mapViewMode != map_view_mode::none && param.showGroundTruth) ||
       (param.excludedRank != taxon_rank::none) )
    {
        groundTruth = ground_truth(db, header);
    }
    //clade exclusion
    if(param.excludedRank != taxon_rank::none && groundTruth.has_taxon()) {
        auto exclTaxid = db.ranks(groundTruth.tax())[int(param.excludedRank)];
        remove_hits_on_rank(db, param.excludedRank, exclTaxid, hits);
    }

    //classify ----------------------------------
    auto numWindows = std::uint_least64_t( 2 + (
        std::max(query1.size() + query2.size(), param.insertSizeMax) /
        db.target_window_stride() ));

    top_matches_in_contiguous_window_range tophits {hits, numWindows};
    auto cls = sequence_classification(db, param, tophits);

    //evaluate classification -------------------
    if(param.testPrecision) {
        auto lowestCorrectRank = lowest_common_rank(db, cls, groundTruth);

        stats.assign_known_correct(cls.rank(), groundTruth.rank(),
                                   lowestCorrectRank);

        //check if taxa of assigned target are covered
        if(param.testCoverage && groundTruth.has_taxon()) {
            update_coverage_statistics(db, cls, groundTruth, stats);
        }

    } else {
        stats.assign(cls.rank());
    }

    //show mapping ------------------------------
    bool showMapping = param.mapViewMode == map_view_mode::all ||
        (param.mapViewMode == map_view_mode::mapped_only && !cls.none());

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
        os << param.outSeparator;

        if(param.showGroundTruth) {
            if(groundTruth.sequence_level()) {
                show_ranks_of_target(os, db, groundTruth.target(),
                    param.showTaxaAs, param.lowestRank,
                    param.showLineage ? param.highestRank : param.lowestRank);
            }
            else if(groundTruth.has_taxon()) {
                show_ranks(os, db, db.ranks(groundTruth.tax()),
                    param.showTaxaAs, groundTruth.tax().rank,
                    param.showLineage ? param.highestRank : groundTruth.tax().rank);
            }
            else {
                os << "n/a";
            }
            os << param.outSeparator;
        }

        //print results
        if(param.showAllHits) {
            show_matches(os, db, hits, param.lowestRank);
            os << param.outSeparator;
        }
        if(param.showTopHits) {
            show_matches(os, db, tophits, param.lowestRank);
            os << param.outSeparator;
        }
        if(param.showLocations) {
            show_candidate_ranges(os, db, tophits);
            os << param.outSeparator;
        }
        show_classification(os, db, param, cls);
    }

    //optional alignment ------------------------------
    if(param.testAlignment && !cls.none()) {
        //check which of the top sequences has a better alignment
        const auto& origin = db.origin_of_target(tophits.target_id(0));

        try {
            //load candidate file and forward to sequence
            auto reader = make_sequence_reader(origin.filename);
            for(std::size_t i = 0; i < origin.index; ++i) reader->next();
            if(reader->has_next()) {
                auto tgtSequ = reader->next().data;
                auto subject = make_view_from_window_range(tgtSequ,
                                    tophits.window(0),
                                    db.target_window_stride());

                auto align = make_alignment(query1, query2, subject,
                                            needleman_wunsch_scheme{});

                stats.register_alignment_score(align.score);

                //print alignment to top candidate
                if(align.score > 0 && param.showAlignment) {
                    const auto w = db.target_window_stride();

                    os << '\n'
                       << param.comment << "  score  " << align.score
                       << "  aligned to "
                       << origin.filename << " #" << origin.index
                       << " in range [" << (w * tophits.window(0).beg)
                       << ',' << (w * tophits.window(0).end) << "]\n"
                       << param.comment << "  query  " << align.query << '\n'
                       << param.comment << "  target " << align.subject;
                }
            }
        } catch(std::exception&) { }
    }

    if(showMapping) os << '\n';
}



/*****************************************************************************
 *
 * @brief classifies single sequences from a single source
 *
 *****************************************************************************/
void classify(parallel_queue& queue,
    const database& db, const query_options& param,
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

                        process_database_answer(db, param,
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
            classify_pairs(queue, db, param,
                           *make_sequence_reader(fname1),
                           *make_sequence_reader(fname2),
                           os, stats);
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

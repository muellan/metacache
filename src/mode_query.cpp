/*****************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * version 1.0
 *
 * Copyright (C) 2016 André Müller (muellan@uni-mainz.de)
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

#include "args_handling.h"
#include "timer.h"
#include "filesys_utility.h"
#include "cmdline_utility.h"
#include "sequence_io.h"
#include "sequence_view.h"
#include "alignment.h"

#include "modes_common.h"


namespace mc {



/*****************************************************************************
 * @brief query options
 *****************************************************************************/
enum class pairing_mode : unsigned char {
    none, files, sequences
};


/*****************************************************************************
 * @brief query options
 *****************************************************************************/
enum class align_mode : unsigned char {
    none, semi_global
};



/*****************************************************************************
 *
 * @brief query options
 *
 *****************************************************************************/
struct query_param
{
    //-----------------------------------------------------
    //output options & formatting
    //-----------------------------------------------------
    //make a separate output file for each input file
    bool splitOutput = false;
    //show top candidate sequences and their associated k-mer hash hit count
    bool showTopHits = false;
    //show all k-mer-hash hits in database for each given read
    bool showAllHits = false;
    //show classification (read mappings), if false, only summary will be shown
    bool showMappings = false;
    //show candidate position(s) in reference sequence(s)
    bool showPosition = false;
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

    //-------------------------------------------
    //classification options
    //-----------------------------------------------------
    pairing_mode pairing = pairing_mode::none;
    //ranks/taxa to classify on and to show
    taxon_rank lowestRank  = taxon_rank::Sequence;
    taxon_rank highestRank = taxon_rank::Domain;
    //ground truth rank to exclude (for clade exclusion test)
    taxon_rank excludedRank = taxon_rank::none;
    //
    int hitsMin  = -1;  //< 0 : deduced from database parameters
    float hitsDiff = 0.5;
    //use full kmer counting (if sequence files available) for ambiguous cases
    int useCommonKmerCount = 0;
    float kmerCountMin = 0.9f;
    float kmerCountDiff = 1.1f;
    //use alignment (if sequence files available) for ambiguous cases
    bool useAlignment = false;
    float alignmentMin  = 0.9f;
    float alignmentDiff = 1.1f;
    //maximum range in sequence that read (pair) is expected to be in
    std::size_t insertSizeMax = 0;

    //-----------------------------------------------------
    //analysis options
    //-----------------------------------------------------
    //test precision (ground truth must be available)
    bool testPrecision = false;
    bool testCoverage = false;
    //additional file with query -> ground truth mapping
    std::string sequ2taxonPreFile;

    //-------------------------------------------
    //query sampling scheme
    //-----------------------------------------------------
    int sketchlen = -1;  //< 0 : use value from database
    int winlen = -1;     //< 0 : use value from database
    int winstride = -1;  //< 0 : use value from database

    //-------------------------------------------
    //database tuning parameters
    //-----------------------------------------------------
    float maxLoadFactor = -1;        //< 0 : use value from database
    int maxGenomesPerSketchVal = -1; //< 0 : use value from database

    //-------------------------------------------
    //filenames
    std::string dbfile;
    std::vector<std::string> infiles;
    std::string outfile;
};




/*****************************************************************************
 *
 * @brief command line args -> query options
 *
 *****************************************************************************/
query_param
get_query_param(const args_parser& args)
{
    const query_param defaults;

    query_param param;

    //files
    param.dbfile = database_name(args);

    param.infiles = sequence_filenames(args);
    if(param.infiles.empty()) {
        throw std::invalid_argument{"No sequence filenames provided."};
    }

    //pairing
    if(args.contains("paired_files") || args.contains("pair_files")) {
        if(param.infiles.size() > 1) {
            param.pairing = pairing_mode::files;
            std::sort(param.infiles.begin(), param.infiles.end());
        }
    }
    else if(args.contains("paired_sequences") || args.contains("pair_sequences")) {
        param.pairing = pairing_mode::sequences;
    }

    //analysis options
    param.testCoverage = args.contains("coverage");
    param.testPrecision = param.testCoverage || args.contains("precision");

    //classification ranks
    auto lowestRank = args.get<std::string>("lowest", "");
    if(!lowestRank.empty()) {
        auto r = taxonomy::rank_from_name(lowestRank);
        if(r < taxonomy::taxon_rank::root) {
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
    param.hitsMin  = args.get<int>("hitmin", defaults.hitsMin);
    param.hitsDiff = args.get<float>("hitdiff", defaults.hitsDiff);

    //kmer counting
    param.useCommonKmerCount = args.get<int>("kmercount", defaults.useCommonKmerCount);
    if(param.useCommonKmerCount < 0)
        param.useCommonKmerCount = 0;
    else if(param.useCommonKmerCount > 32)
        param.useCommonKmerCount = 32;

    param.kmerCountMin  = args.get<float>("kcmin", defaults.kmerCountMin);
    param.kmerCountDiff = args.get<float>("kcdiff", defaults.kmerCountDiff);

    //alignment
    param.useAlignment  = args.contains("align");
    param.alignmentMin  = args.get<float>("alignmin", defaults.alignmentMin);
    param.alignmentDiff = args.get<float>("aligndiff", defaults.alignmentDiff);
    param.showAlignment = args.contains("showalign");

    //output formatting
    param.showLineage = args.contains("lineage");

    param.outfile = args.get<std::string>("out", "");
    param.splitOutput = args.contains("splitout");

    param.outSeparator = args.get<std::string>("separator", "\t|\t");

    param.showPosition = args.contains("positions");

    param.showTopHits = args.contains("tophits");
    param.showAllHits = args.contains("allhits");

    if(args.contains("taxidsonly") ||
       args.contains("taxids_only") ||
       args.contains("taxid_only") ||
       args.contains("taxidonly"))
    {
        param.showTaxaAs = taxon_print_mode::id_only;
    }
    else if(args.contains("taxids") || args.contains("taxid")) {
        param.showTaxaAs = taxon_print_mode::id_name;
    }
    else {
        param.showTaxaAs = taxon_print_mode::name_only;
    }

    param.showMappings =  param.showTopHits ||
                          param.showAllHits ||
                          (!args.contains("nomap") &&
                           !args.contains("noshowmap") &&
                           !args.contains("nomapping") &&
                           !args.contains("nomappings") );

    param.showGroundTruth = args.contains("ground_truth");

    param.insertSizeMax = args.get<std::size_t>("insertsize",
                                                defaults.insertSizeMax);

    //database tuning parameters
    param.maxLoadFactor = args.get<float>("max_load_fac", defaults.maxLoadFactor);

    param.maxGenomesPerSketchVal = args.get<int>("max_genomes_per_feature",
                                                 defaults.maxGenomesPerSketchVal);

    return param;
}



/*****************************************************************************
 *
 * @brief  print genome.window hit statistics from database
 *
 *****************************************************************************/
template<int n>
void show_candidate_ranges(std::ostream& os,
    const database& db,
    const matches_in_contiguous_window_range_top<n>& cand)
{
    const auto w = db.genome_window_stride();

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
 * @param  query1  read
 * @param  query2  paired read
 * @return index of best candidate or -1 if ambiguous or unsatisfactory
 *
 *****************************************************************************/
template<int n>
int
best_kmer_intersection_candidate(std::ostream&,
    const database& db, const query_param& param,
    const sequence& query1, const sequence& query2,
    const matches_in_contiguous_window_range_top<n>& cand)
{
    std::size_t fstVal = 0;
    std::size_t sndVal = 0;
    int fstIdx = -1;
    int sndIdx = -1;

//    std::cout << '\n' << param.comment
//              << "  " << param.fullKmerCount << "-mer counts: ";

    //largest kmer intersection
    for(int i = 0; i < n; ++i) {
        try {
            const auto& origin = db.origin_of_genome(cand.genome_id(i));
            if(!origin.filename.empty()) {
                //open sequence file, forward to index in file
                auto reader = make_sequence_reader(origin.filename);
                for(std::size_t j = 0; j < origin.index; ++j) reader->next();
                if(reader->has_next()) {

                    auto subject = reader->next().data;
                    //make views to candidate window ranges
                    auto w = db.genome_window_stride();
                    auto win = make_view_from_window_range(subject, cand.window(i), w);

                    auto s = kmer_intersection_size(param.useCommonKmerCount, query1, win);
                    if(!query2.empty()) {
                        s += kmer_intersection_size(param.useCommonKmerCount, query2, win);
                    }

//                    std::cout << s << ' ';

                    if(s > fstVal) {
                        sndVal = fstVal;
                        sndIdx = fstIdx;
                        fstVal = s;
                        fstIdx = i;
                    }
                    else if(s > sndVal) {
                        sndVal = s;
                        sndIdx = i;
                    }
                }
            }
        }
        catch(std::exception&) {}
    }
//    std::cout << '\n';

    if(fstIdx < 0 || sndIdx < 0) return -1;

    //isec size must exceed a minimum fraction of the largest possible size
    auto maxIsecSize = query1.size() + query2.size() - param.useCommonKmerCount + 1;
    if(fstVal < param.kmerCountMin * maxIsecSize) return -1;

    //largest isec size significantly greater thean second largest?
    return (fstVal >= param.kmerCountDiff * sndVal) ? fstIdx : sndIdx;
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
template<int n>
int
best_alignment_candidate(std::ostream& os,
    const database& db, const query_param& param,
    const sequence& query1, const sequence& query2,
    const matches_in_contiguous_window_range_top<n>& cand)
{
    //check which of the top sequences has a better alignment
    const auto& origin0 = db.origin_of_genome(cand.genome_id(0));
    if(origin0.filename.empty()) return -1;
    const auto& origin1 = db.origin_of_genome(cand.genome_id(1));
    if(origin1.filename.empty()) return -1;

    //open sequence files & forward to sequence index
    try {
        //load candidate sequences
        auto reader0 = make_sequence_reader(origin0.filename);
        for(std::size_t i = 0; i < origin0.index; ++i) reader0->next();
        if(!reader0->has_next()) return -1;

        auto reader1 = make_sequence_reader(origin1.filename);
        for(std::size_t i = 0; i < origin1.index; ++i) reader1->next();
        if(!reader1->has_next()) return -1;

        auto subject0 = reader0->next().data;
        auto subject1 = reader1->next().data;

        //make views to candidate window ranges
        auto w = db.genome_window_stride();

        auto s0 = make_view_from_window_range(subject0, cand.window(0), w);
        auto s1 = make_view_from_window_range(subject1, cand.window(1), w);

        using scheme_t = needleman_wunsch_scheme;

        std::size_t als0  = 0; std::size_t als0r = 0;
        std::size_t als1  = 0; std::size_t als1r = 0;

        if(!param.showAlignment) {
            //compute scores only
            als0 = semi_global_alignment_score(query1, s0, scheme_t{});
            als1 = semi_global_alignment_score(query1, s1, scheme_t{});
            auto query1r = make_reverse_complement(query1);
            als0r = semi_global_alignment_score(query1r, s0, scheme_t{});
            als1r = semi_global_alignment_score(query1r, s1, scheme_t{});

            //align paired read as well
            if(!query2.empty()) {
                als0 += semi_global_alignment_score(query2, s0, scheme_t{});
                als1 += semi_global_alignment_score(query2, s1, scheme_t{});
                auto query2r = make_reverse_complement(query2);
                als0r += semi_global_alignment_score(query2r, s0, scheme_t{});
                als1r += semi_global_alignment_score(query2r, s1, scheme_t{});
            }
            if(als0r > als0) als0 = als0r;
            if(als1r > als1) als1 = als1r;
        }
        else {
            //compute alignment
            auto align0 = semi_global_alignment(query1, s0, scheme_t{}, true);
            auto align1 = semi_global_alignment(query1, s1, scheme_t{}, true);
            als0 = align0.score;
            als1 = align1.score;
            //reverse complement
            auto query1r = make_reverse_complement(query1);
            auto align0r = semi_global_alignment(query1r, s0, scheme_t{}, true);
            auto align1r = semi_global_alignment(query1r, s1, scheme_t{}, true);
            als0r = align0r.score;
            als1r = align1r.score;

            //align paired read as well
            if(!query2.empty()) {
                als0 += semi_global_alignment(query2, s0, scheme_t{}).score;
                als1 += semi_global_alignment(query2, s1, scheme_t{}).score;
                auto query2r = make_reverse_complement(query2);
                als0r += semi_global_alignment(query2r, s0, scheme_t{}).score;
                als1r += semi_global_alignment(query2r, s1, scheme_t{}).score;
            }

//            os     << '\n'
//                   << param.comment << "Alignement candidates:\n"
//                   << param.comment << "  Query:  " << query1 << '\n';
//            if(!query2.empty()) {
//                os << param.comment << "     +    " << query2 << '\n';
//            }
//            os     << param.comment << "  S1:     " << s0 << '\n'
//                   << param.comment << "  S2:     " << s1 << '\n'
            os     << '\n'
                   << param.comment << "  Alignments:\n";

            if(als0 > als0r) {
                os << param.comment << "    score " << align0.score << '\n'
                   << param.comment << "       Q: " << align0.query << '\n'
                   << param.comment << "      S1: " << align0.subject << '\n';
            } else {
                os << param.comment << "    score " << align0r.score << '\n'
                   << param.comment << "    (r)Q: " << align0r.query << '\n'
                   << param.comment << "      S1: " << align0r.subject << '\n';
            }
            if(als1 > als1r) {
                os << param.comment << "    score " << align1.score << '\n'
                   << param.comment << "       Q: " << align1.query << '\n'
                   << param.comment << "      S2: " << align1.subject << '\n';
            } else {
                os << param.comment << "    score " << align1r.score << '\n'
                   << param.comment << "    (r)Q: " << align1r.query << '\n'
                   << param.comment << "      S2: " << align1r.subject << '\n';
            }

            if(als0r > als0) als0 = als0r;
            if(als1r > als1) als1 = als1r;
        }

        auto minAl = param.alignmentMin * scheme_t::max_score() *
                     (query1.size() + query2.size());

        if(als0 >= minAl && als0 >= param.alignmentDiff * als1) return 0;
        if(als1 >= minAl && als1 >= param.alignmentDiff * als0) return 1;
    }
    catch (std::exception&) {
        //file I/O error
    }
    return -1;
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
    const matches_in_contiguous_window_range_top<maxn>& cand,
    float trustedMajority,
    taxon_rank lowestRank = taxon_rank::subSpecies,
    taxon_rank highestRank = taxon_rank::Domain)
{
    if(maxn < 3 || cand.count() < 3) {
        auto tax = &db.ranked_lca_of_genomes(cand.genome_id(0),
                                             cand.genome_id(1));

        //classify if rank is below or at the highest rank of interest
        if(tax->rank <= highestRank) return tax;
    }
    else {
        if(lowestRank == taxon_rank::Sequence) ++lowestRank;

//        std::cout << "vote (" << trustedMajority << "): ";

        std::unordered_map<taxon_id,int> scores;
        scores.rehash(2*cand.count());

        for(auto r = lowestRank; r <= highestRank; r = taxonomy::next_main_rank(r)) {
            //hash-count taxon id occurrences on rank 'r'
            int totalscore = 0;
            for(int i = 0, n = cand.count(); i < n; ++i) {
                //use genome id instead of taxon if at sequence level
                auto tid = db.ranks_of_genome(cand.genome_id(i))[int(r)];
                if(tid > 0) {
                    auto score = cand.hits(i);
                    totalscore += score;
                    auto it = scores.find(tid);
                    if(it != scores.end()) {
                        it->second += score;
                    } else {
                        scores.insert(it, {tid, score});
                    }
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
            if(topscore >= (totalscore * trustedMajority)) {

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
sequence_classification(std::ostream& os,
    const database& db, const query_param& param,
    const sequence& query1, const sequence& query2,
    const matches_in_contiguous_window_range_top<n>& cand)
{
    //sum of top-2 hits < threshold => considered not classifiable
    if((cand.hits(0) + cand.hits(1)) < param.hitsMin) {
        return classification{};
    }

    //either top 2 are the same sequences with at least 'hitsMin' many hits
    //(checked before) or hit difference between these top 2 is above threshhold
    if( (cand.genome_id(0) == cand.genome_id(1))
        || (cand.hits(0) - cand.hits(1)) >= param.hitsMin)
    {
        //return top candidate
        auto gid = genome_id(cand.genome_id(0));
        return classification{gid, &db.taxon_of_genome(gid)};
    }

    //use kmer counting to disambiguate?
    if(param.useCommonKmerCount > 0) {
        int bac = best_kmer_intersection_candidate(os, db, param,
                                                   query1, query2, cand);
        if(bac >= 0) {
            auto gid = genome_id(cand.genome_id(bac));
            return classification{gid, &db.taxon_of_genome(gid)};
        }
    }

    //do alignment of query to top candidate sequences?
    if(param.useAlignment) {
        int bac = best_alignment_candidate(os, db, param, query1, query2, cand);
        if(bac >= 0) {
            auto gid = genome_id(cand.genome_id(bac));
            return classification{gid, &db.taxon_of_genome(gid)};
        }
    }

    return classification{
        lowest_common_taxon(db, cand, param.hitsDiff,
                               param.lowestRank, param.highestRank) };
}



/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
void show_classification(std::ostream& os,
                         const database& db, const query_param& param,
                         const classification& cls)
{
    if(cls.sequence_level()) {
        auto rmax = param.showLineage ? param.highestRank : param.lowestRank;

        show_ranks_of_genome(os, db, cls.gid(),
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
 * @brief  clade exclusion
           masks out hits of a taxon; can be very slow!
 *
 *****************************************************************************/
void remove_hits_on_rank(const database& db,
                         taxon_rank rank, taxon_id taxid,
                         match_result& res)
{
    auto maskedRes = res;
    for(const auto& orig : res) {
        auto r = db.ranks_of_genome(orig.first.gid)[int(rank)];
        if(r != taxid) {
            maskedRes.insert(orig);
        }
    }
    using std::swap;
    swap(maskedRes,res);
}



/*****************************************************************************
 *
 * @brief process database hit results
 *
 *****************************************************************************/
void process_database_answer(
     const database& db, const query_param& param,
     const std::string& header,
     const sequence& query1, const sequence& query2,
     match_result&& hits, std::ostream& os, rank_statistics& stats)
{
    classification groundTruth;
    if(param.testPrecision ||
       (param.showMappings && param.showGroundTruth) ||
       (param.excludedRank != taxon_rank::none) )
    {
        groundTruth = ground_truth(db, header);
    }

    //preprocess
    if(param.excludedRank != taxon_rank::none && groundTruth.has_taxon()) {
        auto exclTaxid = db.ranks(groundTruth.tax())[int(param.excludedRank)];
        remove_hits_on_rank(db, param.excludedRank, exclTaxid, hits);
    }

    //print query header and ground truth
    if(param.showMappings || param.showTopHits || param.showAllHits) {
        //show first contiguous string only
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
                show_ranks_of_genome(os, db, groundTruth.gid(),
                    param.showTaxaAs, param.lowestRank,
                    param.showLineage ? param.highestRank : param.lowestRank);
            }
            else if(groundTruth.has_taxon()) {
                show_ranks(os, db, db.ranks(groundTruth.tax()),
                    param.showTaxaAs, param.lowestRank,
                    param.showLineage ? param.highestRank : param.lowestRank);
            }
            else {
                os << "n/a";
            }
            os << param.outSeparator;
        }
    }

    //classify
    int numWindows = 2 + (
        std::max(query1.size() + query2.size(), param.insertSizeMax) /
        db.genome_window_stride() );

    top_matches_in_contiguous_window_range tophits {hits, numWindows};
    auto cls = sequence_classification(os, db, param, query1, query2, tophits);

    //print results
    if(param.showMappings) {
        if(param.showAllHits) {
            show_matches(os, db, hits, param.lowestRank);
            os << param.outSeparator;
        }
        if(param.showTopHits) {
            show_matches(os, db, tophits, param.lowestRank);
            os << param.outSeparator;
        }
        if(param.showPosition) {
            show_candidate_ranges(os, db, tophits);
            os << param.outSeparator;
        }
        show_classification(os, db, param, cls);
    }

    if(param.testPrecision) {
        auto lowestCorrectRank = lowest_common_rank(db, cls, groundTruth);

        stats.assign_known_correct(cls.rank(), groundTruth.rank(),
                                   lowestCorrectRank);

        //check if taxa of assigned genome are covered
        if(param.testCoverage && groundTruth.has_taxon()) {
            update_coverage_statistics(db, cls, groundTruth, stats);
        }
    } else {
        stats.assign(cls.rank());
    }

    if(param.showMappings) {
        os << '\n';
    }
}



/*****************************************************************************
 *
 * @brief classifies paired-end sequences in file pairs
 *
 *****************************************************************************/
void classify_on_file_pairs(
    const database& db, const query_param& param,
    const std::vector<std::string>& infilenames,
    std::ostream& os, rank_statistics& stats)
{
    //pair up reads from two consecutive files in the list
    for(size_t i = 0; i < infilenames.size(); i += 2) {
        const auto& fname1 = infilenames[i];
        const auto& fname2 = infilenames[i + 1];
        if(param.showMappings) {
            os << param.comment << fname1 << " + " << fname2 << '\n';
        }
        else if(!param.outfile.empty() && infilenames.size() > 1) {
            show_progress_indicator(i / float(infilenames.size()));
        }

        try {
            auto reader1 = make_sequence_reader(fname1);
            auto reader2 = make_sequence_reader(fname2);
            while(reader1->has_next() && reader2->has_next()) {
                auto query1 = reader1->next();
                auto query2 = reader2->next();
                if(!query1.data.empty() && !query2.data.empty() ) {
                    //merge sketches of both reads
                    auto res = db.matches(query1.data);
                    db.accumulate_matches(query2.data, res);

                    process_database_answer(db, param,
                        query1.header, query1.data, query2.data,
                        std::move(res), os, stats);
                }
            }
        }
        catch(std::exception& e) {
            if(param.showMappings) {
                os << "FAIL: " << e.what() << std::endl;
            }
        }
    }
}



/*****************************************************************************
 *
 * @brief classify sequences per input file
 *
 *****************************************************************************/
void classify_per_file(
    const database& db, const query_param& param,
    const std::vector<std::string>& infilenames,
    std::ostream& os, rank_statistics& stats)
{

    for(size_t i = 0; i < infilenames.size(); ++i) {
        const auto& fname = infilenames[i];
        if(param.showMappings) {
            os << param.comment << fname << '\n';
        }
        else if(!param.outfile.empty() && infilenames.size() > 1) {
            show_progress_indicator(i / float(infilenames.size()));
        }

        try {
            auto reader = make_sequence_reader(fname);
            while(reader->has_next()) {
                auto query = reader->next();
                if(!query.data.empty()) {
                    //make sure we use the first read's header in case of pairing
                    auto res = db.matches(query.data);

                    //pair up with next read?
                    if(param.pairing == pairing_mode::sequences &&
                        reader->has_next())
                    {
                        auto query1 = std::move(query.data);
                        auto header = std::move(query.header);

                        query = reader->next();
                        db.accumulate_matches(query.data, res);

                        process_database_answer(db, param,
                            header, query1, query.data,
                            std::move(res), os, stats);
                    }
                    else {
                        process_database_answer(db, param,
                            query.header, query.data, sequence{},
                            std::move(res), os, stats);
                    }
                }
            }
        }
        catch(std::exception& e) {
            if(param.showMappings) {
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
void classify_sequences(const database& db,
                        const query_param& param,
                        const std::vector<std::string>& infilenames,
                        std::ostream& os)
{
    const auto& comment = param.comment;

    if(param.showMappings) {
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
    if(param.useAlignment > 0) {
        os << comment << "Classification considers top "
           << top_matches_in_contiguous_window_range::max_count() << " matches\n";
    }
    if(param.useCommonKmerCount > 0) {
        os << comment << "Classification uses all "
           << param.useCommonKmerCount << "-mer counting." "\n";
    }
    if(param.useAlignment > 0) {
        os << comment << "Classification uses semi-global alignment.\n";
    }
    if(param.outfile.empty()) os.flush();

    timer time;
    time.start();

    rank_statistics stats;
    if(param.pairing == pairing_mode::files) {
        classify_on_file_pairs(db, param, infilenames, os, stats);
    }
    else {
        classify_per_file(db, param, infilenames, os, stats);
    }
    if(!param.outfile.empty() || !param.showMappings) {
        clear_current_line();
    }

    time.stop();

    //show results
    const auto speed = stats.total() / time.minutes();
    os << comment << "queries: " << stats.total() << '\n'
       << comment << "time:    " << time.milliseconds() << " ms\n"
       << comment << "speed:   " << speed << " queries/min\n";

    if(stats.total() > 0) {
        show_per_rank_statistics(os, stats, comment);
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
void process_input_files(const database& db,
                         const query_param& param,
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
                      << param.outfile << std::endl;

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

    auto param = get_query_param(args);

    auto db = make_database<database>(param.dbfile);


    //configure database
    if(param.maxLoadFactor > 0) {
        db.max_load_factor(param.maxLoadFactor);
    }
    if(param.maxGenomesPerSketchVal > 1) {
        db.erase_features_with_more_genomes_than(param.maxGenomesPerSketchVal);
    }

    //deduced query parameters
    if(param.hitsMin < 1) {
        param.hitsMin = db.genome_sketcher().sketch_size() / 2 - 1;
        if(param.hitsMin < 1) param.hitsMin = 1;
    }

    //use a different sketching scheme for querying?
    if(param.winlen > 0)    db.query_window_size(param.winlen);
    if(param.winstride > 0) db.query_window_stride(param.winstride);

    if(param.sketchlen > 0) {
        auto s = db.query_sketcher();
        s.sketch_size(param.sketchlen);
        db.query_sketcher(std::move(s));
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

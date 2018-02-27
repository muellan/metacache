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
#include "query_options.h"
#include "querying.h"
#include "classification.h"
#include "print_results.h"


namespace mc {

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::flush;


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

        query_with_file_pairs(queue, db, opt, infilenames, os,
            [&](const string& header, const sequence& seq1, const sequence& seq2,
                matches_per_location&& matches, std::ostringstream& obuf)
            {
                process_database_answer(db, opt,
                                        header, seq1, seq2,
                                        std::move(matches), stats, obuf);
            });
    }
    else {
        query_per_file(queue, db, opt, infilenames, os,
            [&](const string& header, const sequence& seq1, const sequence& seq2,
                matches_per_location&& matches, std::ostringstream& obuf)
            {
                process_database_answer(db, opt,
                                        header, seq1, seq2,
                                        std::move(matches), stats, obuf);
            });
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

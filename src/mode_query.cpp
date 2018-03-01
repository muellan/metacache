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

using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::flush;


/*****************************************************************************
 *
 * @brief prints classification parameters
 *
 *****************************************************************************/
void show_classification_parameters(const query_options& opt, std::ostream& os)
{
    const auto& comment = opt.output.comment;
    if(opt.output.mapViewMode != map_view_mode::none) {
        os << comment << "Reporting per-read mappings (non-mapping lines start with '"
           << comment << "').\n";

        os << comment
           << "Output will be constrained to ranks from '"
           << taxonomy::rank_name(opt.classify.lowestRank) << "' to '"
           << taxonomy::rank_name(opt.classify.highestRank) << "'.\n";

        if(opt.output.showLineage) {
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

    if(opt.test.excludedRank != taxon_rank::none) {
        os << comment << "Clade Exclusion on Rank: "
           << taxonomy::rank_name(opt.test.excludedRank);
    }

    if(opt.process.pairing == pairing_mode::files) {
        os << comment << "File based paired-end mode:\n"
           << comment << "  Reads from two consecutive files will be interleaved.\n";
    }
    else if(opt.process.pairing == pairing_mode::sequences) {
        os << comment << "Per file paired-end mode:\n"
           << comment << "  Reads from two consecutive sequences in each file will be paired up.\n";
    }

    os << comment << "  Max insert size considered " << opt.classify.insertSizeMax << ".\n";

    if(opt.test.alignmentScores) {
        os << comment << "Query sequences will be aligned to best candidate target => SLOW!\n";
    }

    os << comment << "Using " << opt.process.numThreads << " threads\n";
}



/*************************************************************************//**
 *
 * @brief
 *
 *****************************************************************************/
void run_classification(const vector<string>& infiles,
                        const database& db, const query_options& opt,
                        std::ostream& output, std::ostream& status)
{
    if(opt.showClassificationParams) {
        show_classification_parameters(opt, output);
    }
    output.flush();
    status.flush();

    timer time;
    time.start();


    classification_statistics statistics;
    //per thread batch object storing / processing query results
    using buffer_t = std::ostringstream;
    const auto bufferSource = [] { return buffer_t{}; };

    //classify and write classification result to buffer
    const auto bufferUpdate = [&] (
            buffer_t& buf, matches_per_location&& matches,
            const string& header,
            const sequence& seq1, const sequence& seq2)
        {
            process_database_answer(db, opt, header, seq1, seq2,
                                    std::move(matches), statistics, buf);
        };

    //write results to output stream when batch is finished
    const auto bufferSink = [&] (const buffer_t& buf) {
        output << buf.str();
    };

    const auto displayStatus = [&] (const std::string& msg, float progress) {
        if(opt.output.mapViewMode != map_view_mode::none) {
            output << opt.output.comment << msg << '\n';
        }
        if(progress >= 0.0f) {
            show_progress_indicator(status, progress);
        }
    };

    query_database(infiles, db, opt.process,
                   bufferSource, bufferUpdate, bufferSink,
                   displayStatus);

    time.stop();

    clear_current_line(status);
    status.flush();

    //show results
    if(opt.showResults) {
        const auto numQueries = (opt.process.pairing == pairing_mode::none)
                                ? statistics.total() : 2 * statistics.total();

        const auto speed = numQueries / time.minutes();
        const auto& comment = opt.output.comment;
        output << comment << "queries: " << numQueries << '\n'
               << comment << "time:    " << time.milliseconds() << " ms\n"
               << comment << "speed:   " << speed << " queries/min\n";

        if(statistics.total() > 0) {
            show_classification_statistics(output, statistics, comment);
        } else {
            output << comment << "No valid query sequences found." << endl;
        }
    }

    output.flush();
}



/*****************************************************************************
 *
 * @brief runs classification on input files and writes the results
 *        to an output file or stdout
 *
 *****************************************************************************/
void process_input_files(const vector<string>& infiles,
                         const database& db, const query_options& opt,
                         const string& outfilename)
{
    if(outfilename.empty()) {
        run_classification(infiles, db, opt, cout, cerr);
    }
    else {
        std::ofstream os {outfilename, std::ios::out};

        if(os.good()) {
            cout << "Output will be redirected to file: " << outfilename << endl;

            run_classification(infiles, db, opt, os, cout);
        }
        else {
            throw file_write_error{
                "Could not write to output file " + outfilename};
        }
    }
}



/*****************************************************************************
 *
 * @brief runs classification on input files;
 *        handle output split
 *
 *****************************************************************************/
void process_input_files(const vector<string>& infiles,
                         const database& db,
                         const query_options& opt)
{
    if(infiles.empty()) {
        cerr << "No input filenames provided!" << endl;
        return;
    }
    else {
        bool anyreadable = false;
        for(const auto& f : infiles) {
            if(file_readable(f)) { anyreadable = true; break; }
        }
        if(!anyreadable) {
            throw std::runtime_error{
                        "None of the query sequence files could be opened"};
        }
    }

    //process files / file pairs separately
    if(opt.splitFiles) {
        string outfile;
        //process each input file pair separately
        if(opt.process.pairing == pairing_mode::files && infiles.size() > 1) {
            for(std::size_t i = 0; i < infiles.size(); i += 2) {
                const auto& f1 = infiles[i];
                const auto& f2 = infiles[i+1];
                if(!opt.outfile.empty()) {
                    outfile = opt.outfile
                            + "_" + extract_filename(f1)
                            + "_" + extract_filename(f2)
                            + ".txt";
                }
                process_input_files(vector<string>{f1,f2}, db, opt, outfile);
            }
        }
        //process each input file separately
        else {
            for(const auto& f : infiles) {
                if(!opt.outfile.empty()) {
                    outfile = opt.outfile + "_"
                            + extract_filename(f) + ".txt";
                }
                process_input_files(vector<string>{f}, db, opt, outfile);
            }
        }
    }
    //process all input files at once
    else {
        process_input_files(infiles, db, opt, opt.outfile);
    }
}



/*****************************************************************************
 *
 * @brief modify database content and sketching scheme according to
 *        command line options
 *
 *****************************************************************************/
void configure_database(database& db, const database_query_options& opt)
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
void adapt_options_to_database(classification_options& opt, const database& db)
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
void run_interactive_query_mode(const vector<string>& initInfiles,
                                const string& dbfile,
                                const database& db,
                                const query_options& initOpt)
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
            vector<string> arglist {"query", dbfile};
            std::istringstream iss(input);
            while(iss >> input) { arglist.push_back(input); }
            args_parser args{std::move(arglist)};

            //read command line options (use initial ones as defaults)
            auto infiles = sequence_filenames(args);
            if(!initInfiles.empty()) {
                infiles.insert(infiles.begin(), initInfiles.begin(),
                                                initInfiles.end());
            }
            auto opt = get_query_options(args, infiles, initOpt);
            adapt_options_to_database(opt.classify, db);

            if(opt.process.pairing == pairing_mode::files) {
                std::sort(infiles.begin(), infiles.end());
            }

            try {
                process_input_files(infiles, db, opt);
            }
            catch(std::exception& e) {
                if(initOpt.output.showErrors) cerr << e.what() << endl;
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
    auto infiles = sequence_filenames(args);
    auto opt = get_query_options(args,infiles);

    if(opt.process.pairing == pairing_mode::files) {
        std::sort(infiles.begin(), infiles.end());
    }

    auto dbfile = database_name(args);
    auto db = make_database<database>(dbfile);
    configure_database(db, get_database_query_options(args));

    if(opt.showDBproperties) {
        print_static_properties(db);
        print_content_properties(db);
    }

    if(!infiles.empty()) {
        cerr << "Classifying query sequences." << endl;

        adapt_options_to_database(opt.classify, db);
        process_input_files(infiles, db, opt);
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

        run_interactive_query_mode(infiles, dbfile, db, opt);
    }
}


} // namespace mc

/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2018 André Müller (muellan@uni-mainz.de)
 *                       & Robin Kobus  (rkobus@uni-mainz.de)
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

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "args_handling.h"
#include "cmdline_utility.h"
#include "filesys_utility.h"
#include "query_options.h"
#include "classification.h"
#include "classification_statistics.h"
#include "printing.h"
#include "config.h"


namespace mc {

using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::flush;



/*************************************************************************//**
 *
 * @brief runs classification on input files
 *
 *****************************************************************************/
void show_summary(const query_options& opt,
                  const classification_results& results)
{
    const auto& statistics = results.statistics;
    const auto numQueries = (opt.process.pairing == pairing_mode::none)
                            ? statistics.total() : 2 * statistics.total();

    const auto speed = numQueries / results.time.minutes();
    const auto& comment = opt.output.format.comment;
    results.mapout
        << comment << "queries: " << numQueries << '\n'
        << comment << "time:    " << results.time.milliseconds() << " ms\n"
        << comment << "speed:   " << speed << " queries/min\n";

    if(statistics.total() > 0) {
        show_taxon_statistics(results.mapout, statistics, comment);
    } else {
        results.status << comment << "No valid query sequences found." << endl;
    }
}



/*************************************************************************//**
 *
 * @brief runs classification on input files; sets output target streams
 *
 *****************************************************************************/
void process_input_files(const vector<string>& infiles,
                         const database& db, const query_options& opt,
                         const string& mapFilename,
                         const string& auxFilename)
{
    std::ostream* mapout = &cout;
    std::ostream* auxout = &cout;
    std::ostream* status = &cerr;

    std::ofstream mapFile;
    if(!mapFilename.empty()) {
        mapFile.open(mapFilename, std::ios::out);

        if(mapFile.good()) {
            cout << "Mappings will be written to file: " << mapFilename << endl;
            mapout = &mapFile;
            //default: auxiliary output same as mappings output
            auxout = mapout;
        }
        else {
            throw file_write_error{"Could not write to file " + mapFilename};
        }
    }

    std::ofstream auxFile;
    if(!auxFilename.empty()) {
        auxFile.open(auxFilename, std::ios::out);

        if(auxFile.good()) {
            cout << "Other output will be written to file: " << auxFilename << endl;
            auxout = &auxFile;
        }
        else {
            throw file_write_error{"Could not write to file " + auxFilename};
        }
    }

    classification_results results {*mapout,*auxout,*status};

    if(opt.showQueryParams) {
        show_query_parameters(results.mapout, opt);
    }

    results.flush_all_streams();

    results.time.start();
    map_queries_to_targets(infiles, db, opt, results);
    results.time.stop();

    clear_current_line(results.status);
    results.status.flush();

    if(opt.showSummary) show_summary(opt, results);

    results.flush_all_streams();
}



/*************************************************************************//**
 *
 * @brief runs classification on input files;
 *        handles output file split
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
        bool noneReadable = std::none_of(infiles.begin(), infiles.end(),
                           [](const auto& f) { return file_readable(f); });
        if(noneReadable) {
            throw std::runtime_error{
                        "None of the query sequence files could be opened"};
        }
    }

    //process files / file pairs separately
    if(opt.splitFiles) {
        string mapfile;
        string auxfile;
        //process each input file pair separately
        if(opt.process.pairing == pairing_mode::files && infiles.size() > 1) {
            for(std::size_t i = 0; i < infiles.size(); i += 2) {
                const auto& f1 = infiles[i];
                const auto& f2 = infiles[i+1];
                if(!opt.outfile.empty()) {
                    mapfile = opt.outfile
                            + "_" + extract_filename(f1)
                            + "_" + extract_filename(f2)
                            + ".txt";
                }
                if(!opt.auxfile.empty() && opt.auxfile != opt.outfile) {
                    auxfile = opt.auxfile
                            + "_" + extract_filename(f1)
                            + "_" + extract_filename(f2)
                            + ".txt";
                }
                process_input_files(vector<string>{f1,f2}, db, opt, mapfile, auxfile);
            }
        }
        //process each input file separately
        else {
            for(const auto& f : infiles) {
                if(!opt.outfile.empty()) {
                    mapfile = opt.outfile + "_"
                            + extract_filename(f) + ".txt";
                }
                if(!opt.auxfile.empty() && opt.auxfile != opt.outfile) {
                    auxfile = opt.auxfile + "_"
                            + extract_filename(f) + ".txt";
                }
                process_input_files(vector<string>{f}, db, opt, mapfile, auxfile);
            }
        }
    }
    //process all input files at once
    else {
        process_input_files(infiles, db, opt, opt.outfile, opt.auxfile);
    }
}



/*************************************************************************//**
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



/*************************************************************************//**
 *
 * @brief primitive REPL mode for repeated querying using the same database
 *
 *****************************************************************************/
void run_interactive_query_mode(const vector<string>& initInfiles,
                                const string& dbname,
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
            vector<string> arglist {"query", dbname};
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



/*************************************************************************//**
 *
 * @brief read database and modify db content and sketching scheme according to
 *        command line options
 *
 *****************************************************************************/
database
read_database(const string& filename, const database_query_options& opt)
{
    database db;

    if(opt.maxLoadFactor > 0.4 && opt.maxLoadFactor < 0.99) {
        db.max_load_factor(opt.maxLoadFactor);
        cerr << "Using custom hash table load factor of "
             << opt.maxLoadFactor << endl;
    }

    cerr << "Reading database from file '" << filename << "' ... " << flush;

    try {
        db.read(filename);
        cerr << "done." << endl;
    }
    catch(const file_access_error& e) {
        cerr << "FAIL" << endl;
        throw file_access_error{"Could not read database file '" + filename + "'"};
    }

    if(opt.removeOverpopulatedFeatures) {
        auto old = db.feature_count();

        auto maxlpf = opt.maxLocationsPerFeature - 1;
        if(maxlpf < 0 || maxlpf >= database::max_supported_locations_per_feature())
            maxlpf = database::max_supported_locations_per_feature() - 1;

        maxlpf = std::min(maxlpf, db.max_locations_per_feature() - 1);
        if(maxlpf > 0) { //always keep buckets with size 1
            cerr << "\nRemoving features with more than "
                 << maxlpf << " locations... " << std::flush;

            auto rem = db.remove_features_with_more_locations_than(maxlpf);

            cerr << rem << " of " << old << " removed." << endl;
        }
        //in case new max is less than the database setting
        db.max_locations_per_feature(opt.maxLocationsPerFeature);
    }
    else if(opt.maxLocationsPerFeature > 1) {
        db.max_locations_per_feature(opt.maxLocationsPerFeature);
        cerr << "max locations per feature set to "
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

    return db;
}



/*************************************************************************//**
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

    auto dbname = database_name(args);
    if(dbname.empty()) throw file_access_error{"No database name given"};

    auto db = read_database(dbname, get_database_query_options(args));

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

        run_interactive_query_mode(infiles, dbname, db, opt);
    }
}


} // namespace mc

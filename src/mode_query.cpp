/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2019 André Müller (muellan@uni-mainz.de)
 *                       & Robin Kobus  (kobus@uni-mainz.de)
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
 * @brief runs classification on input files; sets output target streams
 *
 *****************************************************************************/
void process_input_files(const vector<string>& infiles,
                         const database& db, const query_options& opt,
                         const string& queryMappingsFilename,
                         const string& targetsFilename,
                         const string& abundanceFilename)
{
    std::ostream* perReadOut   = &cout;
    std::ostream* perTargetOut = &cout;
    std::ostream* perTaxonOut  = &cout;
    std::ostream* status       = &cerr;

    std::ofstream mapFile;
    if(!queryMappingsFilename.empty()) {
        mapFile.open(queryMappingsFilename, std::ios::out);

        if(mapFile.good()) {
            cout << "Per-Read mappings will be written to file: " << queryMappingsFilename << endl;
            perReadOut = &mapFile;
            //default: auxiliary output same as mappings output
            perTargetOut = perReadOut;
            perTaxonOut  = perReadOut;
        }
        else {
            throw file_write_error{"Could not write to file " + queryMappingsFilename};
        }
    }

    std::ofstream targetsFile;
    if(!targetsFilename.empty()) {
        targetsFile.open(targetsFilename, std::ios::out);

        if(targetsFile.good()) {
            cout << "Per-Target mappings will be written to file: " << targetsFilename << endl;
            perTargetOut = &targetsFile;
        }
        else {
            throw file_write_error{"Could not write to file " + targetsFilename};
        }
    }

    std::ofstream abundanceFile;
    if(!abundanceFilename.empty()) {
        abundanceFile.open(abundanceFilename, std::ios::out);

        if(abundanceFile.good()) {
            cout << "Per-Taxon mappings will be written to file: " << abundanceFilename << endl;
            perTaxonOut = &abundanceFile;
        }
        else {
            throw file_write_error{"Could not write to file " + abundanceFilename};
        }
    }

    classification_results results {*perReadOut,*perTargetOut,*perTaxonOut,*status};

    if(opt.output.showQueryParams) {
        show_query_parameters(results.perReadOut, opt);
    }

    results.flush_all_streams();

    results.time.start();
    map_queries_to_targets(infiles, db, opt, results);
    results.time.stop();

    clear_current_line(results.status);
    results.status.flush();

    if(opt.output.showSummary) show_summary(opt, results);

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
    if(opt.output.splitFiles) {
        string queryMappingsFile;
        string targetsFile;
        string abundanceFile;
        //process each input file pair separately
        if(opt.process.pairing == pairing_mode::files && infiles.size() > 1) {
            for(std::size_t i = 0; i < infiles.size(); i += 2) {
                const auto& f1 = infiles[i];
                const auto& f2 = infiles[i+1];
                if(!opt.output.queryMappingsFile.empty()) {
                    queryMappingsFile = opt.output.queryMappingsFile
                            + "_" + extract_filename(f1)
                            + "_" + extract_filename(f2)
                            + ".txt";
                }
                if(!opt.output.targetsFile.empty() &&
                    opt.output.targetsFile != opt.output.queryMappingsFile)
                {
                    targetsFile = opt.output.targetsFile
                            + "_" + extract_filename(f1)
                            + "_" + extract_filename(f2)
                            + ".txt";
                }
                if(!opt.output.abundanceFile.empty() &&
                    opt.output.abundanceFile != opt.output.queryMappingsFile)
                {
                    abundanceFile = opt.output.abundanceFile
                            + "_" + extract_filename(f1)
                            + "_" + extract_filename(f2)
                            + ".txt";
                }
                process_input_files(vector<string>{f1,f2}, db, opt,
                                    queryMappingsFile, targetsFile, abundanceFile);
            }
        }
        //process each input file separately
        else {
            for(const auto& f : infiles) {
                if(!opt.output.queryMappingsFile.empty()) {
                    queryMappingsFile = opt.output.queryMappingsFile + "_"
                            + extract_filename(f) + ".txt";
                }
                if(!opt.output.targetsFile.empty() &&
                    opt.output.targetsFile != opt.output.queryMappingsFile)
                {
                    targetsFile = opt.output.targetsFile + "_"
                            + extract_filename(f) + ".txt";
                }
                if(!opt.output.abundanceFile.empty() &&
                    opt.output.abundanceFile != opt.output.queryMappingsFile)
                {
                    abundanceFile = opt.output.abundanceFile + "_"
                            + extract_filename(f) + ".txt";
                }
                process_input_files(vector<string>{f}, db, opt,
                                    queryMappingsFile, targetsFile, abundanceFile);
            }
        }
    }
    //process all input files at once
    else {
        process_input_files(infiles, db, opt,
                            opt.output.queryMappingsFile,
                            opt.output.targetsFile,
                            opt.output.abundanceFile);
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
        db.read_values_only(filename);
        cerr << "done again." << endl;
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
    if((opt.sketchlen > 0) || (opt.winlen > 0) || (opt.winstride > 0)) {
        auto s = db.query_sketcher();
        if(opt.sketchlen > 0) s.sketch_size(opt.sketchlen);
        if(opt.winstride > 0) s.window_size(opt.winlen);
        if(opt.winlen    > 0) s.window_stride(opt.winstride);
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

    if(opt.output.showDBproperties) {
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

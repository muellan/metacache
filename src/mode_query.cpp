/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2020 André Müller (muellan@uni-mainz.de)
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


#include "mode_query.h"

#include "classification.h"
#include "classification_statistics.h"
#include "cmdline_utility.h"
#include "config.h"
#include "filesys_utility.h"
#include "options.h"
#include "printing.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>


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

    std::ofstream targetMappingsFile;
    if(!targetsFilename.empty()) {
        targetMappingsFile.open(targetsFilename, std::ios::out);

        if(targetMappingsFile.good()) {
            cout << "Per-Target mappings will be written to file: " << targetsFilename << endl;
            perTargetOut = &targetMappingsFile;
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
void process_input_files(const database& db,
                         const query_options& opt)
{
    const auto& infiles = opt.infiles;

    if(infiles.empty()) {
        cerr << "No input filenames provided!\n";
        return;
    }
    else {
        bool noneReadable = std::none_of(infiles.begin(), infiles.end(),
                           [](const auto& f) { return file_readable(f); });
        if(noneReadable) {
            string msg = "None of the following query sequence files could be opened:";
            for(const auto& f : infiles) { msg += "\n    " + f; }
            throw std::runtime_error{std::move(msg)};
        }
    }

    //process files / file pairs separately
    const auto& ano = opt.output.analysis;

    if(opt.splitOutputPerInput) {
        string queryMappingsFile;
        string targetMappingsFile;
        string abundanceFile;
        //process each input file pair separately
        if(opt.pairing == pairing_mode::files && infiles.size() > 1) {
            for(std::size_t i = 0; i < infiles.size(); i += 2) {
                const auto& f1 = infiles[i];
                const auto& f2 = infiles[i+1];
                if(!opt.queryMappingsFile.empty()) {
                    queryMappingsFile = opt.queryMappingsFile
                            + "_" + extract_filename(f1)
                            + "_" + extract_filename(f2)
                            + ".txt";
                }
                if(!ano.targetMappingsFile.empty() &&
                    ano.targetMappingsFile != opt.queryMappingsFile)
                {
                    targetMappingsFile = ano.targetMappingsFile
                            + "_" + extract_filename(f1)
                            + "_" + extract_filename(f2)
                            + ".txt";
                }
                if(!ano.abundanceFile.empty() &&
                    ano.abundanceFile != opt.queryMappingsFile)
                {
                    abundanceFile = ano.abundanceFile
                            + "_" + extract_filename(f1)
                            + "_" + extract_filename(f2)
                            + ".txt";
                }
                process_input_files(vector<string>{f1,f2}, db, opt,
                    queryMappingsFile, targetMappingsFile, abundanceFile);
            }
        }
        //process each input file separately
        else {
            for(const auto& f : infiles) {
                if(!opt.queryMappingsFile.empty()) {
                    queryMappingsFile = opt.queryMappingsFile + "_"
                            + extract_filename(f) + ".txt";
                }
                if(!ano.targetMappingsFile.empty() &&
                    ano.targetMappingsFile != opt.queryMappingsFile)
                {
                    targetMappingsFile = ano.targetMappingsFile + "_"
                            + extract_filename(f) + ".txt";
                }
                if(!ano.abundanceFile.empty() &&
                    ano.abundanceFile != opt.queryMappingsFile)
                {
                    abundanceFile = ano.abundanceFile + "_"
                            + extract_filename(f) + ".txt";
                }
                process_input_files(vector<string>{f}, db, opt,
                    queryMappingsFile, targetMappingsFile, abundanceFile);
            }
        }
    }
    //process all input files at once
    else {
        process_input_files(infiles, db, opt,
                            opt.queryMappingsFile,
                            ano.targetMappingsFile,
                            ano.abundanceFile);
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
void run_interactive_query_mode(const database& db,
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
            vector<string> args {initOpt.dbfile};
            std::istringstream iss(input);
            while(iss >> input) { args.push_back(input); }

            //read command line options (use initial ones as defaults)
            try {
                auto opt = get_query_options(args, initOpt);
                adapt_options_to_database(opt.classify, db);
                process_input_files(db, opt);
            }
            catch(std::exception& e) {
                if(initOpt.output.showErrors) cerr << e.what() << '\n';
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
read_database(const string& filename,
              const database_storage_options& dbopt,
              const sketching_options& skopt)
{
    database db;

    if(dbopt.maxLoadFactor > 0.4 && dbopt.maxLoadFactor < 0.99) {
        db.max_load_factor(dbopt.maxLoadFactor);
        cerr << "Using custom hash table load factor of "
             << dbopt.maxLoadFactor << '\n';
    }

    try {
        db.read(filename, dbopt.numParts);
    }
    catch(const file_access_error& e) {
        cerr << "FAIL" << endl;
        throw;
    }

    if(dbopt.removeOverpopulatedFeatures) {
        auto old = db.feature_count();

        auto maxlpf = dbopt.maxLocationsPerFeature - 1;
        if(maxlpf < 0 || maxlpf >= database::max_supported_locations_per_feature())
            maxlpf = database::max_supported_locations_per_feature() - 1;

        maxlpf = std::min(maxlpf, db.max_locations_per_feature() - 1);
        if(maxlpf > 0) { //always keep buckets with size 1
            cerr << "\nRemoving features with more than "
                 << maxlpf << " locations... " << std::flush;

            auto rem = db.remove_features_with_more_locations_than(maxlpf);

            cerr << rem << " of " << old << " removed.\n";
        }
        //in case new max is less than the database setting
        db.max_locations_per_feature(dbopt.maxLocationsPerFeature);
    }
    else if(dbopt.maxLocationsPerFeature > 1) {
        db.max_locations_per_feature(dbopt.maxLocationsPerFeature);
        cerr << "Max locations per feature set to "
             << dbopt.maxLocationsPerFeature << '\n';
    }

    //use a different sketching scheme for querying?
    if((skopt.sketchlen > 0) || (skopt.winlen > 0) || (skopt.winstride > 0)) {
        auto s = db.query_sketcher();
        if(skopt.sketchlen > 0) s.sketch_size(skopt.sketchlen);
        if(skopt.winlen > 0)    s.window_size(skopt.winlen);
        // if no custom window stride requested => set to w-k+1
        auto winstride = skopt.winstride;
        if(winstride < 0) winstride = s.window_size() - s.kmer_size() + 1;
        s.window_stride(winstride);

        cerr << "custom query sketching settings:"
             << " -winlen "    << s.window_size()
             << " -winstride " << s.window_stride()
             << " -sketchlen " << s.sketch_size() << '\n';

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
void main_mode_query(const cmdline_args& args)
{
    auto opt = get_query_options(args);

    auto db = read_database(opt.dbfile, opt.dbconfig, opt.sketching);

    if(!opt.infiles.empty()) {
        cerr << "Classifying query sequences.\n";

        adapt_options_to_database(opt.classify, db);
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

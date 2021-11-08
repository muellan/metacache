/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2021 André Müller (muellan@uni-mainz.de)
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


#include "querying.h"

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
        const size_t stride = (opt.pairing == pairing_mode::files) &&
                              (infiles.size() > 1) ? 2 : 1;

        for(std::size_t i = 0; i < infiles.size(); i += stride) {
            string suffix;
            vector<string> input;

            if(stride == 2) {
                //process each input file pair separately
                const auto& f1 = infiles[i];
                const auto& f2 = infiles[i+1];
                suffix = "_" + extract_filename(f1)
                       + "_" + extract_filename(f2)
                       + ".txt";
                input = vector<string>{f1,f2};
            }
            else {
                //process each input file separately
                const auto& f = infiles[i];
                suffix = "_" + extract_filename(f) + ".txt";
                input = vector<string>{f};
            }

            if(!opt.queryMappingsFile.empty()) {
                queryMappingsFile = opt.queryMappingsFile + suffix;
            }
            if(!ano.targetMappingsFile.empty() &&
                ano.targetMappingsFile != opt.queryMappingsFile)
            {
                targetMappingsFile = ano.targetMappingsFile + suffix;
            }
            if(!ano.abundanceFile.empty() &&
                ano.abundanceFile != opt.queryMappingsFile)
            {
                abundanceFile = ano.abundanceFile + suffix;
            }
            process_input_files(input, db, opt,
                queryMappingsFile, targetMappingsFile, abundanceFile);
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
void adapt_options_to_database(query_options& opt, const database& db)
{
    sketching_opt& skopt = opt.sketching;

    //use sketching scheme from database?
    const auto& dbsk = db.target_sketching();

    skopt.kmerlen = dbsk.kmerlen;
    if(skopt.sketchlen < 1)
        skopt.sketchlen = dbsk.sketchlen;
    if(skopt.winlen < 1)
        skopt.winlen = dbsk.winlen;
    // if no custom window stride requested => set to w-k+1
    if(skopt.winstride < 1)
        skopt.winstride = skopt.winlen - skopt.kmerlen + 1;

    if((skopt.sketchlen != dbsk.sketchlen) ||
       (skopt.winlen    != dbsk.winlen)    ||
       (skopt.winstride != dbsk.winstride))
    {
        cerr << "custom query sketching settings:"
             << " -sketchlen " << skopt.sketchlen
             << " -winlen "    << skopt.winlen
             << " -winstride " << skopt.winstride
             << '\n';
    }

    classification_options& clopt = opt.classify;

    //deduce hit threshold from database?
    if(clopt.hitsMin < 1) {
        auto sks = db.target_sketching().sketchlen;
        if(sks >= 6) {
            clopt.hitsMin = static_cast<int>(sks / 3.0);
        } else if (sks >= 4) {
            clopt.hitsMin = 2;
        } else {
            clopt.hitsMin = 1;
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
    cout << "Running in interactive mode:\n"
            " - Enter input file name(s) and command line options and press return.\n"
            " - The initially given command line options will be used as defaults.\n"
            " - All command line options that would modify the database are ignored.\n"
            " - Each line will be processed separately.\n"
            " - Lines starting with '#' will be ignored.\n"
            " - Enter an empty line or press Ctrl-D to quit MetaCache.\n"
            << endl;

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
                adapt_options_to_database(opt, db);
                process_input_files(db, opt);
            }
            catch(std::exception& e) {
                if(initOpt.output.showErrors) cerr << e.what() << '\n';
            }
        }
    }
}


} // namespace mc

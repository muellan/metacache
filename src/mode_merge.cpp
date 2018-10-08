/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2018 André Müller (muellan@uni-mainz.de)
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
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include "filesys_utility.h"
#include "args_handling.h"
#include "config.h"
#include "query_options.h"
#include "taxonomy_io.h"
#include "classification.h"
#include "io_error.h"
#include "candidates.h"
#include "printing.h"


namespace mc {

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;


/*************************************************************************//**
 *
 * @brief merge parameters
 *
 *****************************************************************************/
struct merge_options
{
    taxonomy_options taxonomy;

    info_level infoLevel = info_level::moderate;

	query_options query;
};


/*************************************************************************//**
 *
 * @brief command line args -> merge options
 *
 *****************************************************************************/
merge_options
get_merge_options(const args_parser& args)
{
    merge_options opt;

    if(args.contains("silent")) {
        opt.infoLevel = info_level::silent;
    } else if(args.contains("verbose")) {
        opt.infoLevel = info_level::verbose;
    }

    opt.taxonomy = get_taxonomy_options(args);

    opt.query = get_query_options(args, {});

    return opt;
}



/*************************************************************************//**
 *
 * @brief forward stream past the next occurence of a specific character
 *
 *****************************************************************************/
inline void
forward(std::istream& is, char c)
{
    is.ignore(std::numeric_limits<std::streamsize>::max(), c);
}



/*************************************************************************//**
 *
 * @brief merge classification result files
 *
 *****************************************************************************/
void merge_result_files(const vector<string>& infiles,
                        const database& db,
						const query_options& opt,
                        classification_results& results) {
    vector<string> queryNames;
    //get query headers

    vector<classification_candidates> queryCandidates;

    candidate_generation_rules rules;

    rules.mergeBelow    = opt.classify.lowestRank;
    rules.maxCandidates = opt.classify.maxNumCandidatesPerQuery;

    bool getHeader = true;

    for(const auto& file : infiles) {
        std::ifstream ifs(file);
        string line;

        //check classification rank
        while(ifs.good()) {
            getline(ifs, line);
            if(line[0] != '#') {
                throw io_format_error("classificaion ranks not found in file "+file);
            }
            if(line.compare(0,16,"# Classification") == 0) {
                auto pos = line.find("sequence");
                if(pos != string::npos)
                    throw io_format_error("cannot merge results on sequence level");
                break;
            }
        }

        //get layout
        int colTop = 0;
        while(ifs.good()) {
            getline(ifs, line);
            if(line[0] != '#') {
                throw io_format_error("TABLE_LAYOUT not found in file "+file);
            }
            if(line.compare(0,15,"# TABLE_LAYOUT:") == 0) {
                // cout << line.c_str() << endl;
                std::stringstream lineStream(line.substr(15));
                string column;
                lineStream >> column;
                // cout << column << endl;
                if(column != "query_id") {
                    throw io_format_error("no query_id in file "+file);
                }
                int col = 0;
                while(lineStream.good()) {
                    forward(lineStream, '|');
                    lineStream >> column;
                    ++col;
                    if(column == "top_hits") {
                        colTop = col;
                        break;
                    }
                }
                break;
            }
        }
        if(!colTop)
            throw io_format_error("no top_hits in file "+file);
        // cout << "colTop " << colTop << endl;

        char lineBegin = ifs.peek();
        //skip comments
        while(ifs.good() && lineBegin == '#') {
            forward(ifs, '\n');
            lineBegin = ifs.peek();
        }
        //skip results
        auto resultsBegin = ifs.tellg();
        while(ifs.good() && lineBegin != '#') {
            forward(ifs, '\n');
            lineBegin = ifs.peek();
        }
        auto resultsEnd = ifs.tellg();

        //get number of queries
        size_t numQueries = 0;
        while(ifs.good() && numQueries == 0) {
            if(lineBegin == '#') {
                getline(ifs, line);
                if(line.compare(0,10,"# queries:") == 0) {
                    // cout << line.c_str() << endl;
                    std::stringstream lineStream(line.substr(10));
                    lineStream >> numQueries;
                    // cout << numQueries << endl;

                    queryNames.resize(numQueries);
                    queryCandidates.resize(numQueries);
                }
            }
            else {
                forward(ifs, '\n');
            }
            lineBegin = ifs.peek();
        }

        //read results
        //TODO: in parallel
        ifs.seekg(resultsBegin);
        lineBegin = ifs.peek();
        while(ifs.good() && lineBegin != '#') {
            size_t queryId;
            ifs >> queryId;
            cout << queryId << endl;
            // if(queryNames.size() <= queryId) {
            //     queryNames.resize(queryId+1);
            //     queryCandidates.resize(queryId+1);
            // }
            forward(ifs, '|');
            if(getHeader) {
                string header;
                ifs >> header;
                cout << header << endl;
                queryNames[queryId] = std::move(header);
            }

            //skip to tophits
            for(int i = 1; i < colTop; ++i)
                forward(ifs, '|');
            forward(ifs, '\t');

            //get tophits
            char nextChar = ifs.peek();
            while(nextChar != '\t') {
                taxon_id taxid;
                ifs >> taxid;
                cout << taxid << endl;
                const taxon* tax = db.taxon_with_id(taxid);
                forward(ifs, ':');
                match_candidate::count_type hits;
                ifs >> hits;
                cout << hits << endl;

                queryCandidates[queryId].insert(match_candidate{tax, hits}, db, rules);

                nextChar = ifs.get();
                //if(nextChar == ',') continue;
            }
            //end of tophits

            forward(ifs, '\n');
            lineBegin = ifs.peek();

            string dummy;
            std::cin >> dummy;
        }

        getHeader = false;
    }

    //for all queries
    //classify
    //output
    map_candidates_to_targets(queryNames, queryCandidates, db, opt, results);
}


/*************************************************************************//**
 *
 * @brief runs merge on input files; sets output target streams
 *
 *****************************************************************************/
void process_result_files(const vector<string>& infiles,
                         const database& db, const query_options& opt,
                         const string& queryMappingsFilename,
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
    merge_result_files(infiles, db, opt, results);
    results.time.stop();

    clear_current_line(results.status);
    results.status.flush();

    if(opt.output.showSummary) show_summary(opt, results);

    results.flush_all_streams();
}


/*************************************************************************//**
 *
 * @brief runs merge on input files;
 *        handles output file split
 *
 *****************************************************************************/
void process_result_files(const vector<string>& infiles,
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

    //TODO: group input by pre/suffix


    //process files / file pairs separately
    if(opt.output.splitFiles) {
        string queryMappingsFile;
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
                if(!opt.output.abundanceFile.empty() && 
                    opt.output.abundanceFile != opt.output.queryMappingsFile)
                {
                    abundanceFile = opt.output.abundanceFile
                            + "_" + extract_filename(f1)
                            + "_" + extract_filename(f2)
                            + ".txt";
                }
                process_result_files(vector<string>{f1,f2}, db, opt,
                                    queryMappingsFile, abundanceFile);
            }
        }
    }
    //process all input files at once
    else {
        process_result_files(infiles, db, opt,
                            opt.output.queryMappingsFile,
                            opt.output.abundanceFile);
    }
}


/*************************************************************************//**
 *
 * @brief merge classification result files
 *
 *****************************************************************************/
void main_mode_merge(const args_parser& args)
{
    const auto nargs = args.non_prefixed_count();
    if(nargs < 2) {
    	cerr << "no result filenames provided. abort" << endl;
    	return;
    }
    auto infiles = result_filenames(args);
    //TODO: sort & group files by prefix or suffix

    merge_options opt = get_merge_options(args);

    database db;

    if(!opt.taxonomy.path.empty()) {
        db.reset_taxa_above_sequence_level(
            make_taxonomic_hierarchy(opt.taxonomy.nodesFile,
                                     opt.taxonomy.namesFile,
                                     opt.taxonomy.mergeFile,
                                     opt.infoLevel) );

        if(opt.infoLevel != info_level::silent) {
            cout << "Taxonomy applied to database." << endl;
        }
    }

    //TODO update ranks cache at species level

    if(!infiles.empty()) {
        cerr << "Merging result files." << endl;

        process_result_files(infiles, db, opt.query, "", "");
    }
    else {
        cout << "No input files provided. Abort." << endl;
        //interactive mode?
    }
}

} // namespace mc

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
#include <string>
#include <vector>

#include "args_handling.h"
#include "config.h"
#include "query_options.h"
#include "taxonomy_io.h"


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
 * @brief merge classification result files
 *
 *****************************************************************************/
void merge_result_files(const vector<string>& infiles,
						const taxonomy& tax,
						const merge_options& opt) {

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

	taxonomy taxa = make_taxonomic_hierarchy(opt.taxonomy.nodesFile,
                                     opt.taxonomy.namesFile,
                                     opt.taxonomy.mergeFile,
                                     opt.infoLevel);

	ranked_lineages_cache ranksCache{taxa;


    if(!infiles.empty()) {
        cerr << "Merging result files." << endl;

        merge_result_files(infiles, taxa,  opt);
    }
    else {
        cout << "No input files provided. Abort." << endl;
        //interactive mode?
    }
}

} // namespace mc

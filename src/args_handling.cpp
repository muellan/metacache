/*****************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
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

#include <stdexcept>

#include "filesys_utility.h"
#include "args_handling.h"


namespace mc {


//-------------------------------------------------------------------
std::string
database_name(const args_parser& args)
{
    //first 2 non-prefixed args are mode and database name
    if(args.non_prefixed_count() > 1) {
        auto filename = args.non_prefixed(1);
        if(filename.find(".db") == std::string::npos) {
            filename += ".db";
        }
        if(filename != "") return filename;
    }
    throw std::invalid_argument{"No database filename provided"};
}



//-------------------------------------------------------------------
std::vector<std::string>
sequence_filenames(const args_parser& args)
{
    //first 2 non-prefixed args are mode and database name
    //all other non-prefixed args that are not
    //preceded by options ("-xyz") should be sequence file or folder names
    auto n = args.non_prefixed_count();
    if(n < 3) {
        throw std::invalid_argument{"No input sequence filenames provided"};
    }

    auto files = std::vector<std::string>{};
    files.reserve(n-2);
    for(std::size_t i = 2; i < n; ++i) {
        if(args.is_preceded_by_prefixed_arg(i)) break;
        auto name = args.non_prefixed(i);

        auto fnames = files_in_directory(name);

        if(!fnames.empty()) {
            files.insert(files.end(), std::make_move_iterator(fnames.begin()),
                                      std::make_move_iterator(fnames.end()));
        } else {
            files.push_back(name);
        }
    }

    return files;
}



//-------------------------------------------------------------------
taxonomy_param
get_taxonomy_param(const args_parser& args)
{
    taxonomy_param param;

    param.path = args.get<std::string>("taxonomy", std::string(""));
    if(!param.path.empty() &&
        param.path.back() != '/') param.path += '/';

    param.nodesFile = param.path + "nodes.dmp";
    param.namesFile = param.path + "names.dmp";
    param.mergeFile = param.path + "merged.dmp";

    param.mappingPreFiles.push_back("assembly_summary.txt");

    //manually added accession to taxon map file names
    auto postm = args.get<std::string>("taxpostmap", "");
    if(!postm.empty()) {
        param.mappingPostFiles.push_back(postm);
    }

    //default NCBI accession to taxon map file names
    param.mappingPostFiles.push_back(param.path + "nucl_gb.accession2taxid");
    param.mappingPostFiles.push_back(param.path + "nucl_wgs.accession2taxid");
    param.mappingPostFiles.push_back(param.path + "nucl_est.accession2taxid");
    param.mappingPostFiles.push_back(param.path + "nucl_gss.accession2taxid");

    //find additional maps by file extension ".accession2taxid"
    for(const auto f : files_in_directory(param.path)) {
        if(f.find(".accession2taxid") != std::string::npos) {
            param.mappingPostFiles.push_back(f);
        }
    }

    return param;
}



} // namespace mc

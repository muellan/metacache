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

#ifndef MC_ARGS_HANDLING_H_
#define MC_ARGS_HANDLING_H_


#include <string>
#include <vector>

#include "args_parser.h"
#include "taxonomy.h"


namespace mc {



/*****************************************************************************
 *
 * @brief database configuration parameters
 *
 *****************************************************************************/
struct taxonomy_param
{
    std::string path;
    std::string nodesFile;
    std::string namesFile;
    std::string mergeFile;
    std::vector<std::string> mappingPreFiles;
    std::vector<std::string> mappingPostFiles;
};




/*****************************************************************************
 *
 * @brief extracts database name from command line arguments
 *
 *****************************************************************************/
std::string database_name(const args_parser&);



/*****************************************************************************
 *
 * @brief extracts input sequence filenames from command line arguments
 *
 *****************************************************************************/
std::vector<std::string> sequence_filenames(const args_parser&);



/*****************************************************************************
 *
 * @brief extracts database config parameters from cmd line args
 *
 *****************************************************************************/
taxonomy_param get_taxonomy_param(const args_parser&);



} // namespace mc


#endif

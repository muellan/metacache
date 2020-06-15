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

#ifndef MC_TAXONOMY_TOOLS_H_
#define MC_TAXONOMY_TOOLS_H_


#include <vector>
#include <map>
#include <string>

#include "taxonomy.h"
#include "io_options.h"


namespace mc {


/*************************************************************************//**
 *
 * @brief reads taxonomic tree + names from NCBI's taxnonomy files
 *
 *****************************************************************************/
taxonomy
make_taxonomic_hierarchy(const std::string& taxNodesFile,
                         const std::string& taxNamesFile = "",
                         const std::string& mergeTaxFile = "",
                         info_level info = info_level::moderate);



/*************************************************************************//**
 *
 * @brief Reads a text file that maps sequences (ids or filenames) to taxon ids.
 *     Originally intended for use with the NCBI RefSeq "assembly_summary.txt"
 *     files but can be used with any text file provided that:
 *       - All lines in the file have the identical number of
 *         tab separated columns.
 *       - Either the first column contains the key (sequence id or file name).
 *       - Either the first lines contain the token "taxid" in the
 *         column which holds the taxids or the file has only two columns.
 *
 *****************************************************************************/
void read_sequence_to_taxon_id_mapping(
    const std::string& mappingFile,
    std::map<std::string,taxonomy::taxon_id>&,
    info_level info = info_level::moderate);



/*************************************************************************//**
 *
 * @brief reads file -> taxon id mappings files from the same directories
 *        as the input files
 *
 *****************************************************************************/
std::map<std::string,taxonomy::taxon_id>
make_sequence_to_taxon_id_map(const std::vector<std::string>& mappingFilenames,
                              const std::vector<std::string>& globalMappingFilenames,
                              const std::vector<std::string>& infilenames,
                              info_level infoLvl);


} // namespace mc


#endif

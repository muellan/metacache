/*****************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * version 1.0
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

#ifndef MC_MODE_BUILD_H_
#define MC_MODE_BUILD_H_

#include "args_parser.h"
#include "sketch_database.h"
#include "min_hasher.h"


namespace mc {


/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
using sequence_type = std::string;
using sketcher_type = single_function_min_hasher;
//using sketcher_type = multi_function_min_hasher;
using database_type = sketch_database<sequence_type,sketcher_type>;
using taxon_rank    = database_type::taxon_rank;
using genome_id     = database_type::genome_id;



/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
enum class sequence_id_type {
    custom, acc, acc_ver, gi
};



/*****************************************************************************
 *
 * @brief prints sequence to property mapping
 *
 *****************************************************************************/
void main_mode_annotate(const args_parser&);



/*****************************************************************************
 *
 * @brief builds a database from reference input sequences
 *
 *****************************************************************************/
void main_mode_build(const args_parser&);



/*****************************************************************************
 *
 * @brief adds reference sequences to an existing database
 *
 *****************************************************************************/
void main_mode_build_add(const args_parser&);



/*****************************************************************************
 *
 * @brief generates random reads (with errors) from sequences
 *
 *****************************************************************************/
void main_mode_help(const args_parser&);
void show_help_for_topic(const std::string&);



/*****************************************************************************
 *
 * @brief shows database properties
 *
 *****************************************************************************/
void main_mode_info(const args_parser&);



/*****************************************************************************
 *
 * @brief generates random reads (with errors) from sequences
 *
 *****************************************************************************/
void main_mode_make_reads(const args_parser&);



/*****************************************************************************
 *
 * @brief run query reads against pre-built database
 *
 *****************************************************************************/
void main_mode_query(const args_parser&);


} // namespace mc


#endif

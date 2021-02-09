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

#ifndef MC_QUERYING_H_
#define MC_QUERYING_H_


#include "database.h"
#include "options.h"


namespace mc {



/*************************************************************************//**
 *
 * @brief runs classification on input files;
 *        handles output file split
 *
 *****************************************************************************/
void process_input_files(const database& db,
                         const query_options& opt);



/*************************************************************************//**
 *
 * @brief sets up some query options according to database parameters
 *        or command line options
 *
 *****************************************************************************/
void adapt_options_to_database(classification_options& opt, const database& db);



/*************************************************************************//**
 *
 * @brief primitive REPL mode for repeated querying using the same database
 *
 *****************************************************************************/
void run_interactive_query_mode(const database& db,
                                const query_options& initOpt);


} // namespace mc

#endif

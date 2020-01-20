/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2020 André Müller (muellan@uni-mainz.de)
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

#ifndef MC_MODES_H_
#define MC_MODES_H_


#include "cmdline_utility.h"


namespace mc {


/*************************************************************************//**
 *
 * @brief builds a database from reference input sequences (= targets)
 *
 *****************************************************************************/
void main_mode_build(const cmdline_args&);



/*************************************************************************//**
 *
 * @brief adds reference sequences (= targets) to an existing database
 *
 *****************************************************************************/
void main_mode_modify(const cmdline_args&);



/*************************************************************************//**
 *
 * @brief run query reads against pre-built database
 *
 *****************************************************************************/
void main_mode_query(const cmdline_args&);



/*************************************************************************//**
 *
 * @brief merge classification result files
 *
 *****************************************************************************/
void main_mode_merge(const cmdline_args&);



/*************************************************************************//**
 *
 * @brief shows database properties
 *
 *****************************************************************************/
void main_mode_info(const cmdline_args&);



/*************************************************************************//**
 *
 * @brief help
 *
 *****************************************************************************/
void main_mode_help(const cmdline_args&);


} // namespace mc


#endif

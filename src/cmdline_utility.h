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

#ifndef MC_CMDLINE_TOOLS_H_
#define MC_CMDLINE_TOOLS_H_

#include <iosfwd>
#include <vector>
#include <string>


namespace mc {


using cmdline_args = std::vector<std::string>;



/*************************************************************************//**
 *
 * @brief make args C-array into vector of strings
 *
 *****************************************************************************/
cmdline_args make_args_list(char** first, char** last);



/*************************************************************************//**
 *
 * @brief prints a progress indicator like this:  [====>   ] 50%
 *
 *****************************************************************************/
void show_progress_indicator(std::ostream&, float done, int totalLength = 80);

void clear_current_line(std::ostream&, int length = 80);



} // namespace mc


#endif

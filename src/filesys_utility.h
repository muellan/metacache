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

#ifndef MC_FS_TOOLS_H_
#define MC_FS_TOOLS_H_


#include <string>
#include <vector>
#include <set>
#include <fstream>


namespace mc {


/*************************************************************************//**
 *
 * @brief
 * @details this is the only POSIX dependency (working with gcc or MinGW)
 *
 *****************************************************************************/
std::vector<std::string>
files_in_directory(std::string parentDir, int recurse = 10);



/*************************************************************************//**
 *
 * @brief
 *
 *****************************************************************************/
std::set<std::string>
unique_directories(const std::vector<std::string>& filenames);



/*************************************************************************//**
 *
 * @brief
 *
 *****************************************************************************/
std::string
extract_filename(const std::string& filepath);



/*************************************************************************//**
 *
 * @brief returns the size of a file in bytes
 *
 *****************************************************************************/
std::ifstream::pos_type file_size(const std::string& filename);



/*************************************************************************//**
 *
 * @return true, if file with name 'filename' could be opened for reading
 *
 *****************************************************************************/
bool file_readable(const std::string& filename);


} // namespace mc


#endif

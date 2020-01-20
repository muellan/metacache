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
#include <dirent.h> //POSIX header
#include <cstring>
#include <iterator>

#include "filesys_utility.h"


namespace mc {


//-------------------------------------------------------------------
std::vector<std::string>
files_in_directory(std::string parentDir, int recurse)
{
    if(!parentDir.empty() &&
       (parentDir.back() == '/' || parentDir.back() == '\\'))
    {
        parentDir.pop_back();
    }

    auto files = std::vector<std::string>{};
    DIR* dir;
    struct dirent *entry;
    dir = opendir(parentDir.c_str());
    if(dir) {
        while((entry = readdir(dir))) {
            if(strcmp(entry->d_name, ".") != 0 &&
               strcmp(entry->d_name, "..") != 0)
            {
                bool single = true;

                if(recurse > 0) {
                    auto ff = files_in_directory(
                                  parentDir + "/" + std::string(entry->d_name),
                                  recurse-1);

                    if(!ff.empty()) {
                        files.insert(files.end(),
                                     std::make_move_iterator(ff.begin()),
                                     std::make_move_iterator(ff.end()));
                        single = false;
                    }
                }

                if(single) {
                    files.push_back(parentDir + "/" + std::string(entry->d_name));
                }
            }
        }
        closedir(dir);
    }
    return files;
}



//-------------------------------------------------------------------
std::set<std::string>
unique_directories(const std::vector<std::string>& filenames)
{
    auto dirs = std::set<std::string>{};

    for(const auto& filename : filenames) {
        auto i = filename.find_last_of("/\\");
        dirs.insert(filename.substr(0,i));
    }

    return dirs;
}



//-------------------------------------------------------------------
std::string
extract_filename(const std::string& filepath)
{
    if(filepath.empty()) return std::string("");

    auto i = filepath.find_last_of("/\\");
    return filepath.substr(i+1);
}



//-------------------------------------------------------------------
std::ifstream::pos_type
file_size(const std::string& filename)
{
    std::ifstream is{filename, std::ifstream::ate | std::ifstream::binary};
    if(!is.good()) return 0;
    return is.tellg();
}



//-------------------------------------------------------------------
bool file_readable(const std::string& filename)
{
    return std::ifstream{filename}.good();
}


} // namespace mc


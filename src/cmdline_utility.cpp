/*****************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2017 André Müller (muellan@uni-mainz.de)
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

#include "cmdline_utility.h"


namespace mc {


//-------------------------------------------------------------------
void show_progress_indicator(float done, int totalLength)
{
    auto m = int((totalLength - 7) * done);
    std::cout << "\r[";
    for(int j = 0; j < m; ++j) std::cout << '=';
    std::cout << ">";
    m = totalLength - 7 - m;
    for(int j = 0; j < m; ++j) std::cout << ' ';
    std::cout << "] " << int(100 * done) << "%" << std::flush;
}


//-------------------------------------------------------------------
void clear_current_line(int length)
{
    std::cout << '\r';
    for(; length > 0; --length) std::cout << ' ';
    std::cout << '\r';
}


} // namespace mc

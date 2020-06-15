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

#include <iostream>

#include "cmdline_utility.h"


namespace mc {


//-------------------------------------------------------------------
cmdline_args make_args_list(char** first, char** last)
{
    cmdline_args args;
    if(last < first) return args;

    args.reserve(last-first);

    for(; first != last; ++first) {
        args.push_back(*first);
    }

    return args;
}



//-------------------------------------------------------------------
void show_progress_indicator(std::ostream& os, float done, int totalLength)
{
    if(done < 0.f) done = 0.f;
    auto m = int((totalLength - 7) * done);
    os << "\r[";
    for(int j = 0; j < m; ++j) os << '=';
    os << ">";
    m = totalLength - 7 - m;
    for(int j = 0; j < m; ++j) os << ' ';
    os << "] " << int(100 * done) << "%" << std::flush;
}



//-------------------------------------------------------------------
void clear_current_line(std::ostream& os, int length)
{
    os << '\r';
    for(; length > 0; --length) os << ' ';
    os << '\r' << std::flush;
}


} // namespace mc

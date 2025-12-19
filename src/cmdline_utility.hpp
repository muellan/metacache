/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2026 André Müller (github.com/muellan)
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

#ifndef MC_CMDLINE_TOOLS_HPP_
#define MC_CMDLINE_TOOLS_HPP_


#include <atomic>
#include <future>
#include <iosfwd>
#include <string>
#include <vector>


namespace mc {


using cmdline_args = std::vector<std::string>;



//-----------------------------------------------------------------------------
/**
 * @brief make args C-array into vector of strings
 */
cmdline_args make_args_list(char** first, char** last);



//-----------------------------------------------------------------------------
/**
 * @brief prints a progress indicator like this:  [====>   ] 50%
 */
void show_progress_indicator (std::ostream&, float done, int totalLength = 80);

void clear_current_line (std::ostream&, int length = 80);



//-----------------------------------------------------------------------------
/**
 * @brief show progress on single thread, updateable from multiple threads
 */
struct concurrent_progress
{
    std::atomic_size_t counter{0};
    std::atomic_size_t total{0};
    bool initialized{false};

    float progress () const noexcept {
        return std::min(1.0f, float(counter)/total);
    }

    void show (std::ostream& os) {
        initialized = true;
        if (total > 0) {
            show_progress_indicator(os, progress());
        }
        else {
            show_progress_indicator(os, 0);
        }
    }

    void clear_line (std::ostream& os) const {
        if (initialized) clear_current_line(os);
    }
};



//-----------------------------------------------------------------------------
/**
 * @brief show progress on single thread until all futures are ready
 */
void show_progress_until_ready (std::ostream&, concurrent_progress&,
                                std::vector<std::future<void>>&);

} // namespace mc

#endif

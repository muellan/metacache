/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2021 Robin Kobus  (kobus@uni-mainz.de)
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

#ifndef MC_QUERY_BATCH_H_
#define MC_QUERY_BATCH_H_




namespace mc {


//-------------------------------------------------------------------
/** @brief used for query result storage/accumulation
 */
template<class Location>
class matches_sorter {

public:
    using location = Location;
    using match_locations = std::vector<location>;

    void sort() {
        merge_sort(locs_, offsets_, temp_);
    }

    void clear() {
        locs_.clear();
        offsets_.clear();
    }

    void next() {
        offsets_.clear();
        offsets_.emplace_back(locs_.size());
    }

    bool empty() const noexcept { return locs_.empty(); }
    auto size()  const noexcept { return locs_.size(); }

    auto begin() const noexcept { return locs_.begin(); }
    auto end()   const noexcept { return locs_.end(); }

    const match_locations&
    locations() const noexcept { return locs_; }

    template<class InputIt>
    void append(InputIt first, InputIt last)
    {
        locs_.insert(locs_.end(), first, last);
        offsets_.emplace_back(locs_.size());
    }

private:
    static void
    merge_sort(match_locations& inout,
               std::vector<size_t>& offsets,
               match_locations& temp)
    {
        if(offsets.size() < 3) return;
        temp.resize(inout.size());

        int numChunks = offsets.size()-1;
        for(int s = 1; s < numChunks; s *= 2) {
            for(int i = 0; i < numChunks; i += 2*s) {
                auto begin = offsets[i];
                auto mid = i + s <= numChunks ? offsets[i + s] : offsets[numChunks];
                auto end = i + 2*s <= numChunks ? offsets[i + 2*s] : offsets[numChunks];
                std::merge(inout.begin()+begin, inout.begin()+mid,
                            inout.begin()+mid, inout.begin()+end,
                            temp.begin()+begin);
            }
            std::swap(inout, temp);
        }
        if(numChunks % 2) {
            std::copy(inout.begin()+offsets.front(),
                      inout.begin()+offsets.back(),
                      temp.begin()+offsets.front());
            std::swap(inout, temp);
        }
    }

    match_locations locs_; // match locations from hashmap
    std::vector<std::size_t> offsets_;  // bucket sizes for merge sort
    match_locations temp_; // temp buffer for merge sort
};


//-------------------------------------------------------------------
/** @brief used for query result storage/accumulation
 */
template<class Location>
struct query_handler {

    using location = Location;
    using sorter = matches_sorter<location>;
    using match_locations = typename sorter::match_locations;

    explicit query_handler(const sketcher& querySketcher_) :
        querySketcher{querySketcher_},
        matchesSorter{}
    {}

    sketcher querySketcher;
    sorter matchesSorter;
};


}

#endif

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

#ifndef MC_CANDIDATES_H_
#define MC_CANDIDATES_H_


#include "config.h"


namespace mc {


/*****************************************************************************
 *
 * @brief
 *
 *****************************************************************************/
template<class ValueT>
struct index_range
{
    using value_type = ValueT;

    constexpr
    index_range() noexcept :
        beg(0), end(0)
    {}
    constexpr
    index_range(value_type first, value_type last) noexcept :
        beg{first}, end{last}
    {}

    constexpr value_type size()  const noexcept { return beg - end; }
    constexpr value_type empty() const noexcept { return beg == end; }

    value_type beg;
    value_type end;
};




/*****************************************************************************
*
* @brief processes a database hit list and
*        stores the top (in terms of accumulated hits) 'maxNo' contiguous
*        window ranges of *distinct* targets
*
*****************************************************************************/
template<int maxNo, class Database>
class top_matches_in_contiguous_window_range
{
    static_assert(maxNo > 1, "no must be > 1");


public:

    //---------------------------------------------------------------
    using matches_per_location = typename Database::matches_per_location;
    using target_id     = typename Database::target_id;
    using window_id     = typename Database::window_id;
    using hit_count     = std::uint_least64_t;
    using window_range  = index_range<window_id>;

    //---------------------------------------------------------------
    static constexpr int max_count() noexcept { return maxNo; }

    static constexpr target_id invalid_target_id() noexcept {
        return Database::invalid_target_id();
    }


    /****************************************************************
     * @pre matches must be sorted by target (first) and window (second)
     */
    top_matches_in_contiguous_window_range(
        const matches_per_location& matches, window_id numWindows = 3)
    :
        tgt_{}, hits_{}, pos_{}
    {
        using std::begin;
        using std::end;

        for(int i = 0; i < maxNo; ++i) {
            tgt_[i] = invalid_target_id();
            hits_[i] = 0;
        }

        target_id tgt = invalid_target_id();
        window_id win = 0;
        window_id maxWinBeg = 0;
        window_id maxWinEnd = 0;
        hit_count hits = 0;
        hit_count maxHits = 0;

        //check hits per query sequence
        auto fst = begin(matches);
        auto lst = fst;
        while(lst != end(matches)) {
            //look for neighboring windows with the highest total hit count
            //as long as we are in the same target and the windows are in a
            //contiguous range
            if(lst->first.tgt == tgt) {
                //add new hits to the right
                hits += lst->second;
                //subtract hits to the left that fall out of range
                while(fst != lst &&
                     (lst->first.win - fst->first.win) >= numWindows)
                {
                    hits -= fst->second;
                    //move left side of range
                    ++fst;
                    win = fst->first.win;
                }
                //track best of the local sub-ranges
                if(hits > maxHits) {
                    maxHits = hits;
                    maxWinBeg = win;
                    maxWinEnd = win + distance(fst,lst);
                }
            }
            else {
                //reset to new target
                win = lst->first.win;
                tgt = lst->first.tgt;
                hits = lst->second;
                maxHits = hits;
                maxWinBeg = win;
                maxWinEnd = win;
                fst = lst;
            }
            //keep track of 'maxNo' largest hits for *distinct* targets
            for(int i = 0; i < maxNo; ++i) {
                if(maxHits >= hits_[i]) {
                    //check if tgt already in list
                    int pos = maxNo-1;
                    for(int j = i; j < maxNo; ++j) {
                        if(tgt == tgt_[j]) {
                            pos = j;
                            break;
                        }
                    }
                    //shift targets left of tgt to the right
                    // for(int j = maxNo-1; j > i; --j) {
                    for(int j = pos; j > i; --j) {
                        hits_[j] = hits_[j-1];
                        tgt_[j] = tgt_[j-1];
                        pos_[j] = pos_[j-1];
                    }
                    //set hits & associated sequence (position)
                    hits_[i] = maxHits;
                    tgt_[i] = tgt;
                    pos_[i].beg = maxWinBeg;
                    pos_[i].end = maxWinEnd;
                    break;
                }
            }
            ++lst;
        }
    }

    int count() const noexcept {
        for(int i = 0; i < maxNo; ++i) {
            if(hits_[i] < 1) return i;
        }
        return maxNo;
    }

    target_id target(int rank) const noexcept { return tgt_[rank];  }
    hit_count hits(int rank)   const noexcept { return hits_[rank]; }

    hit_count total_hits() const noexcept {
        hit_count h = 0;
        for(int i = 0; i < maxNo; ++i) h += hits_[i];
        return h;
    }

    const window_range&
    window(int rank) const noexcept { return pos_[rank]; }

    window_id window_length(int rank) const noexcept {
        return pos_[rank].end - pos_[rank].beg;
    }

private:
    target_id tgt_[maxNo];
    hit_count hits_[maxNo];
    window_range pos_[maxNo];
};



} // namespace mc


#endif

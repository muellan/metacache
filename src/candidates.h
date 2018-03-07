/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2018 André Müller (muellan@uni-mainz.de)
 *                       & Robin Kobus  (rkobus@uni-mainz.de)
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


#include <cstdint>
#include <limits>

#include "config.h"


namespace mc {


/*************************************************************************//**
 *
 * @brief inclusive index range: [beg,end]
 *
 *****************************************************************************/
struct window_range
{
    using window_id = mc::window_id;

    constexpr
    window_range() noexcept = default;

    constexpr
    window_range(window_id first, window_id last) noexcept :
        beg{first}, end{last}
    {}

    constexpr window_id size()  const noexcept { return end - beg; }
    constexpr window_id empty() const noexcept { return beg == end; }

    window_id beg = 0;
    window_id end = 0;
};



/*************************************************************************//**
 *
 * @brief hit count and position in candidate taxon
 *
 *****************************************************************************/
struct match_candidate
{
    using taxon        = mc::taxon;
    using window_id    = mc::window_id;
    using window_range = mc::window_range;
    using count_type   = std::uint_least64_t;

    constexpr
    match_candidate() noexcept = default;

    const taxon* tax = nullptr;
    count_type   hits = 0;
    window_range pos;
};



/*************************************************************************//**
 *
 * @brief produces all contiguous window ranges of matches
 *        that are at most 'numWindows' long
 *
 * @pre   matches must be sorted by taxon (first) and window (second)
 *
 *****************************************************************************/
template<class Consumer>
void for_all_contiguous_window_ranges(const matches_per_location& matches,
                                      window_id numWindows,
                                      Consumer&& consume)
{
    using hit_count = match_candidate::count_type;

    using std::begin;
    using std::end;

    auto fst = begin(matches);

    //list empty?
    if(fst == end(matches)) return;

    //first entry in list
    hit_count hits = fst->hits;
    match_candidate curBest;
    curBest.tax = fst->loc.tax;
    curBest.hits = hits;
    curBest.pos.beg = fst->loc.win;;
    curBest.pos.end = fst->loc.win;;
    auto lst = fst;
    ++lst;

    //rest of list: check hits per query sequence
    while(lst != end(matches)) {
        //look for neighboring windows with the highest total hit count
        //as long as we are in the same target and the windows are in a
        //contiguous range
        if(lst->loc.tax == curBest.tax) {
            //add new hits to the right
            hits += lst->hits;
            //subtract hits to the left that fall out of range
            while(fst != lst &&
                (lst->loc.win - fst->loc.win) >= numWindows)
            {
                hits -= fst->hits;
                //move left side of range
                ++fst;
            }
            //track best of the local sub-ranges
            if(hits > curBest.hits) {
                curBest.hits = hits;
                curBest.pos.beg = fst->loc.win;
                curBest.pos.end = lst->loc.win;
            }
        }
        else { //end of current target
            consume(curBest);
            //reset to new target
            fst = lst;
            hits = fst->hits;
            curBest.tax  = fst->loc.tax;
            curBest.hits = hits;
            curBest.pos.beg = fst->loc.win;
            curBest.pos.end = fst->loc.win;
        }

        ++lst;
    }
    consume(curBest);
}



/*************************************************************************//**
*
* @brief processes a database hit list and
*        stores the top (in terms of accumulated hits) 'maxNo' contiguous
*        window ranges of *distinct* targets
*
*****************************************************************************/
template<int maxNo>
class top_distinct_matches_in_contiguous_window_range
{
    static_assert(maxNo > 1, "no must be > 1");


public:
    using const_iterator = const match_candidate*;


    /****************************************************************
     * @pre matches must be sorted by taxon (first) and window (second)
     */
    top_distinct_matches_in_contiguous_window_range(
        const database& db,
        const matches_per_location& matches,
        //max. allowed number of windows in a contiguous range
        window_id numWindows = 3,
        //list only the best candidate of a taxon on rank 'mergeOn'
        taxon_rank mergeOn = taxon_rank::Sequence)
    :
        top_{}
    {
        for_all_contiguous_window_ranges(matches, numWindows,
            [&db,mergeOn,this] (match_candidate cand)
        {
            if(!cand.tax) return;

            if(mergeOn > taxon_rank::Sequence) {
                auto ancestor = db.ancestor(cand.tax, mergeOn);
                if(ancestor) cand.tax = ancestor;
            }

            if(cand.tax->rank() == taxon_rank::Sequence) {
                //note: sequence-level duplicates are not possible
                for(int i = 0; i < maxNo; ++i) {
                    if(cand.hits >= top_[i].hits) {
                        //shift targets left of tax to the right
                        for(int j = maxNo-1; j > i; --j) top_[j] = top_[j-1];
                        top_[i] = cand;
                        return;
                    }
                }
            }
            else {
                for(int i = 0; i < maxNo; ++i) {
                    if(cand.hits >= top_[i].hits) {
                        //check if tgt already in list
                        int pos = maxNo-1;
                        for(int j = i; j < maxNo; ++j) {
                            if(cand.tax == top_[j].tax) {
                                pos = j;
                                break;
                            }
                        }
                        //shift targets left of tgt to the right
                        for(int j = pos; j > i; --j) top_[j] = top_[j-1];
                        top_[i] = cand;
                        break;
                    }
                    else if(cand.tax == top_[i].tax) {
                        break;
                    }
                }
            }
        });
    }


    //---------------------------------------------------------------
    const_iterator begin() const noexcept {
        using std::begin;
        return begin(top_);
    }
    const_iterator end() const noexcept {
        return begin() + size();
    }

    bool empty() const noexcept { return !top_[0].tax; }

    int size() const noexcept {
        for(int i = 0; i < maxNo; ++i) if(!top_[i].tax) return i;
        return maxNo;
    }

    const match_candidate&
    operator [] (int rank) const noexcept {
        return top_[rank];
    }


private:
    match_candidate top_[maxNo];
};




/*************************************************************************//**
 *
 * @brief
 *
 *****************************************************************************/
class all_distinct_matches_in_contiguous_window_range
{
    using candidates_list = std::vector<match_candidate>;

public:
    using const_iterator = candidates_list::const_iterator;


    /****************************************************************
     * @pre matches must be sorted by taxon (first) and window (second)
     */
    all_distinct_matches_in_contiguous_window_range(
        const database&, //not needed here
        const matches_per_location& matches,
        //max. allowed number of windows in a contiguous range
        window_id numWindows = 3,
        taxon_rank = taxon_rank::Sequence) // will be ignored
    :
        top_{}
    {
        for_all_contiguous_window_ranges(matches, numWindows,
            [this] (const match_candidate& cand) {
                top_.push_back(cand);
            });
    }


    //---------------------------------------------------------------
    const_iterator begin() const noexcept { return top_.begin(); }
    const_iterator end()   const noexcept { return top_.end(); }

    bool empty() const noexcept { return top_.empty(); }
    int  size()  const noexcept { return int(top_.size()); }

    const match_candidate&
    operator [] (int rank) const noexcept { return top_[rank]; }


private:
    candidates_list top_;
};


} // namespace mc


#endif

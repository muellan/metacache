/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2019 André Müller (muellan@uni-mainz.de)
 *                       & Robin Kobus  (kobus@uni-mainz.de)
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

#ifndef MC_CANDIDATE_GENERATION_H_
#define MC_CANDIDATE_GENERATION_H_


#include <cstdint>
#include <limits>

#include "config.h"
#include "taxonomy.h"

#include "../dep/cudahelpers/cuda_helpers.cuh"

namespace mc {


template<class Value>
struct span {
    using value_type = Value;

    value_type * begin() const noexcept { return begin_; }

    value_type * end() const noexcept { return end_; }

    size_t size() const noexcept { return end_ - begin_; }

    bool empty() const noexcept { return begin_ == end_; }

    value_type * begin_ = nullptr;
    value_type * end_   = nullptr;
};


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

    constexpr window_id size()  const noexcept { return end - beg + 1; }

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
    using count_type   = std::uint32_t;

    // constexpr
    match_candidate() noexcept = default;

    match_candidate(const taxon* tax, count_type hits) :
        tax{tax},
        tgt{std::numeric_limits<target_id>::max()},
        hits{hits},
        pos{}
    {};

    // const taxon* tax = nullptr;
    const taxon* tax;
    // target_id    tgt = std::numeric_limits<target_id>::max();
    target_id    tgt;
    // count_type   hits = 0;
    count_type   hits;
    window_range pos;
};



/*************************************************************************//**
 *
 * @brief candidate generation parameters
 *
 *****************************************************************************/
struct candidate_generation_rules
{
    constexpr candidate_generation_rules() noexcept = default;

    //maximum length of contiguous window range
    window_id maxWindowsInRange = 3;

    //maximum number of candidates to be generated
    std::size_t maxCandidates = std::numeric_limits<std::size_t>::max();

    //list only the best candidate of a taxon on rank
    taxon_rank mergeBelow = taxon_rank::Sequence;
};



/*************************************************************************//**
 *
 * @brief  produces all contiguous window ranges of matches
 *         that are at most 'numWindows' long;
 *         the main loop is aborted if 'consume' returns false
 *
 * @pre    matches must be sorted by taxon (first) and window (second)
 *
 * @tparam Consumer : function object that processes one match_candidate
 *                    and returns a bool indicating, if loop should be continued
 *
 *****************************************************************************/
#ifdef __CUDACC__
    #pragma hd_warning_disable
#endif
template<class InputIterator, class Consumer>
HOSTDEVICEQUALIFIER
void for_all_contiguous_window_ranges(
    InputIterator fst,
    InputIterator end,
    window_id numWindows,
    Consumer&& consume)
{
    using hit_count = match_candidate::count_type;

    //list empty?
    if(fst == end) return;

    //first entry in list
    hit_count hits = 1;
    match_candidate curBest;
    curBest.tgt     = fst->tgt;
    curBest.hits    = hits;
    curBest.pos.beg = fst->win;
    curBest.pos.end = fst->win;
    auto lst = fst;
    ++lst;

    //rest of list: check hits per query sequence
    while(lst != end) {
        //look for neighboring windows with the highest total hit count
        //as long as we are in the same target and the windows are in a
        //contiguous range
        if(lst->tgt == curBest.tgt) {
            //add new hits to the right
            hits++;
            //subtract hits to the left that fall out of range
            while(fst != lst &&
                (lst->win - fst->win) >= numWindows)
            {
                hits--;
                //move left side of range
                ++fst;
            }
            //track best of the local sub-ranges
            if(hits > curBest.hits) {
                curBest.hits    = hits;
                curBest.pos.beg = fst->win;
                curBest.pos.end = lst->win;
            }
        }
        else { //end of current target
            if(!consume(curBest)) return;
            //reset to new target
            fst = lst;
            hits = 1;
            curBest.tgt     = fst->tgt;
            curBest.hits    = hits;
            curBest.pos.beg = fst->win;
            curBest.pos.end = fst->win;
        }

        ++lst;
    }
    if(!consume(curBest)) return;
}


/*************************************************************************//**
 *
 * @brief  produces all contiguous window ranges of matches
 *         that are at most 'numWindows' long;
 *         the main loop is aborted if 'consume' returns false
 *
 * @pre    matches must be sorted by taxon (first) and window (second)
 *
 * @tparam Consumer : function object that processes one match_candidate
 *                    and returns a bool indicating, if loop should be continued
 *
 *****************************************************************************/
template<class Locations, class Consumer>
void for_all_contiguous_window_ranges(
    const Locations& matches,
    window_id numWindows,
    Consumer&& consume)
{
    using std::begin;
    using std::end;

    for_all_contiguous_window_ranges(
        begin(matches), end(matches), numWindows,
        std::forward<Consumer>(consume));
}

} // namespace mc


#endif
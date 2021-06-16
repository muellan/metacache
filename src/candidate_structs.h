/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2021 André Müller (muellan@uni-mainz.de)
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

#ifndef MC_CANDIDATE_STRUCTS_H_
#define MC_CANDIDATE_STRUCTS_H_


#include "config.h"
#include "options.h"
#include "taxonomy.h"

#include <cstdint>
#include <limits>


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

    friend bool
    operator < (const window_range& a, const window_range& b) noexcept {
        if(a.beg < b.beg) return true;
        if(a.beg > b.beg) return false;
        return (a.end < b.end);
    }
    friend bool
    operator > (const window_range& a, const window_range& b) noexcept {
        if(a.beg > b.beg) return true;
        if(a.beg < b.beg) return false;
        return (a.end > b.end);
    }

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
 * @brief make rules for candidate generation
 *
 *****************************************************************************/
template<class Query>
candidate_generation_rules
make_candidate_generation_rules(const sketcher& targetSketcher,
                                const classification_options& opt,
                                const Query& query)
{
    candidate_generation_rules rules;

    rules.maxWindowsInRange = window_id( 2 + (
        std::max(query.seq1.size() + query.seq2.size(), opt.insertSizeMax) /
        targetSketcher.window_stride() ));

    rules.mergeBelow    = opt.lowestRank;
    rules.maxCandidates = opt.maxNumCandidatesPerQuery;

    return rules;
}



} // namespace mc


#endif
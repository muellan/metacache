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
            if(!consume(curBest)) return;
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
    if(!consume(curBest)) return;
}



/*************************************************************************//**
*
* @brief processes a database match list and
*        stores contiguous window ranges of *distinct* targets
*        as a list *sorted* by accumulated hits
*
*****************************************************************************/
class best_distinct_matches_in_contiguous_window_ranges
{
    using candidates_list = std::vector<match_candidate>;

public:
    using size_type      = std::size_t;
    using const_iterator = candidates_list::const_iterator;


    /****************************************************************
     * @pre matches must be sorted by taxon (first) and window (second)
     */
    best_distinct_matches_in_contiguous_window_ranges(
        const database& db,
        const matches_per_location& matches,
        const candidate_generation_rules& rules = candidate_generation_rules{})
    :
        top_{}
    {
        using hit_count = match_candidate::count_type;

        for_all_contiguous_window_ranges(matches, rules.maxWindowsInRange,
            [&,this] (match_candidate cand)
        {
            if(!cand.tax) return true;

            if(rules.mergeBelow > taxon_rank::Sequence) {
                auto ancestor = db.ancestor(cand.tax, rules.mergeBelow);
                if(ancestor) cand.tax = ancestor;
            }

            if(cand.tax->rank() == taxon_rank::Sequence) {
                auto i = std::lower_bound(top_.begin(), top_.end(), cand.hits,
                    [&] (const match_candidate& c, hit_count h) {
                        return c.hits > h;
                    });

                if(i != top_.end() || top_.size() < rules.maxCandidates) {
                    top_.insert(i, cand);
                    if(top_.size() > rules.maxCandidates)
                        top_.resize(rules.maxCandidates);
                }
            }
            //above sequence level, taxa can occur more than once
            else {
                auto i = std::find_if(top_.begin(), top_.end(),
                    [&] (const match_candidate& c) {
                        return c.tax == cand.tax;
                    });

                if(i != top_.end()) {
                    //taxon already in list, update, if more hits
                    if(cand.hits > i->hits) *i = cand;
                    std::sort(top_.begin(), i+1,
                        [] (const match_candidate& a, const match_candidate& b) {
                            return a.hits > b.hits;
                        });
                }
                //taxon not in list yet
                else  {
                    auto j = std::lower_bound(top_.begin(), top_.end(), cand.hits,
                        [&] (const match_candidate& c, hit_count h) {
                            return c.hits > h;
                        });

                    if(j != top_.end() || top_.size() < rules.maxCandidates) {
                        top_.insert(j, cand);
                        if(top_.size() > rules.maxCandidates)
                            top_.resize(rules.maxCandidates);
                    }
                }

            }
            return true;
        });
    }


    //---------------------------------------------------------------
    const_iterator begin() const noexcept { return top_.begin(); }
    const_iterator end()   const noexcept { return top_.end(); }

    bool empty()       const noexcept { return top_.empty(); }
    size_type size() const noexcept { return top_.size(); }

    const match_candidate&
    operator [] (size_type i) const noexcept { return top_[i]; }


private:
    candidates_list top_;
};




/*************************************************************************//**
*
* @brief processes a database match list and
*        stores contiguous window ranges of *distinct* targets
*
*****************************************************************************/
class distinct_matches_in_contiguous_window_ranges
{
    using candidates_list = std::vector<match_candidate>;

public:
    using size_type      = std::size_t;
    using const_iterator = candidates_list::const_iterator;


    /****************************************************************
     * @pre matches must be sorted by taxon (first) and window (second)
     */
    distinct_matches_in_contiguous_window_ranges(
        const database&, //not needed here
        const matches_per_location& matches,
        const candidate_generation_rules& rules = candidate_generation_rules{})
    :
        top_{}
    {
        for_all_contiguous_window_ranges(matches, rules.maxWindowsInRange,
            [&,this] (const match_candidate& cand) {
                if(top_.size() >= rules.maxWindowsInRange) return false;
                top_.push_back(cand);
                return true;
            });
    }


    //---------------------------------------------------------------
    const_iterator begin() const noexcept { return top_.begin(); }
    const_iterator end()   const noexcept { return top_.end(); }

    bool empty() const noexcept { return top_.empty(); }

    size_type size()  const noexcept { return top_.size(); }

    const match_candidate&
    operator [] (size_type i) const noexcept { return top_[i]; }


private:
    candidates_list top_;
};


} // namespace mc


#endif

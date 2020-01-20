/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2020 André Müller (muellan@uni-mainz.de)
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

#ifndef MC_CANDIDATES_H_
#define MC_CANDIDATES_H_


#include <cstdint>
#include <limits>

#include "sketch_database.h"


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
    using count_type   = std::uint_least64_t;

    constexpr
    match_candidate() noexcept = default;

    match_candidate(const taxon* tax, count_type hits) :
        tax{tax},
        tgt{std::numeric_limits<target_id>::max()},
        hits{hits},
        pos{}
    {};

    const taxon* tax = nullptr;
    target_id    tgt = std::numeric_limits<target_id>::max();
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
template<class Locations, class Consumer>
void for_all_contiguous_window_ranges(
    const Locations& matches,
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
    hit_count hits = 1;
    match_candidate curBest;
    curBest.tgt     = fst->tgt;
    curBest.hits    = hits;
    curBest.pos.beg = fst->win;
    curBest.pos.end = fst->win;
    auto lst = fst;
    ++lst;

    //rest of list: check hits per query sequence
    while(lst != end(matches)) {
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
* @brief processes a database match list and
*        stores contiguous window ranges of *distinct* targets
*        as a list *sorted* by accumulated hits
*
*****************************************************************************/
class best_distinct_matches_in_contiguous_window_ranges
{
    using candidates_list  = std::vector<match_candidate>;

public:
    using size_type      = std::size_t;
    using iterator       = candidates_list::iterator;
    using const_iterator = candidates_list::const_iterator;


    /****************************************************************
     */
    best_distinct_matches_in_contiguous_window_ranges() = default;

    /****************************************************************
     * @pre matches must be sorted by taxon (first) and window (second)
     */
    template<class Locations>
    best_distinct_matches_in_contiguous_window_ranges(
        const database& db,
        const Locations& matches,
        const candidate_generation_rules& rules = candidate_generation_rules{})
    :
        top_{}
    {
        for_all_contiguous_window_ranges(matches, rules.maxWindowsInRange,
            [&,this] (match_candidate& cand) {
                return insert(cand, db, rules);
            });
    }


    //---------------------------------------------------------------
    auto begin() const noexcept { return top_.begin(); }
    auto end()   const noexcept { return top_.end(); }

    bool empty()     const noexcept { return top_.empty(); }
    size_type size() const noexcept { return top_.size(); }

    const match_candidate&
    operator [] (size_type i) const noexcept { return top_[i]; }

    iterator erase(const_iterator pos) { return top_.erase(pos); }

    /****************************************************************
     * @brief insert candidate and keep list sorted
     */
    bool insert(match_candidate cand,
                const database& db,
                const candidate_generation_rules& rules = candidate_generation_rules{})
    {
        if(rules.mergeBelow > taxon_rank::Sequence)
            cand.tax = db.lowest_ranked_ancestor(cand.tgt, rules.mergeBelow);
        else
            cand.tax = db.taxon_of_target(cand.tgt);

        if(!cand.tax) return true;

        auto greater = [] (const match_candidate& a, const match_candidate& b) {
                           return a.hits > b.hits;
                       };

        if(cand.tax->rank() == taxon_rank::Sequence) {
            auto i = std::upper_bound(top_.begin(), top_.end(), cand, greater);

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
                if(cand.hits > i->hits) {
                    *i = cand;
                    std::sort(top_.begin(), i+1, greater);
                }
            }
            //taxon not in list yet
            else  {
                auto j = std::upper_bound(top_.begin(), top_.end(), cand, greater);

                if(j != top_.end() || top_.size() < rules.maxCandidates) {
                    top_.insert(j, cand);

                    if(top_.size() > rules.maxCandidates)
                        top_.resize(rules.maxCandidates);
                }
            }

        }
        return true;
    };


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
    using candidates_list  = std::vector<match_candidate>;

public:
    using size_type      = std::size_t;
    using const_iterator = candidates_list::const_iterator;


    distinct_matches_in_contiguous_window_ranges() = default;


    /****************************************************************
     * @pre matches must be sorted by taxon (first) and window (second)
     */
    distinct_matches_in_contiguous_window_ranges(
        const database& db,
        const match_locations& matches,
        const candidate_generation_rules& rules = candidate_generation_rules{})
    :
        cand_{}
    {
        for_all_contiguous_window_ranges(matches, rules.maxWindowsInRange,
            [&,this] (match_candidate& cand) {
                return insert(cand, db, rules);
            });
    }


    /****************************************************************
     * @brief insert candidate and keep list sorted
     */
    bool insert(match_candidate cand,
                const database& db,
                const candidate_generation_rules& = candidate_generation_rules{})
    {
        cand.tax = db.taxon_of_target(cand.tgt);
        cand_.push_back(cand);
        return true;
    }


    //---------------------------------------------------------------
    auto begin() const noexcept { return cand_.begin(); }
    auto end()   const noexcept { return cand_.end(); }

    bool empty() const noexcept { return cand_.empty(); }

    size_type size()  const noexcept { return cand_.size(); }

    const match_candidate&
    operator [] (size_type i) const noexcept { return cand_[i]; }

private:
    candidates_list cand_;
};


} // namespace mc


#endif

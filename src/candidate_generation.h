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

#ifndef MC_CANDIDATE_GENERATION_H_
#define MC_CANDIDATE_GENERATION_H_


#include "candidate_structs.h"
#include "options.h"

#include <cstdint>
#include <limits>


namespace mc {


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
template<class InputIterator, class Consumer>
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
    curBest.tax     = nullptr;
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
            curBest.tax     = nullptr;
            curBest.tgt     = fst->tgt;
            curBest.hits    = hits;
            curBest.pos.beg = fst->win;
            curBest.pos.end = fst->win;
        }

        ++lst;
    }
    if(!consume(curBest)) return;
}


//-------------------------------------------------------------------
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
     * @brief copy candidates from span
     */
    void assign(const span<const match_candidate> cand) {
        top_.assign(cand.begin(), cand.end());
    }


    /****************************************************************
     * @pre matches must be sorted by taxon (first) and window (second)
     */
    template<class Locations>
    void insert(
        const taxonomy_cache& taxonomy,
        const Locations& matches,
        const candidate_generation_rules& rules = candidate_generation_rules{})
    {
        if(rules.maxCandidates < std::numeric_limits<std::size_t>::max())
            top_.reserve(rules.maxCandidates+1);

        for_all_contiguous_window_ranges(matches, rules.maxWindowsInRange,
            [&,this] (match_candidate& cand) {
                return insert(cand, taxonomy, rules);
            });
    }


    /****************************************************************
     * @brief insert candidate and keep list sorted
     */
    bool insert(match_candidate cand,
                const taxonomy_cache& taxonomy,
                const candidate_generation_rules& rules = candidate_generation_rules{})
    {
        // early exit
        if(top_.size() == rules.maxCandidates && top_.back().hits >= cand.hits) return true;

        if(!cand.tax) {
            if(rules.mergeBelow > taxon_rank::Sequence)
                cand.tax = taxonomy.lowest_ranked_ancestor(cand.tgt, rules.mergeBelow);
            else
                cand.tax = taxonomy.cached_taxon_of_target(cand.tgt);
        }

        if(!cand.tax) return true;

        auto greater = [] (const match_candidate& a, const match_candidate& b) {
                           return a.hits > b.hits;
                       };

        if(rules.mergeBelow == taxon_rank::Sequence) {
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
            else {
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


    //---------------------------------------------------------------
    auto begin() const noexcept { return top_.begin(); }
    auto end()   const noexcept { return top_.end(); }

    bool empty()     const noexcept { return top_.empty(); }
    size_type size() const noexcept { return top_.size(); }

    void clear() { top_.clear(); }

    const match_candidate&
    operator [] (size_type i) const noexcept { return top_[i]; }

    iterator erase(const_iterator pos) { return top_.erase(pos); }

    auto view() const noexcept { return span<const match_candidate>(top_); }

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
    using iterator       = candidates_list::iterator;
    using const_iterator = candidates_list::const_iterator;


    /****************************************************************
     * @brief copy candidates from span
     */
    void assign(const span<const match_candidate> cand) {
        cand_.assign(cand.begin(), cand.end());
    }


    /****************************************************************
     * @pre matches must be sorted by taxon (first) and window (second)
     */
    template<class Locations>
    void insert(
        const taxonomy_cache& taxonomy,
        const Locations& matches,
        const candidate_generation_rules& rules = candidate_generation_rules{})
    {
        for_all_contiguous_window_ranges(matches, rules.maxWindowsInRange,
            [&,this] (match_candidate& cand) {
                return insert(cand, taxonomy, rules);
            });
    }


    /****************************************************************
     * @brief insert candidate
     */
    bool insert(match_candidate cand,
                const taxonomy_cache& taxonomy,
                const candidate_generation_rules& = candidate_generation_rules{})
    {
        cand.tax = taxonomy.cached_taxon_of_target(cand.tgt);
        cand_.push_back(cand);
        return true;
    }


    //---------------------------------------------------------------
    auto begin() const noexcept { return cand_.begin(); }
    auto end()   const noexcept { return cand_.end(); }

    bool empty() const noexcept { return cand_.empty(); }

    size_type size()  const noexcept { return cand_.size(); }

    void clear() { cand_.clear(); }

    const match_candidate&
    operator [] (size_type i) const noexcept { return cand_[i]; }

    iterator erase(const_iterator pos) { return cand_.erase(pos); }

    auto view() const noexcept { return span<const match_candidate>(cand_); }

private:
    candidates_list cand_;
};



} // namespace mc


#endif

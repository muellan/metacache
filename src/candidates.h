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

#ifndef MC_CANDIDATES_H_
#define MC_CANDIDATES_H_


#include <cstdint>
#include <limits>

#include "sketch_database.h"
#include "candidate_generation.h"


namespace mc {


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

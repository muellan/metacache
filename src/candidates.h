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


#include <cstdint>
#include <limits>


namespace mc {


/*****************************************************************************
 *
 * @brief
 *
 *****************************************************************************/
template<class WindowId>
struct window_range
{
    using value_type = WindowId;

    constexpr
    window_range() noexcept = default;

    constexpr
    window_range(value_type first, value_type last) noexcept :
        beg{first}, end{last}
    {}

    constexpr value_type size()  const noexcept { return end - beg; }
    constexpr value_type empty() const noexcept { return beg == end; }

    value_type beg = 0;
    value_type end = 0;
};



/*****************************************************************************
 *
 * @brief
 *
 *****************************************************************************/
template<class Database>
struct match_candidate
{
    using taxon        = typename Database::taxon;
    using window_id    = typename Database::window_id;
    using count_type   = std::uint_least64_t;
    using window_range = mc::window_range<window_id>;


    match_candidate() = default;

    const taxon* tax = nullptr;
    count_type   hits = 0;
    window_range pos;
};



/*****************************************************************************
*
* @brief processes a database hit list and
*        stores the top (in terms of accumulated hits) 'maxNo' contiguous
*        window ranges of *distinct* targets
*
*****************************************************************************/
template<class Database, int maxNo>
class top_distinct_matches_in_contiguous_window_range
{
    static_assert(maxNo > 1, "no must be > 1");


public:
    //---------------------------------------------------------------
    using matches_per_location = typename Database::matches_per_location;
    using target_id    = typename Database::target_id;
    using window_id    = typename Database::window_id;
    using taxon        = typename Database::taxon;
    using taxon_rank   = typename Database::taxon_rank;
    using candidate    = match_candidate<Database>;
    using window_range = typename candidate::window_range;
    using hit_count    = typename candidate::count_type;


    //---------------------------------------------------------------
    static constexpr int max_count() noexcept { return maxNo; }


    /****************************************************************
     * @pre matches must be sorted by taxon (first) and window (second)
     */
    top_distinct_matches_in_contiguous_window_range(
        const Database& db,
        const matches_per_location& matches,
        //max. allowed number of windows in a contiguous range
        window_id numWindows = 3,
        //list only the best candidate of a taxon on rank 'mergeOn'
        taxon_rank mergeOn = taxon_rank::Sequence)
    :
        top_{}
    {
        using std::begin;
        using std::end;

        candidate curBest;
        window_id win = 0;
        hit_count hits = 0;

        //check hits per query sequence
        auto fst = begin(matches);
        auto lst = fst;
        while(lst != end(matches)) {
            //look for neighboring windows with the highest total hit count
            //as long as we are in the same target and the windows are in a
            //contiguous range
            if(lst->first.tax == curBest.tax) {
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
                if(hits > curBest.hits) {
                    curBest.hits = hits;
                    curBest.pos.beg = win;
                    curBest.pos.end = win + distance(fst,lst);
                }
            }
            else { //end of current target
                if(curBest.tax) {
                    if(mergeOn > taxon_rank::Sequence) {
                        auto ancestor = db.ancestor(curBest.tax, mergeOn);
                        if(ancestor) curBest.tax = ancestor;
                    }
                    update_with(curBest);
                }

                //reset to new target
                win  = lst->first.win;
                hits = lst->second;
                curBest.tax  = lst->first.tax;
                curBest.hits = hits;
                curBest.pos.beg = win;
                curBest.pos.end = win;
                fst = lst;
            }

            ++lst;
        }
    }

    int size() const noexcept {
        for(int i = 0; i < maxNo; ++i) if(!top_[i].tax) return i;
        return maxNo;
    }

    const candidate& operator [] (int rank) const noexcept {
        return top_[rank];
    }

    hit_count total_hits() const noexcept {
        hit_count h = 0;
        for(int i = 0; i < maxNo; ++i) h += top_[i].hits;
        return h;
    }

private:
    //---------------------------------------------------------------
    /** @brief keeps track of 'maxNo' taxa with largest hits */
    void update_with(const candidate& latest)
    {
        if(latest.tax->rank() == taxon_rank::Sequence) {
            //note: sequence-level duplicates are not possible
            for(int i = 0; i < maxNo; ++i) {
                if(latest.hits >= top_[i].hits) {
                    top_[i] = latest;
                    //shift targets left of tax to the right
                    for(int j = maxNo-1; j > i; --j) top_[j] = top_[j-1];
                    return;
                }
            }
        }
        else {
            for(int i = 0; i < maxNo; ++i) {
                if(top_[i].tax == latest.tax) { //same taxon already in
                    //more hits -> replace
                    if(latest.hits > top_[i].hits) {
                        top_[i] = latest;
                        return;
                    }
                    //less hits -> done
                    if(latest.hits <= top_[i].hits) return;
                }
                //not same taxon and more hits
                else if(latest.hits > top_[i].hits) {
                    top_[i] = latest;
                    //look for duplicates (can only be behind the current one)
                    for(int j = i+1; j < maxNo; ++j) {
                        if(top_[j].tax == latest.tax) {
                            //shift targets right of j to the left
                            for(int k = j+1; k < maxNo; ++k) top_[k-1] = top_[k];
                            //blank out last one
                            top_[maxNo-1].tax = nullptr;
                            top_[maxNo-1].hits = 0;
                            return; //can only happen once
                        }
                    }
                }
            }
        }
    }


    //---------------------------------------------------------------
    candidate top_[maxNo];
};



} // namespace mc


#endif

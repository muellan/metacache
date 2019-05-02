/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2019 André Müller (muellan@uni-mainz.de)
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

#ifndef MC_MATCHES_PER_TARGET_H_
#define MC_MATCHES_PER_TARGET_H_


#include <vector>
#include <unordered_map>
#include <algorithm> // move, sort, lower_bound
#include <iterator>  // make_move_iterator

#include "candidates.h"


namespace mc {


/*************************************************************************//**
 *
 * @brief records matches (and their query origin) per classification target
 *
 *****************************************************************************/
class matches_per_target
{

public:
    //---------------------------------------------------------------
    using window_id        = database::window_id;
    using match_count_type = database::match_count_type;
    using location         = database::match_locations::value_type;
    using taxon            = database::taxon;
    using taxon_rank       = database::taxon_rank;

    struct window_matches {
        window_matches() = default;

        constexpr
        window_matches(window_id w, match_count_type c) noexcept :
            win{w}, hits{c}
        {}

        friend bool
        operator < (const window_matches& a, const window_matches& b) noexcept {
            return (a.win < b.win);
        }
        friend bool
        operator > (const window_matches& a, const window_matches& b) noexcept {
            return (a.win > b.win);
        }

        window_id win = 0;
        match_count_type hits = 0;
    };

    using matches_per_window = std::vector<window_matches>;

    struct candidate {
        candidate(query_id qid, matches_per_window mpw):
            qeryid{qid}, matches{std::move(mpw)}
        {}
        query_id qeryid = 0;
        matches_per_window matches;
    };

private:
    using hits_per_target =
        std::unordered_map<const taxon*, std::vector<candidate>>;

public:
    using const_iterator  = hits_per_target::const_iterator;


    //-------------------------------------------------------------------
    bool empty() const noexcept {
        return hitsPerTarget_.empty();
    }

    std::size_t size() const noexcept {
        return hitsPerTarget_.size();
    }

    const_iterator find(const taxon* tax) const {
        return hitsPerTarget_.find(tax);
    }

    const_iterator begin() const noexcept { return hitsPerTarget_.begin(); }
    const_iterator end()   const noexcept { return hitsPerTarget_.end(); }


    //---------------------------------------------------------------
    void insert(query_id qid,
                const match_locations& matches,
                const classification_candidates& candidates,
                match_count_type minHitsPerCandidate = 0)
    {
        auto lowTaxIter = candidates.begin_lowest_taxa();

        for(const auto& cand : candidates) {

            if(cand.tax && cand.hits >= minHitsPerCandidate) {

                auto tax = cand.tax;
                // try to get sequence-level taxon if neccessary
                if(tax->rank() != taxon_rank::Sequence) tax = *lowTaxIter;

                if(tax->rank() == taxon_rank::Sequence) {
                    // find candidate in matches
                    location lm{tax, cand.pos.beg};
                    auto it = std::lower_bound(matches.begin(), matches.end(), lm,
                        [] (const location& a, const location& b) {
                            //assumes that sequence level taxa have negative taxids
                            if(a.tax->id() < b.tax->id()) return false;
                            if(a.tax->id() > b.tax->id()) return true;
                            return (a.win < b.win);
                        });
                    // fill window vector
                    if(it == matches.end()) return;

                    // create new window vector
                    matches_per_window mpw;
                    mpw.reserve(cand.pos.end - cand.pos.beg + 1);

                    while(it != matches.end() &&
                          it->tax == tax &&
                          it->win <= cand.pos.end)
                    {
                        if(mpw.size() > 0 && mpw.back().win == it->win)
                            mpw.back().hits++;
                        else
                            mpw.emplace_back(it->win, 1);
                        ++it;
                    }
                    // insert into map
                    hitsPerTarget_[tax].emplace_back(qid, std::move(mpw));
                }
            }
        }
        ++lowTaxIter;
    }

    //---------------------------------------------------------------
    void merge(matches_per_target&& other)
    {
        for(auto& mapping : other) {

            auto& source = mapping.second;
            auto& target = hitsPerTarget_[mapping.first];

            target.insert(target.end(),
                          std::make_move_iterator(source.begin()),
                          std::make_move_iterator(source.end()) );
        }
    }

    //---------------------------------------------------------------
    void sort_match_lists()
    {
        for(auto& mapping : hitsPerTarget_) {
            std::sort(mapping.second.begin(), mapping.second.end(),
                [] (const candidate& a, const candidate& b) {
                    if(a.matches.front() < b.matches.front()) return true;
                    if(a.matches.front() > b.matches.front()) return false;
                    if(a.matches.back() < b.matches.back()) return true;
                    if(a.matches.back() > b.matches.back()) return false;
                    return (a.qeryid < b.qeryid);
                });
        }
    }

private:
    hits_per_target hitsPerTarget_;
};


} // namespace mc

#endif

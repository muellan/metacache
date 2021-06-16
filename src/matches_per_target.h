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

#ifndef MC_MATCHES_PER_TARGET_H_
#define MC_MATCHES_PER_TARGET_H_


#include "candidate_generation.h"
#include "database.h"

#include <algorithm> // move, sort, lower_bound
#include <iterator>  // make_move_iterator
#include <unordered_map>
#include <vector>


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
    using window_range     = match_candidate::window_range;
    using count_type       = match_candidate::count_type;

    struct cover_candidate {
        cover_candidate(query_id qid, window_range pos, count_type hits):
            qid{qid}, hits{hits}, pos{pos}
        {}

        friend bool
        operator < (const cover_candidate& a, const cover_candidate& b) noexcept {
            if(a.pos < b.pos) return true;
            if(a.pos > b.pos) return false;
            return (a.qid < b.qid);
        }

        query_id     qid;
        count_type   hits;
        window_range pos;
    };

private:
    using hits_per_target =
        std::unordered_map<target_id, std::vector<cover_candidate>>;

public:
    using iterator        = hits_per_target::iterator;
    using const_iterator  = hits_per_target::const_iterator;


    //-------------------------------------------------------------------
    bool empty() const noexcept {
        return hitsPerTarget_.empty();
    }

    std::size_t size() const noexcept {
        return hitsPerTarget_.size();
    }

    const_iterator find(target_id tgt) const {
        return hitsPerTarget_.find(tgt);
    }

    const_iterator begin() const noexcept { return hitsPerTarget_.begin(); }
    const_iterator end()   const noexcept { return hitsPerTarget_.end(); }

    iterator erase(const_iterator pos) { return hitsPerTarget_.erase(pos); }

    hits_per_target::size_type erase(const hits_per_target::key_type& key) {
        return hitsPerTarget_.erase(key);
    }

    //---------------------------------------------------------------
    void insert(query_id qid,
                const classification_candidates& candidates,
                match_count_type minHitsPerCandidate = 0)
    {
        for(const auto& cand : candidates) {
            // insert valid candidates into map
            if(cand.tax && cand.hits >= minHitsPerCandidate) {
                hitsPerTarget_[cand.tgt].emplace_back(qid, cand.pos, cand.hits);
            }
        }
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
            std::sort(mapping.second.begin(), mapping.second.end());
        }
    }

private:
    hits_per_target hitsPerTarget_;
};


} // namespace mc

#endif

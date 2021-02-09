/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2020 Robin Kobus  (kobus@uni-mainz.de)
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

#ifndef MC_GPU_RESULT_PROCESSING_H_
#define MC_GPU_RESULT_PROCESSING_H_

#include "config.h"
#include "candidate_structs.h"
#include "dna_encoding.h"

#include "../dep/hpc_helpers/include/cuda_helpers.cuh"

namespace mc {


/*************************************************************************//**
 *
 * @brief  remove memory gaps between results of different windows;
 *         determine segment offsets and store them in global memory
 *
 *****************************************************************************/
template<class Id, class Result>
__global__
void compact_results(
    uint32_t numWindows,
    const int * resultOffsets,
    uint32_t maxResultsPerWindow,
    const Result * results_in,
    Result * results_out,
    const Id * segmentIds,
    int * segmentOffsets
)
{
    using id_type = Id;

    for(int bid = blockIdx.x; bid < numWindows; bid += gridDim.x) {
        const int begin = (bid > 0) ? resultOffsets[bid-1] : 0;
        const int end   = resultOffsets[bid];
        const int numResults = end - begin;

        const ulonglong2 * inPtr = reinterpret_cast<const ulonglong2 *>(results_in + bid*maxResultsPerWindow);
        ulonglong2 * outPtr = reinterpret_cast<ulonglong2 *>(results_out + begin);

        for(int tid = threadIdx.x; tid < numResults/2; tid += blockDim.x) {
            outPtr[tid] = inPtr[tid];
        }

        // determine end offset for this query
        if(threadIdx.x == 0) {
            const id_type segmentId = segmentIds[bid];

            // intermediate windows are ignored (marked max)
            // last window of segment sets end offset
            if(segmentId != std::numeric_limits<id_type>::max())
                segmentOffsets[segmentId+1] = end;
        }
    }

    // for(int tid = threadIdx.x + blockIdx.x * blockDim.x; tid < numQueries; tid += blockDim.x * gridDim.x) {
    //     const id_type segmentId = segmentIds[tid];
    //     const int end = resultOffsets[tid];

    //     // intermediate windows are ignored (marked max)
    //     // last query of segment sets end offset
    //     if(segmentId != std::numeric_limits<id_type>::max())
    //         segmentOffsets[segmentId+1] = end;
    // }
}


//-----------------------------------------------------------------------------
__device__
const taxon*
lowest_ranked_ancestor(const ranked_lineage * targetLineages,
                       target_id tgt,
                       taxon_rank lowest)
{
    const auto& lineage = targetLineages[tgt];
    for(int rank = int(lowest); rank < int(taxon_rank::none); ++rank) {
        if(lineage[rank])
            return lineage[rank];
    }
    return nullptr;
}

//-----------------------------------------------------------------------------
__device__
const taxon*
cached_taxon_of_target(const ranked_lineage * targetLineages, target_id tgt)
{
    return targetLineages[tgt][0];
}


/*************************************************************************//**
 * @brief  a location with one or more hits
 *****************************************************************************/
struct match {
    location loc;
    uint32_t hits;
};


/*************************************************************************//**
 *
 * @brief  calculate length of runs of identical locations to create matches
 *
 *****************************************************************************/
template<class Iterator, class Matches>
__device__
uint32_t reduce_locations(
    const Iterator begin, const Iterator end,
    Matches& matches, const uint32_t matchesBeginOffset, const uint32_t numLocations)
{
    using hit_count = match_candidate::count_type;

    const int tid = threadIdx.x;

    const auto myPointer = begin + tid;

    const uint64_t myLocation = myPointer < end ?
        *reinterpret_cast<const uint64_t*>(myPointer) : uint64_t(~0);
    hit_count hits = 1;

    // sum hits of same location run
    for(int i = 1; i < 32; i *= 2) {
        uint64_t otherLocation;
        otherLocation = __shfl_down_sync(0xFFFFFFFF, myLocation, i);
        hit_count otherHits = __shfl_down_sync(0xFFFFFFFF, hits, i);
        if(tid < (32-i) && myLocation == otherLocation)
            hits += otherHits;
    }
    // find last thread of same location run
    uint64_t otherLocation;
    // tid == 0: otherLocation = myLocation
    otherLocation = __shfl_up_sync(0xFFFFFFFF, myLocation, 1);

    // find threads which have to write
    bool predicate = !(myLocation == otherLocation);

    if(tid == 0) {
        // predicate is false because myLocation == otherLocation
        // check if first location is the same as previous last location
        uint32_t pos = (matchesBeginOffset+numLocations-1) % 96;
        location otherLocation = matches.locs[pos];
        if(numLocations > 0 && myLocation == ((uint64_t(otherLocation.tgt) << 32) + otherLocation.win)) {
            matches.hits[pos] += hits;
        }
        else {
            // different location -> thread has to write
            predicate = true;
        }
    }

    predicate = predicate && (myLocation != uint64_t{~0});

    // find write positions
    const uint32_t mask = __ballot_sync(0xFFFFFFFF, predicate);
    const uint32_t count = __popc(mask);
    const uint32_t numPre = __popc(mask & ((1 << tid) -1));
    const uint32_t locPos = (matchesBeginOffset + numLocations + numPre) % 96;

    // store location and hits sum
    if(predicate) {
        matches.locs[locPos] = {uint32_t(myLocation), uint32_t(myLocation >> 32)};
        matches.hits[locPos] = {hits};

        // printf("tid: %d tgt: %d win: %d hits: %d\n",
            // tid, myLocation.tgt, myLocation.win, hits);
    }

    return numLocations + count;
}


/*************************************************************************//**
 *
 * @brief  insert candidate into thread local tophits list
 *
 *****************************************************************************/
 template<int MAX_CANDIDATES>
__device__
bool insert_into_tophits(
    match_candidate& cand, match_candidate (&top)[MAX_CANDIDATES],
    const ranked_lineage * lineages, const taxon_rank mergeBelow)
{
    // early exit
    if(top[MAX_CANDIDATES-1].hits >= cand.hits) return true;

    if(mergeBelow > taxon_rank::Sequence)
        cand.tax = lowest_ranked_ancestor(lineages, cand.tgt, mergeBelow);
    else
        cand.tax = cached_taxon_of_target(lineages, cand.tgt);

    if(!cand.tax) return true;

    int insertPos = MAX_CANDIDATES;
    int taxPos = MAX_CANDIDATES-1;
    // find insert position of cand
    #pragma unroll
    for(int i = 0; i < MAX_CANDIDATES; ++i) {
        if(cand.hits > top[i].hits && i < insertPos) {
            insertPos = i;
        }
    }
    // taxa can occur more than once, but only one candidate per taxon allowed
    // find candidate with same taxon
    #pragma unroll
    for(int i = 0; i < MAX_CANDIDATES; ++i) {
        if(cand.tax == top[i].tax && i < taxPos) {
            taxPos = i;
        }
    }

    // insert except if the same taxon with more hits already exists
    if(taxPos >= insertPos) {
        // move smaller candidates backwards until taxPos
        #pragma unroll
        for(int y = MAX_CANDIDATES-1; y > 0; y--){
            if(y <= taxPos && y > insertPos){
                top[y] = top[y-1];
            }
        }
        // insert cand
        #pragma unroll
        for(int y = 0 ; y < MAX_CANDIDATES; y++){
            if(y == insertPos){
                top[y] = cand;
            }
        }
    }

    return true;
}


/*************************************************************************//**
 *
 * @brief  produces contiguous window ranges of matches
 *         that are at most 'maxWin' long;
 *
 *****************************************************************************/
 template<int MAX_CANDIDATES, class Matches>
 __device__
uint32_t process_matches(
    const Matches& matches, const uint32_t matchesBeginOffset, const uint32_t numLocations,
    const window_id maxWin,
    match_candidate (&top)[MAX_CANDIDATES],
    const ranked_lineage * lineages, const taxon_rank mergeBelow)
{
    using hit_count = match_candidate::count_type;

    const int tid = threadIdx.x;
    const uint32_t numProcessed = min(32, numLocations);

    match_candidate candidate;
    if(tid < numProcessed) {
        location loc = matches.locs[(matchesBeginOffset+tid) % 96];
        uint32_t hits = matches.hits[(matchesBeginOffset+tid) % 96];
        candidate.tgt = loc.tgt;
        candidate.hits = hits;
        candidate.pos.beg = loc.win;
        candidate.pos.end = loc.win;

        for(int i = 1; i < maxWin && (tid+i < numLocations); ++i) {
            location loc = matches.locs[(matchesBeginOffset+tid+i) % 96];
            uint32_t hits = matches.hits[(matchesBeginOffset+tid+i) % 96];
            // check if m in same range
            if(candidate.tgt == loc.tgt && (candidate.pos.beg+maxWin > loc.win)) {
                candidate.hits += hits;
                candidate.pos.end = loc.win;
            }
            else {
                break;
            }
        }
        // printf("tid: %d tgt: %d beg: %d end: %d hits: %d\n",
            // tid, candidate.tgt, candidate.pos.beg, candidate.pos.end, candidate.hits);

        // top[0] = candidate;
    }

    if(tid < numProcessed) {
        // printf("max -- tid: %d tgt: %d beg: %d end: %d hits: %d tid bits: %d\n",
            //    tid, candidate.tgt, candidate.pos.beg, candidate.pos.end, candidate.hits, (maxHits & ((1 << 5) - 1)));

        // save max candidates
        insert_into_tophits<MAX_CANDIDATES>(candidate, top, lineages, mergeBelow);
    }

    return numProcessed;
}


template<
    int MAX_CANDIDATES,
    class Result>
__global__
void generate_top_candidates(
    const uint32_t numSegments,
    const int * segmentOffsets,
    const Result * locations,
    const window_id * maxWindowsInRange,
    const ranked_lineage * lineages,
    const taxon_rank mergeBelow,
    const uint32_t maxCandidatesPerQuery,
    match_candidate * topCandidates,
    const bool updateTopCandidates)
{
    using hit_count = match_candidate::count_type;

    //TODO reduce size
    __shared__ struct {
        location locs[96];
        uint32_t hits[96];
    } matches;
    __shared__ match_candidate topShared[MAX_CANDIDATES];

    for(int bid = blockIdx.x;
            bid < numSegments;
            bid += gridDim.x)
    {
        const int tid = threadIdx.x;

        // if(tid == 0)
        //     printf("block %d offsets %d %d\n", bid, segmentOffsets[bid], segmentOffsets[bid+1]);

        auto begin = locations + segmentOffsets[bid];
        const auto end = locations + segmentOffsets[bid+1];

        match_candidate top[MAX_CANDIDATES];

        #pragma unroll
        for(int i = 0; i < MAX_CANDIDATES; ++i) {
            top[i].hits = 0;
        }

        if(updateTopCandidates && tid < maxCandidatesPerQuery) {
            top[0] = topCandidates[bid*maxCandidatesPerQuery+tid];
            // if(top[0].hits > 0)
            //     printf("bid: %d top: %d tgt: %d hits: %d tax: %llu\n", bid, tid, top[0].tgt, top[0].hits, top[0].tax);
        }

        if(tid < maxCandidatesPerQuery) {
            topShared[tid].tax = nullptr;
            // topShared[tid].tgt = std::numeric_limits<target_id>::max();
            topShared[tid].hits = 0;
            // topShared[tid].pos = {0, 0};
        }

        uint32_t matchesBeginOffset = 0;
        uint32_t numLocations = 0;
        const window_id maxWin = maxWindowsInRange[bid];
        // if(tid == 0) printf("bid: %d maxWin: %d\n", bid, maxWin);

        while(begin < end || numLocations > 0) {
            for(; (begin < end) && (numLocations <= 32+maxWin-1); begin += 32) {
                numLocations = reduce_locations(begin, end, matches, matchesBeginOffset, numLocations);

                // for(int i = tid; i < numLocations; i += 32) {
                //     printf("bid %d tid: %d tgt: %d win: %d hits: %d\n",
                //         bid, tid, matches[tid].loc.tgt, matches[tid].loc.win, matches[tid].hits);
                // }
            }

            // process matches
            const uint32_t numProcessed = process_matches<MAX_CANDIDATES>(
                matches, matchesBeginOffset, (begin < end) ? numLocations-1 : numLocations,
                maxWin, top, lineages, mergeBelow);

            // if(tid == 0)
                // printf("#locs: %d #processed: %d\n", numLocations, numProcessed);

            // remove processed matches from shared memory
            numLocations -= numProcessed;
            matchesBeginOffset += numProcessed;
        }

        // int tophitsCount = 0;
        // #pragma unroll
        // for(int i = 0; i < MAX_CANDIDATES; ++i) {
        //     if(top[i].hits > 0)
        //         ++tophitsCount;
        // }
        // printf("bid: %d tid: %d: tophits: %d\n", bid, tid, tophitsCount);


        int tophitsUsed = 0;
        int tophitsTotal = 0;

        match_candidate candidate = top[0];

        while(__ballot_sync(0xFFFFFFFF, candidate.hits > 0) && tophitsTotal < maxCandidatesPerQuery) {
            // add thread id to make hit count unique
            hit_count maxHits = (candidate.hits << 5) + 32-1 - tid;

            // find top candidates in whole warp
            for(int i = 1; i < 32; i *= 2) {
                hit_count otherHits = __shfl_xor_sync(0xFFFFFFFF, maxHits, i);
                if(maxHits < otherHits)
                    maxHits = otherHits;
            }

            //TODO use whole warp for insertion
            bool insert = false;
            if(32-1 - tid == (maxHits & ((1 << 5) - 1))) {
                insert = true;
                // if same tax already present, do not insert
                for(int i = 0; i < tophitsTotal; ++i) {
                    if(candidate.tax == topShared[i].tax) {
                        insert = false;
                        break;
                    }
                }
                if(insert) {
                    topShared[tophitsTotal] = candidate;
                }
                ++tophitsUsed;

                if(tophitsUsed < MAX_CANDIDATES) {
                    #pragma unroll
                    for(int y = 0; y < MAX_CANDIDATES; y++){
                        if(y == tophitsUsed){
                            candidate = top[y];
                        }
                    }
                }
                else {
                    candidate.hits = 0;
                }
            }
            if(__ballot_sync(0xFFFFFFFF, insert))
                ++tophitsTotal;
        }

        if(tid < maxCandidatesPerQuery) {
            topCandidates[bid*maxCandidatesPerQuery+tid] = topShared[tid];
            // if(topShared[tid].hits > 0)
            //     printf("bid: %d top: %d tgt: %d hits: %d tax: %llu\n", bid, tid, topShared[tid].tgt, topShared[tid].hits, topShared[tid].tax);
        }
    }
}




} // namespace mc


#endif
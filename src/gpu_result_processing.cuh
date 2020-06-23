#ifndef MC_GPU_RESULT_PROCESSING_H_
#define MC_GPU_RESULT_PROCESSING_H_

#include "config.h"
#include "candidate_structs.h"
#include "dna_encoding.h"

#include "../dep/cudahelpers/cuda_helpers.cuh"

namespace mc {


/*************************************************************************//**
 *
 * @brief  remove memory gaps between results of different queries;
 *         determine segment offsets and store them in global memory
 *
 *****************************************************************************/
template<class Id, class Result>
__global__
void compact_kernel(
    uint32_t numQueries,
    const int * resultPrefixScan,
    uint32_t maxResultsPerQuery,
    const Result * results_in,
    Result * results_out,
    const Id * segmentIds,
    int * segmentOffsets
)
{
    using id_type = Id;

    for(int bid = blockIdx.x; bid < numQueries; bid += gridDim.x) {
        const int begin = (bid > 0) ? resultPrefixScan[bid-1] : 0;
        const int end   = resultPrefixScan[bid];
        const int numResults = end - begin;

        const ulonglong2 * inPtr = reinterpret_cast<const ulonglong2 *>(results_in + bid*maxResultsPerQuery);
        ulonglong2 * outPtr = reinterpret_cast<ulonglong2 *>(results_out + begin);

        for(int tid = threadIdx.x; tid < numResults/2; tid += blockDim.x) {
            outPtr[tid] = inPtr[tid];
        }

        // determine end offset for segment of this query
        if(threadIdx.x == 0) {
            const id_type segmentId = segmentIds[bid];
            const id_type nextId    = (bid+1 < numQueries) ? segmentIds[bid+1] : id_type(~0);

            // last query of segment sets end offset
            if(segmentId != nextId)
                segmentOffsets[segmentId+1] = end;
        }
    }

    // for(int tid = threadIdx.x + blockIdx.x * blockDim.x; tid < numQueries; tid += blockDim.x * gridDim.x) {
    //     const id_type segmentId = segmentIds[tid];
    //     const id_type nextId    = (tid+1 < numQueries) ? segmentIds[tid+1] :  id_type(~0);;

    //     const int end = resultPrefixScan[tid];

    //     // last query of segment sets end offset
    //     if(segmentId != nextId)
    //         segmentOffsets[segmentId+1] = end;
    // }
}


/*************************************************************************//**
 *
 * @brief  determine segment offsets and store them in global memory
 *
 *****************************************************************************/
template<class Id>
__global__
void segment_kernel(
    uint32_t numQueries,
    const int * resultPrefixScan,
    const Id * segmentIds,
    int * segmentOffsets
)
{
    using id_type = Id;

    for(int tid = threadIdx.x + blockIdx.x * blockDim.x; tid < numQueries; tid += blockDim.x * gridDim.x) {
        const id_type segmentId = segmentIds[tid];
        const id_type nextId    = (tid+1 < numQueries) ? segmentIds[tid+1] :  id_type(~0);;

        const int end = resultPrefixScan[tid];

        // last query of segment sets end offset
        if(segmentId != nextId)
            segmentOffsets[segmentId+1] = end;
    }
}


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

__device__
const taxon*
taxon_of_target(const ranked_lineage * targetLineages, target_id tgt)
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
template<class Iterator>
__device__
int reduce_locations(Iterator begin, Iterator end, match * matches, int numLocations) {
    using hit_count = match_candidate::count_type;

    const int tid = threadIdx.x;

    const auto myPointer = begin + tid;

    const location myLocation = myPointer < end ? *myPointer : location{~0, ~0};
    hit_count hits = 1;

    // sum hits of same location run
    for(int i = 1; i < 32; i *= 2) {
        location otherLocation;
        otherLocation.win = __shfl_down_sync(0xFFFFFFFF, myLocation.win, i);
        otherLocation.tgt = __shfl_down_sync(0xFFFFFFFF, myLocation.tgt, i);
        hit_count otherHits = __shfl_down_sync(0xFFFFFFFF, hits, i);
        if(tid < (32-i) && myLocation == otherLocation)
            hits += otherHits;
    }
    // find last thread of same location run
    location otherLocation;
    // tid == 0: otherLocation = myLocation
    otherLocation.win = __shfl_up_sync(0xFFFFFFFF, myLocation.win, 1);
    otherLocation.tgt = __shfl_up_sync(0xFFFFFFFF, myLocation.tgt, 1);

    // find threads which have to write
    bool predicate = !(myLocation == otherLocation);

    if(tid == 0) {
        // predicate is false because myLocation == otherLocation
        // check if first location is the same as previous last location
        if(numLocations > 0 && myLocation == matches[numLocations-1].loc) {
            matches[numLocations-1].hits += hits;
        }
        else {
            // different location -> thread has to write
            predicate = true;
        }
    }

    predicate = predicate && (myLocation.tgt != target_id{~0});

    // find write positions
    const uint32_t mask = __ballot_sync(0xFFFFFFFF, predicate);
    const uint32_t count = __popc(mask);
    const uint32_t numPre = __popc(mask & ((1 << tid) -1));
    const uint32_t locPos = numLocations + numPre;

    // store location and hits sum
    if(predicate) {
        matches[locPos] = {myLocation, hits};

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
    match_candidate& cand, match_candidate * top,
    const ranked_lineage * lineages, taxon_rank mergeBelow)
{
    if(mergeBelow > taxon_rank::Sequence)
        cand.tax = lowest_ranked_ancestor(lineages, cand.tgt, mergeBelow);
    else
        cand.tax = taxon_of_target(lineages, cand.tgt);

    if(!cand.tax) return true;

    int insertPos = MAX_CANDIDATES;
    int taxPos = MAX_CANDIDATES-1;
    // find insert position of cand
    for(int i = 0; i < MAX_CANDIDATES; ++i) {
        if(cand.hits > top[i].hits) {
            insertPos = i;
            break;
        }
    }
    // above sequence level, taxa can occur more than once
    if(mergeBelow != taxon_rank::Sequence) {
        // find same taxon
        for(int i = 0; i < MAX_CANDIDATES; ++i) {
            if(cand.tax == top[i].tax) {
                taxPos = i;
                break;
            }
        }
    }
    // insert except if the same taxon with more hits already exists
    if(taxPos >= insertPos) {
        // move smaller candidates backwards until taxPos
        for(int i = taxPos; i > insertPos; --i) {
            top[i] = top[i-1];
        }
        // insert cand
        top[insertPos] = cand;
    }

    return true;
}


/*************************************************************************//**
 *
 * @brief  produces contiguous window ranges of matches
 *         that are at most 'maxWin' long;
 *
 *****************************************************************************/
 template<int MAX_CANDIDATES>
 __device__
int process_matches(
    match * matches, int numLocations, window_id maxWin,
    match_candidate * top,
    const ranked_lineage * lineages, taxon_rank mergeBelow)
{
    using hit_count = match_candidate::count_type;

    const int tid = threadIdx.x;
    const int numProcessed = min(32, numLocations);

    match_candidate candidate;
    if(tid < numProcessed) {
        match m = matches[tid];
        candidate.tgt = m.loc.tgt;
        candidate.hits = m.hits;
        candidate.pos.beg = m.loc.win;
        candidate.pos.end = m.loc.win;

        for(int i = 1; i < maxWin && (tid+i < numLocations); ++i) {
            match m = matches[tid+i];
            // check if m in same range
            if(candidate.tgt == m.loc.tgt && (candidate.pos.beg+maxWin > m.loc.win)) {
                candidate.hits += m.hits;
                candidate.pos.end = m.loc.win;
            }
            else {
                break;
            }
        }
        // printf("tid: %d tgt: %d beg: %d end: %d hits: %d\n",
            // tid, candidate.tgt, candidate.pos.beg, candidate.pos.end, candidate.hits);

        // top[0] = candidate;
    }


    const target_id myTarget = (tid < numProcessed) ? candidate.tgt : std::numeric_limits<target_id>::max();
    // add thread id to make hit count unique
    hit_count maxHits = (candidate.hits << 5) + 32-1 - tid;

    // segmented max reduction to find candidate with max hits for each target
    for(int i = 1; i < 32; i *= 2) {
        target_id otherTarget = __shfl_up_sync(0xFFFFFFFF, myTarget, i);
        hit_count otherHits = __shfl_up_sync(0xFFFFFFFF, maxHits, i);
        if(myTarget == otherTarget && maxHits < otherHits)
            maxHits = otherHits;

        // propagate max hits
        otherTarget = __shfl_down_sync(0xFFFFFFFF, myTarget, i);
        otherHits = __shfl_down_sync(0xFFFFFFFF, maxHits, i);
        if(myTarget == otherTarget && maxHits < otherHits)
            maxHits = otherHits;
    }
    const bool predicate = (tid < numProcessed) && (32-1 - tid == (maxHits & ((1 << 5) - 1)));

    if(predicate) {
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
    uint32_t numSegments,
    const int * segmentOffsets,
    const Result * locations,
    const window_id * maxWindowsInRange,
    const ranked_lineage * lineages,
    taxon_rank mergeBelow,
    uint32_t maxCandidatesPerQuery,
    match_candidate * topCandidates,
    bool updateTopCandidates)
{
    using hit_count = match_candidate::count_type;

    //TODO reduce size
    __shared__ match matches[96];
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

        int numLocations = 0;
        const window_id maxWin = maxWindowsInRange[bid];
        // if(tid == 0) printf("bid: %d maxWin: %d\n", bid, maxWin);

        while(begin < end || numLocations > 0) {
            for(; (begin < end) && (numLocations <= 32+maxWin-1); begin += 32) {
                numLocations = reduce_locations(begin, end, matches, numLocations);

                // for(int i = tid; i < numLocations; i += 32) {
                //     printf("bid %d tid: %d tgt: %d win: %d hits: %d\n",
                //         bid, tid, matches[tid].loc.tgt, matches[tid].loc.win, matches[tid].hits);
                // }
            }

            // process matches
            const int numProcessed = process_matches<MAX_CANDIDATES>(
                matches, (begin < end) ? numLocations-1 : numLocations, maxWin,
                top, lineages, mergeBelow);

            // if(tid == 0)
                // printf("#locs: %d #processed: %d\n", numLocations, numProcessed);

            // remove processed matches from shared memory
            numLocations -= numProcessed;
            for(int i = tid; i < numLocations; i += 32) {
                matches[i] = matches[numProcessed + i];
            }
        }

        // int tophitsCount = 0;
        // for(int i = 0; i < maxCandidatesPerQuery; ++i) {
        //     if(top[i].hits > 0)
        //         ++tophitsCount;
        // }
        // printf("bid: %d tid: %d: tophits: %d\n", bid, tid, tophitsCount);


        int tophitsUsed = 0;
        int tophitsTotal = 0;

        hit_count hits = top[0].hits;

        while(__ballot_sync(0xFFFFFFFF, hits > 0) && tophitsTotal < maxCandidatesPerQuery) {
            const match_candidate& candidate = top[tophitsUsed];

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
            }
            if(__ballot_sync(0xFFFFFFFF, insert))
                ++tophitsTotal;

            hits = top[tophitsUsed].hits;
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
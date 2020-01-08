#ifndef MC_GPU_RESULT_PROCESSING_H_
#define MC_GPU_RESULT_PROCESSING_H_

#include "config.h"
#include "candidate_generation.h"
#include "dna_encoding.h"

#include "../dep/cudahelpers/cuda_helpers.cuh"

namespace mc {


template<class id_type, class result_type>
__global__
void compact_kernel(
    uint32_t numQueries,
    const int * resultPrefixScan,
    uint32_t maxResultsPerQuery,
    const result_type * results_in,
    result_type * results_out,
    const id_type * segmentIds,
    int * segmentOffsets
)
{
    for(int bid = blockIdx.x; bid < numQueries; bid += gridDim.x) {
        const int begin = (bid > 0) ? resultPrefixScan[bid-1] : 0;
        const int end   = resultPrefixScan[bid];
        const int numResults = end - begin;

        const result_type * inPtr  = results_in + bid*maxResultsPerQuery;
        result_type * outPtr = results_out + begin;

        for(int tid = threadIdx.x; tid < numResults; tid += blockDim.x) {
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
    //     const id_type nextId    = (tid < numQueries) ? segmentIds[tid+1] : -1;

    //     const int end = resultPrefixScan[tid];

    //     // last query of segment sets end offset
    //     if(segmentId != nextId)
    //         segmentOffsets[segmentId+1] = end;
    // }
}


template<class id_type>
__global__
void segment_kernel(
    uint32_t numQueries,
    const int * resultPrefixScan,
    const id_type * segmentIds,
    int * segmentOffsets
)
{
    for(int tid = threadIdx.x + blockIdx.x * blockDim.x; tid < numQueries; tid += blockDim.x * gridDim.x) {
        const id_type segmentId = segmentIds[tid];
        const id_type nextId    = (tid < numQueries) ? segmentIds[tid+1] : -1;

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


template<
    int MAX_CANDIDATES,
    class result_type>
__global__
void generate_top_candidates(
    uint32_t numSegments,
    const int * segmentOffsets,
    const result_type * locations,
    const window_id * maxWindowsInRange,
    const ranked_lineage * lineages,
    taxon_rank mergeBelow,
    uint32_t maxCandidatesPerQuery,
    match_candidate * topCandidates)
{
    using hit_count = match_candidate::count_type;

    match_candidate top[MAX_CANDIDATES];

    for(int tid = threadIdx.x + blockIdx.x * blockDim.x;
            tid < numSegments;
            tid += blockDim.x * gridDim.x)
    {
        auto begin = locations + segmentOffsets[tid];
        auto end = locations + segmentOffsets[tid+1];

        for(int i = 0; i < MAX_CANDIDATES; ++i) {
            top[i].hits = 0;
        }

        for_all_contiguous_window_ranges(
            begin, end, maxWindowsInRange[tid], [&] (match_candidate& cand) {
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
        });

        for(int i = 0; i < maxCandidatesPerQuery; ++i) {
            topCandidates[tid*maxCandidatesPerQuery+i] = top[i];
            // if(top[i].hits > 0)
            //     printf("top: %d tid: %d tgt: %d hits: %d tax: %llu\n", i, tid, top[i].tgt, top[i].hits, top[i].tax);
            // else
            //     break;
        }
    }

}




} // namespace mc


#endif
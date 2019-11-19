#ifndef MC_GPU_ENGINE_H_
#define MC_GPU_ENGINE_H_

#include "config.h"
#include "dna_encoding.h"
#include "sketch_database.h"

#include "../dep/cudahelpers/cuda_helpers.cuh"
#include "../dep/cub/cub/block/block_radix_sort.cuh"
#include "../dep/cub/cub/block/block_store.cuh"
#include "../dep/cub/cub/block/block_scan.cuh"

namespace mc {

// template<bool CanonicalKmers = true>
// __global__
// void extract_kmers(
//     std::uint64_t* seq_in,
//     std::uint32_t* ambig_in,
//     std::uint64_t  size_in,
//     std::uint8_t   k,
//     std::uint64_t* kmers_out,
//     std::uint64_t* size_out)
// {
//     //bases stored from high to low bits
//     //masks get highest bits
//     const std::uint64_t k_mask64 = std::uint64_t(~0) << (64 - 2*k);
//     const std::uint32_t k_mask32 = std::uint32_t(~0) << (16 -   k);

//     for (std::uint64_t tid  = blockIdx.x * blockDim.x + threadIdx.x;
//                        tid  < ((size_in*32)-k+1);
//                        tid += blockDim.x * gridDim.x)
//     {
//         const std::uint64_t seq_slot    = tid / 32;
//         const std::uint8_t  kmer_slot   = tid & 31;
//         const std::uint64_t seq_left    = seq_in[seq_slot] << (2*kmer_slot);
//         const std::uint64_t seq_right   = seq_in[seq_slot+1] >> (64-(2*kmer_slot));
//         const std::uint32_t ambig_left  = ambig_in[seq_slot] << kmer_slot;
//         const std::uint32_t ambig_right = ambig_in[seq_slot+1] >> (32-kmer_slot);

//         //continue only if no bases ambiguous
//         if(!((ambig_left+ambig_right) & k_mask16))
//         {
//             //get highest bits
//             std::uint64_t kmer = (seq_left+seq_right) & k_mask64;
//             //shift kmer to lowest bits
//             kmer >>= (64 - 2*k);

//             if(CanonicalKmers)
//             {
//                 kmer = make_canonical_2bit(kmer);
//             }

//             kmers_out[atomicAggInc(size_out)] = kmer;
//         }
//     }
// }

// template<bool CanonicalKmers = true>
// __global__
// void extract_kmers(
//     std::uint32_t* seq_in,
//     std::uint16_t* ambig_in,
//     std::uint64_t  size_in,
//     std::uint8_t   k,
//     std::uint32_t* kmers_out,
//     std::uint64_t* size_out)
// {
//     //bases stored from high to low bits
//     //masks get highest bits
//     const std::uint32_t k_mask32 = std::uint32_t(~0) << (32 - 2*k);
//     const std::uint16_t k_mask16 = std::uint16_t(~0) << (16 -   k);

//     for (std::uint64_t tid  = blockIdx.x * blockDim.x + threadIdx.x;
//                        tid  < ((size_in*16)-k+1);
//                        tid += blockDim.x * gridDim.x)
//     {
//         const std::uint32_t seq_slot    = tid / 16;
//         const std::uint8_t  kmer_slot   = tid & 15;
//         const std::uint32_t seq_left    = seq_in[seq_slot] << (2*kmer_slot);
//         const std::uint32_t seq_right   = seq_in[seq_slot+1] >> (32-(2*kmer_slot));
//         const std::uint16_t ambig_left  = ambig_in[seq_slot] << kmer_slot;
//         const std::uint16_t ambig_right = ambig_in[seq_slot+1] >> (16-kmer_slot);

//         //continue only if no bases ambiguous
//         if(!((ambig_left+ambig_right) & k_mask16))
//         {
//             //get highest bits
//             std::uint32_t kmer = (seq_left+seq_right) & k_mask32;
//             //shift kmer to lowest bits
//             kmer >>= (32 - 2*k);

//             if(CanonicalKmers)
//             {
//                 kmer = make_canonical_2bit(kmer);
//             }

//             kmers_out[atomicAggInc(size_out)] = kmer;
//         }
//     }
// }



//---------------------------------------------------------------
// BLOCK_THREADS has to be 32
template<
    int BLOCK_THREADS = 32,
    int ITEMS_PER_THREAD = 4>
__global__
void extract_features(
    size_t numTargets,
    const target_id * targetIds,
    const window_id * winOffsets,
    const encodinglen_t * encodeOffsets,
    const encodedseq_t * encodedSeq,
    const encodedambig_t * encodedAmbig,
    numk_t k,
    uint32_t sketchSize,
    size_t windowSize,
    size_t windowStride,
    feature * features_out,
    location * locations_out,
    uint32_t * size_out
)
{
    for(size_t tgt = blockIdx.y; tgt < numTargets; tgt += gridDim.y) {
        const encodinglen_t encodingLength = encodeOffsets[tgt+1] - encodeOffsets[tgt];

        typedef cub::BlockRadixSort<feature, BLOCK_THREADS, ITEMS_PER_THREAD> BlockRadixSortT;

        __shared__ typename BlockRadixSortT::TempStorage temp_storage;

        constexpr uint8_t encBits   = sizeof(encodedseq_t)*CHAR_BIT;
        constexpr uint8_t ambigBits = sizeof(encodedambig_t)*CHAR_BIT;

        //letters stored from high to low bits
        //masks get highest bits
        const encodedseq_t   kmerMask  = encodedseq_t(~0) << (encBits - 2*k);
        const encodedambig_t ambigMask = encodedambig_t(~0) << (ambigBits - k);

        const size_t last = encodingLength * ambigBits - k + 1;
        const window_id numWindows = (last + windowStride - 1) / windowStride;

        const encodedseq_t   * const seq   = encodedSeq   + encodeOffsets[tgt];
        const encodedambig_t * const ambig = encodedAmbig + encodeOffsets[tgt];

        //each block processes one window
        for(window_id win = blockIdx.x; win < numWindows; win += gridDim.x)
        {
            feature items[ITEMS_PER_THREAD];
            uint8_t numItems = 0;
            for(int i=0; i<ITEMS_PER_THREAD; ++i)
            {
                items[i] = feature(~0);
            }

            //each thread extracts one feature
            const size_t first = win * windowStride;
            for(size_t tid  = first + threadIdx.x;
                       tid  < min(first + windowSize - k + 1, last);
                       tid += blockDim.x)
            {
                const std::uint32_t  seq_slot    = tid / (ambigBits);
                const std::uint8_t   kmer_slot   = tid & (ambigBits-1);
                const encodedambig_t ambig_left  = ambig[seq_slot] << kmer_slot;
                const encodedambig_t ambig_right = ambig[seq_slot+1] >> (ambigBits-kmer_slot);
                // cuda-memcheck save version
                // const encodedambig_t ambig_right = kmer_slot ? ambig[seq_slot+1] >> (ambigBits-kmer_slot) : 0;

                //continue only if no bases ambiguous
                if(!((ambig_left+ambig_right) & ambigMask))
                {
                    const encodedseq_t seq_left  = seq[seq_slot] << (2*kmer_slot);
                    const encodedseq_t seq_right = seq[seq_slot+1] >> (encBits-(2*kmer_slot));
                    // cuda-memcheck save version
                    // const encodedseq_t seq_right = kmer_slot ? seq[seq_slot+1] >> (encBits-(2*kmer_slot)) : 0;

                    //get highest bits
                    kmer_type kmer = (seq_left+seq_right) & kmerMask;
                    //shift kmer to lowest bits
                    kmer >>= (encBits - 2*k);

                    kmer = make_canonical_2bit(kmer);

                    items[numItems++] = sketching_hash{}(kmer);
                }
                // else
                // {
                //     printf("ambiguous\n");
                // }
            }

            BlockRadixSortT(temp_storage).SortBlockedToStriped(items);

            //output <sketchSize> unique features
            uint32_t kmerCounter = 0;
            for(int i=0; i*BLOCK_THREADS<windowStride && kmerCounter<sketchSize; ++i)
            {
                //load candidate and successor
                const feature kmer = items[i];
                feature nextKmer = threadIdx.x ? items[i] : items[(i+1) % ITEMS_PER_THREAD];
                nextKmer = __shfl_sync(0xFFFFFFFF, nextKmer, threadIdx.x+1);

                //find valid candidates
                const bool predicate = (kmer != feature(~0)) && (kmer != nextKmer);

                const uint32_t mask = __ballot_sync(0xFFFFFFFF, predicate);
                const uint32_t count = __popc(mask);
                const uint32_t capped = min(sketchSize-kmerCounter, count);

                //get output position in global memory
                uint64_t pos;
                if (threadIdx.x == 0) {
                    // printf("kmerCounter: %d count: %d capped: %d\n", kmerCounter, count, capped);
                    pos = atomicAdd(size_out, capped);
                }
                const uint32_t numPre = __popc(mask & ((1 << threadIdx.x) -1));
                pos = __shfl_sync(0xFFFFFFFF, pos, 0) + numPre;

                //write to global memory
                if(predicate && (numPre < capped)) {
                    const target_id targetId = targetIds[tgt];
                    const window_id winOffset = winOffsets[tgt];

                    features_out[pos]  = kmer;
                    locations_out[pos] = location{winOffset+win, targetId};
                }

                kmerCounter += capped;
            }
        }
    }
}



//---------------------------------------------------------------
// BLOCK_THREADS has to be 32
template<
    int BLOCK_THREADS,
    int ITEMS_PER_THREAD,
    class Hashtable,
    class bucket_size_type,
    class result_type>
__global__
void gpu_hahstable_query(
    Hashtable hashtable,
    size_t probingLength,
    uint32_t numQueries,
    const encodinglen_t  * encodeOffsets,
    const encodedseq_t   * encodedSeq,
    const encodedambig_t * encodedAmbig,
    numk_t k,
    uint32_t sketchSize,
    uint32_t windowSize,
    uint32_t windowStride,
    const location * locations,
    bucket_size_type maxLocationsPerFeature,
    result_type * results_out,
    int         * resultCounts
)
{
    //each block processes one query
    for(size_t queryId = blockIdx.x; queryId < numQueries; queryId += gridDim.x) {
        const encodinglen_t encodingLength = encodeOffsets[queryId+1] - encodeOffsets[queryId];

        typedef cub::BlockRadixSort<feature, BLOCK_THREADS, ITEMS_PER_THREAD> BlockRadixSortT;

        __shared__ union {
            typename BlockRadixSortT::TempStorage sort;
            uint64_t sketch[16];
        } tempStorage;

        constexpr uint8_t encBits   = sizeof(encodedseq_t)*CHAR_BIT;
        constexpr uint8_t ambigBits = sizeof(encodedambig_t)*CHAR_BIT;

        //letters stored from high to low bits
        //masks get highest bits
        const encodedseq_t   kmerMask  = encodedseq_t(~0) << (encBits - 2*k);
        const encodedambig_t ambigMask = encodedambig_t(~0) << (ambigBits - k);

        const size_t last = encodingLength * ambigBits - k + 1;

        const encodedseq_t   * const seq   = encodedSeq   + encodeOffsets[queryId];
        const encodedambig_t * const ambig = encodedAmbig + encodeOffsets[queryId];

        //each block processes one window
        {
            feature items[ITEMS_PER_THREAD];
            uint8_t numItems = 0;
            for(int i=0; i<ITEMS_PER_THREAD; ++i)
            {
                items[i] = feature(~0);
            }

            //each thread extracts one feature
            for(size_t tid  = threadIdx.x;
                       tid  < last;
                       tid += BLOCK_THREADS)
            {
                const std::uint32_t  seq_slot    = tid / (ambigBits);
                const std::uint8_t   kmer_slot   = tid & (ambigBits-1);
                const encodedambig_t ambig_left  = ambig[seq_slot] << kmer_slot;
                const encodedambig_t ambig_right = ambig[seq_slot+1] >> (ambigBits-kmer_slot);
                // cuda-memcheck save version
                // const encodedambig_t ambig_right = kmer_slot ? ambig[seq_slot+1] >> (ambigBits-kmer_slot) : 0;

                //continue only if no bases ambiguous
                if(!((ambig_left+ambig_right) & ambigMask))
                {
                    const encodedseq_t seq_left  = seq[seq_slot] << (2*kmer_slot);
                    const encodedseq_t seq_right = seq[seq_slot+1] >> (encBits-(2*kmer_slot));
                    // cuda-memcheck save version
                    // const encodedseq_t seq_right = kmer_slot ? seq[seq_slot+1] >> (encBits-(2*kmer_slot)) : 0;

                    //get highest bits
                    kmer_type kmer = (seq_left+seq_right) & kmerMask;
                    //shift kmer to lowest bits
                    kmer >>= (encBits - 2*k);

                    kmer = make_canonical_2bit(kmer);

                    items[numItems++] = sketching_hash{}(kmer);
                }
                // else
                // {
                //     printf("ambiguous\n");
                // }
            }

            BlockRadixSortT(tempStorage.sort).SortBlockedToStriped(items);

            //output <sketchSize> unique features
            uint32_t kmerCounter = 0;
            for(int i=0; i<ITEMS_PER_THREAD && kmerCounter<sketchSize; ++i)
            {
                //load candidate and successor
                const feature kmer = items[i];
                feature nextKmer = threadIdx.x ? items[i] : items[(i+1) % ITEMS_PER_THREAD];
                nextKmer = __shfl_sync(0xFFFFFFFF, nextKmer, threadIdx.x+1);

                //find valid candidates
                const bool predicate = (kmer != feature(~0)) && (kmer != nextKmer);

                const uint32_t mask = __ballot_sync(0xFFFFFFFF, predicate);
                const uint32_t count = __popc(mask);

                //get position
                const uint32_t numPre = __popc(mask & ((1 << threadIdx.x) -1));
                uint32_t sketchPos = kmerCounter + numPre;

                //write features to shared memory
                if(predicate && (sketchPos < sketchSize)) {
                    tempStorage.sketch[sketchPos] = kmer;
                    // printf("feature: %d sketchPos: %d\n", kmer, sketchPos);
                }

                kmerCounter += count;
            }

            //cap kmer count
            kmerCounter = min(kmerCounter, sketchSize);

            //query hashtable for bucket locations & sizes
            namespace cg = cooperative_groups;

            const auto group =
                cg::tiled_partition<Hashtable::cg_size>(cg::this_thread_block());

            for(uint32_t i = threadIdx.x; i < kmerCounter*Hashtable::cg_size; i += BLOCK_THREADS) {
                const uint8_t gid = i / Hashtable::cg_size;
                typename Hashtable::value_type valuesOffset = 0;
                //if key not found valuesOffset stays 0
                const auto status = hashtable.retrieve(tempStorage.sketch[gid], valuesOffset, group, probingLength);


                if(group.thread_rank() == 0) {
                    tempStorage.sketch[gid] = valuesOffset;

                    //if key not found numValues = 0
                    // bucket_size_type numValues = bucket_size_type(valuesOffset);
                    // valuesOffset >>= sizeof(bucket_size_type)*CHAR_BIT;

                    // printf("status %d\n", status.has_any());
                    // printf("status %d\n", status.has_key_not_found());
                    // printf("status %d offset: %llu numValues: %d\n", status.has_any(), valuesOffset, uint32_t(numValues));
                }
            }
            __syncthreads();

            int totalNumValues = 0;

            result_type * outPtr = results_out + queryId*sketchSize*maxLocationsPerFeature;

            //copy locations of found features into results_out
            for(uint32_t s = 0; s < kmerCounter; ++s) {
                typename Hashtable::value_type valuesOffset = tempStorage.sketch[s];

                bucket_size_type numValues = bucket_size_type(valuesOffset);
                valuesOffset >>= sizeof(bucket_size_type)*CHAR_BIT;

                for(uint32_t i = threadIdx.x; i < numValues; i += BLOCK_THREADS) {
                    //copy found locations
                    outPtr[totalNumValues + i] = locations[valuesOffset + i];
                }

                totalNumValues += numValues;
            }

            if(threadIdx.x == 0) {
                //store number of results of this query
                resultCounts[queryId] = totalNumValues;
            }
        }
    }
}


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


} // namespace mc


#endif
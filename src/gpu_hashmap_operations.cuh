#ifndef MC_GPU_HASHMAP_OPERATIONS_H_
#define MC_GPU_HASHMAP_OPERATIONS_H_

#include "config.h"
#include "dna_encoding.h"

#include "../dep/cudahelpers/cuda_helpers.cuh"
#include "cub/block/block_radix_sort.cuh"

namespace mc {


/****************************************************************
 * @brief output (at most) <sketchSize> unique features
 */
template<
    int ITEMS_PER_THREAD,
    class Feature,
    class Feature2>
__device__ __inline__
uint32_t unique_sketch(
    const Feature items[ITEMS_PER_THREAD],
    Feature2 * sketch,
    uint32_t sketchSize)
{
        uint32_t uniqueFeatureCounter = 0;

        for(int i=0; i<ITEMS_PER_THREAD && uniqueFeatureCounter<sketchSize; ++i)
        {
            // load candidate and successor
            const Feature kmer = items[i];
            Feature nextKmer = threadIdx.x ? items[i] : items[(i+1) % ITEMS_PER_THREAD];
            nextKmer = __shfl_sync(0xFFFFFFFF, nextKmer, threadIdx.x+1);

            // find valid candidates
            const bool predicate = (kmer != Feature(~0)) && (kmer != nextKmer);

            const uint32_t mask = __ballot_sync(0xFFFFFFFF, predicate);
            const uint32_t count = __popc(mask);

            const uint32_t numPre = __popc(mask & ((1 << threadIdx.x) -1));
            uint32_t sketchPos = uniqueFeatureCounter + numPre;

            // write features to shared memory
            if(predicate && (sketchPos < sketchSize)) {
                sketch[sketchPos] = kmer;
                // printf("sketchPos: %2d feature: %d\n", sketchPos, kmer);
            }

            uniqueFeatureCounter += count;
        }
        // cap counter
        uniqueFeatureCounter = min(uniqueFeatureCounter, sketchSize);

        return uniqueFeatureCounter;
}


/****************************************************************
 * @brief for each feature in sketch insert location into hash table;
 *        a location consists of of target id and window id
 */
template<
    int BLOCK_THREADS,
    class Hashtable,
    class Feature>
__device__ __inline__
void insert_into_hashtable(
    Hashtable& hashtable,
    const Feature * sketch,
    uint32_t sketchCounter,
    target_id targetId,
    window_id winOffset)
{
    using location = typename Hashtable::value_type;

    namespace cg = cooperative_groups;

    const auto group =
        cg::tiled_partition<Hashtable::cg_size()>(cg::this_thread_block());

    for(uint32_t i = threadIdx.x; i < sketchCounter*Hashtable::cg_size(); i += BLOCK_THREADS) {
        const uint8_t gid = i / Hashtable::cg_size();

        const auto status = hashtable.insert(
            sketch[gid], location{winOffset, targetId}, group);

        if(group.thread_rank() == 0 && status.has_any()) {
            printf("status %d\n", status.has_any());
        }
    }
}


/****************************************************************
 * @brief query the hash table with each feature in sketch;
 *        write the result back into sketch in place of the feature
 */
template<
    int BLOCK_THREADS,
    class Hashtable,
    class Feature>
__device__ __inline__
void query_hashtable(
    Hashtable& hashtable,
    Feature * sketch,
    uint32_t sketchCounter)
{
    namespace cg = cooperative_groups;

    const auto group =
        cg::tiled_partition<Hashtable::cg_size()>(cg::this_thread_block());

    for(uint32_t i = threadIdx.x; i < sketchCounter*Hashtable::cg_size(); i += BLOCK_THREADS) {
        const uint8_t gid = i / Hashtable::cg_size();
        typename Hashtable::value_type valuesOffset = 0;
        //if key not found valuesOffset stays 0
        const auto status = hashtable.retrieve(
            sketch[gid], valuesOffset, group);

        if(group.thread_rank() == 0) {
            if(status.has_any()) {
                printf("status %d\n", status.has_any());
                // printf("status %d\n", status.has_key_not_found());
            }

            sketch[gid] = valuesOffset;
        }
    }
}


/****************************************************************
 * @brief for each feature in sketch copy its loctions to the output
 */
template<
    int BLOCK_THREADS,
    class Feature,
    class Location,
    class BucketSizeT>
__device__ __inline__
uint32_t copy_loctions(
    Feature * sketch,
    uint32_t sketchSize,
    const Location * locations,
    BucketSizeT maxLocationsPerFeature,
    Location * out)
{
    using bucket_size_type = BucketSizeT;

    uint32_t totalLocations = 0;

    for(uint32_t s = 0; s < sketchSize; ++s) {
        auto bucketOffset = sketch[s];
        bucket_size_type bucketSize = min(bucket_size_type(bucketOffset), maxLocationsPerFeature);
        bucketOffset >>= sizeof(bucket_size_type)*CHAR_BIT;

        //copy locations
        for(uint32_t i = threadIdx.x; i < bucketSize; i += BLOCK_THREADS) {
            out[totalLocations + i] = locations[bucketOffset + i];
        }

        totalLocations += bucketSize;
    }

    return totalLocations;
}


//---------------------------------------------------------------
// BLOCK_THREADS has to be 32
template<
    int BLOCK_THREADS = 32,
    int ITEMS_PER_THREAD = 4,
    class Hashtable>
__global__
__launch_bounds__(BLOCK_THREADS)
void insert_features(
    Hashtable hashtable,
    uint32_t numTargets,
    const target_id * targetIds,
    const window_id * winOffsets,
    const encodinglen_t * encodeOffsets,
    const encodedseq_t * encodedSeq,
    const encodedambig_t * encodedAmbig,
    numk_t k,
    uint32_t maxSketchSize,
    size_t windowSize,
    size_t windowStride
)
{
    using BlockRadixSortT = cub::BlockRadixSort<feature, BLOCK_THREADS, ITEMS_PER_THREAD>;

    __shared__ union {
        typename BlockRadixSortT::TempStorage sort;
        //TODO MAX_SKETCH_SIZE
        feature sketch[16];
    } tempStorage;

    for(size_t tgt = blockIdx.y; tgt < numTargets; tgt += gridDim.y) {
        const encodinglen_t encodingLength = encodeOffsets[tgt+1] - encodeOffsets[tgt];

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
            uint8_t numItems = 0;
            feature items[ITEMS_PER_THREAD];
            //initialize thread local feature store
            for(int i=0; i<ITEMS_PER_THREAD; ++i)
                items[i] = feature(~0);

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
                const encodedseq_t seq_left  = seq[seq_slot] << (2*kmer_slot);
                const encodedseq_t seq_right = seq[seq_slot+1] >> (encBits-(2*kmer_slot));
                // cuda-memcheck save version
                // const encodedseq_t seq_right = kmer_slot ? seq[seq_slot+1] >> (encBits-(2*kmer_slot)) : 0;

                //get highest bits
                kmer_type kmer = (seq_left+seq_right) & kmerMask;
                //shift kmer to lowest bits
                kmer >>= (encBits - 2*k);

                kmer = make_canonical_2bit(kmer);

                bool pred = !((ambig_left+ambig_right) & ambigMask);
                items[numItems++] = pred ? sketching_hash{}(kmer) : feature(~0);

                // if(!((ambig_left+ambig_right) & ambigMask))
                // {
                //     items[numItems++] = sketching_hash{}(kmer);
                // }
                // else
                // {
                //     printf("ambiguous\n");
                // }
            }

            BlockRadixSortT(tempStorage.sort).SortBlockedToStriped(items);

            uint32_t sketchCounter =
                unique_sketch<ITEMS_PER_THREAD>(items, tempStorage.sketch, maxSketchSize);

            __syncthreads();

            const target_id targetId = targetIds[tgt];
            const window_id winOffset = winOffsets[tgt]+win ;

            insert_into_hashtable<BLOCK_THREADS>(
                hashtable, tempStorage.sketch, sketchCounter, targetId, winOffset);
        }
    }
}



//---------------------------------------------------------------
// BLOCK_THREADS has to be 32
template<
    int BLOCK_THREADS,
    int ITEMS_PER_THREAD,
    class Hashtable,
    class BucketSizeT,
    class Location>
__global__
void gpu_hahstable_query(
    Hashtable hashtable,
    uint32_t numQueries,
    const encodinglen_t * sequenceOffsets,
    const char * sequences,
    numk_t k,
    uint32_t maxSketchSize,
    uint32_t windowSize,
    uint32_t windowStride,
    const Location * locations,
    BucketSizeT maxLocationsPerFeature,
    Location * locations_out,
    int      * locationCounts
)
{
    using location = Location;
    using bucket_size_type = BucketSizeT;
    using encodedseq_t = uint32_t;
    using encodedambig_t = uint16_t;
    using BlockRadixSortT = cub::BlockRadixSort<feature, BLOCK_THREADS, ITEMS_PER_THREAD>;

    constexpr uint8_t encBits   = sizeof(encodedseq_t)*CHAR_BIT;
    constexpr uint8_t ambigBits = sizeof(encodedambig_t)*CHAR_BIT;
    constexpr uint8_t encodedBlocksPerWindow = 128 / ambigBits;

    __shared__ union {
        typename BlockRadixSortT::TempStorage sort;
        //TODO MAX_SKETCH_SIZE
        uint64_t sketch[16];
        struct {
            encodedseq_t seq[encodedBlocksPerWindow];
            encodedambig_t ambig[encodedBlocksPerWindow];
        } encodedWindow;
    } tempStorage;

    //each block processes one query (= one window)
    for(int queryId = blockIdx.x; queryId < numQueries; queryId += gridDim.x) {
        const encodinglen_t sequenceLength = sequenceOffsets[queryId+1] - sequenceOffsets[queryId];

        //only process non-empty queries
        if(sequenceLength) {
            if(threadIdx.x < encodedBlocksPerWindow)
                tempStorage.encodedWindow.ambig[threadIdx.x] = encodedambig_t(~0);
            __syncthreads();

            const char * sequenceBegin = sequences + sequenceOffsets[queryId];
            const encodinglen_t sequenceLengthPadded = (sequenceLength+BLOCK_THREADS-1) / BLOCK_THREADS * BLOCK_THREADS;
            const int lane = ambigBits - threadIdx.x % ambigBits - 1;
            constexpr int numSubwarps = BLOCK_THREADS / ambigBits;
            const int subwarpId = threadIdx.x / ambigBits;
            int blockCounter = subwarpId;
            for(int tid = threadIdx.x; tid < sequenceLengthPadded; tid += BLOCK_THREADS) {
                encodedseq_t seq = 0;
                encodedambig_t ambig = 0;
                const char c = (tid < sequenceLength) ? sequenceBegin[tid] : 'N';
                switch(c) {
                    case 'A': case 'a': break;
                    case 'C': case 'c': seq |= 1; break;
                    case 'G': case 'g': seq |= 2; break;
                    case 'T': case 't': seq |= 3; break;
                    default: ambig |= 1; break;
                }
                seq <<= 2*lane;
                ambig <<= lane;

                for(int stride = 1; stride < ambigBits; stride <<= 1) {
                    seq   |= __shfl_xor_sync(0xFFFFFFFF, seq, stride);
                    ambig |= __shfl_xor_sync(0xFFFFFFFF, ambig, stride);
                }

                if(lane == 0) {
                    tempStorage.encodedWindow.seq[blockCounter]   = seq;
                    tempStorage.encodedWindow.ambig[blockCounter] = ambig;
                    // printf("query: %u tid: %d seq: %u ambig: %u\n", queryId, tid, seq, ambig);
                }
                blockCounter += numSubwarps;
            }
            __syncthreads();

            //letters stored from high to low bits
            //masks get highest bits
            const encodedseq_t   kmerMask  = encodedseq_t(~0) << (encBits - 2*k);
            const encodedambig_t ambigMask = encodedambig_t(~0) << (ambigBits - k);

            const size_t last = sequenceLength - k + 1;

            uint8_t numItems = 0;
            feature items[ITEMS_PER_THREAD];
            //initialize thread local feature store
            for(int i=0; i<ITEMS_PER_THREAD; ++i)
                items[i] = feature(~0);

            //each thread extracts one feature
            for(int tid  = threadIdx.x;
                    tid  < last;
                    tid += BLOCK_THREADS)
            {
                const std::uint32_t  seq_slot    = tid / (ambigBits);
                const std::uint8_t   kmer_slot   = tid & (ambigBits-1);
                const encodedambig_t ambig_left  = tempStorage.encodedWindow.ambig[seq_slot] << kmer_slot;
                const encodedambig_t ambig_right = tempStorage.encodedWindow.ambig[seq_slot+1] >> (ambigBits-kmer_slot);
                // cuda-memcheck save version
                // const encodedambig_t ambig_right = kmer_slot ? tempStorage.encodedWindow.ambig[seq_slot+1] >> (ambigBits-kmer_slot) : 0;

                //continue only if no bases ambiguous
                const encodedseq_t seq_left  = tempStorage.encodedWindow.seq[seq_slot] << (2*kmer_slot);
                const encodedseq_t seq_right = tempStorage.encodedWindow.seq[seq_slot+1] >> (encBits-(2*kmer_slot));
                // cuda-memcheck save version
                // const encodedseq_t seq_right = kmer_slot ? tempStorage.encodedWindow.seq[seq_slot+1] >> (encBits-(2*kmer_slot)) : 0;

                //get highest bits
                kmer_type kmer = (seq_left+seq_right) & kmerMask;
                //shift kmer to lowest bits
                kmer >>= (encBits - 2*k);

                kmer = make_canonical_2bit(kmer);

                bool pred = !((ambig_left+ambig_right) & ambigMask);
                items[numItems++] = pred ? sketching_hash{}(kmer) : feature(~0);

                // if(!((ambig_left+ambig_right) & ambigMask))
                // {
                //     items[numItems++] = sketching_hash{}(kmer);
                // }
                // else
                // {
                //     printf("ambiguous\n");
                // }
            }
            __syncthreads();

            BlockRadixSortT(tempStorage.sort).SortBlockedToStriped(items);

            uint32_t realSketchSize =
                unique_sketch<ITEMS_PER_THREAD>(items, tempStorage.sketch, maxSketchSize);

            __syncthreads();

            query_hashtable<BLOCK_THREADS>(
                hashtable, tempStorage.sketch, realSketchSize);

            __syncthreads();

            location * out = locations_out + queryId*maxSketchSize*maxLocationsPerFeature;

            uint32_t numLocations = copy_loctions<BLOCK_THREADS>(
                tempStorage.sketch, realSketchSize, locations, maxLocationsPerFeature, out);

            if(threadIdx.x == 0) {
                //store number of results of this query
                locationCounts[queryId] = numLocations;
            }

            __syncthreads();
        }
    }
}


} // namespace mc

#endif

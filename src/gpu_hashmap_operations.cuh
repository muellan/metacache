#ifndef MC_GPU_HASHMAP_OPERATIONS_H_
#define MC_GPU_HASHMAP_OPERATIONS_H_

#include "config.h"
#include "dna_encoding.h"

#include "../dep/hpc_helpers/include/cuda_helpers.cuh"
#include "../dep/bb_segsort/src/bb_exch_keys.cuh"

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
    const int warpLane = threadIdx.x % WARPSIZE;

    uint32_t uniqueFeatureCounter = 0;

    #pragma unroll
    for(int i=0; i<ITEMS_PER_THREAD; ++i)
    {
        // load candidate and successor
        const Feature kmer = items[i];
        Feature nextKmer = warpLane ? items[i] : items[(i+1) % ITEMS_PER_THREAD];
        nextKmer = __shfl_sync(0xFFFFFFFF, nextKmer, warpLane+1);

        // find valid candidates
        const bool predicate = (kmer != Feature(~0)) && (kmer != nextKmer);

        const uint32_t mask = __ballot_sync(0xFFFFFFFF, predicate);
        const uint32_t count = __popc(mask);

        const uint32_t numPre = __popc(mask & ((1 << warpLane) -1));
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

    const int warpLane = threadIdx.x % WARPSIZE;
    const int groupId = warpLane / Hashtable::cg_size();
    constexpr int groupsPerWarp = WARPSIZE / Hashtable::cg_size();

    for(int i = groupId; i < sketchCounter; i += groupsPerWarp)
    {
        const auto status = hashtable.insert(
            sketch[i], location{winOffset, targetId}, group);

        // if(group.thread_rank() == 0 && status.has_any()) {
        //     printf("status %d\n", status.has_any());
        // }
    }
}


/****************************************************************
 * @brief query the hash table with each feature in sketch;
 *        write the result back into sketch in place of the feature
 */
template<
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

    const int warpLane = threadIdx.x % WARPSIZE;
    const int groupId = warpLane / Hashtable::cg_size();
    constexpr int groupsPerWarp = WARPSIZE / Hashtable::cg_size();

    for(int i = groupId; i < sketchCounter; i += groupsPerWarp)
    {
        typename Hashtable::value_type valuesOffset = 0;
        //if key not found valuesOffset stays 0
        const auto status = hashtable.retrieve(
            sketch[i], valuesOffset, group);

        if(group.thread_rank() == 0) {
            // if(status.has_any()) {
                // printf("status %d\n", status.has_any());
                // printf("status %d\n", status.has_key_not_found());
            // }

            sketch[i] = valuesOffset;
        }
    }
}


/****************************************************************
 * @brief for each feature in sketch copy its loctions to the output
 */
template<
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
    using location = Location;
    using bucket_size_type = BucketSizeT;

    const int warpLane = threadIdx.x % WARPSIZE;

    uint32_t totalLocations = 0;

    for(uint32_t s = 0; s < sketchSize; ++s) {
        auto bucketOffset = sketch[s];
        bucket_size_type bucketSize = min(bucket_size_type(bucketOffset), maxLocationsPerFeature);
        bucketOffset >>= sizeof(bucket_size_type)*CHAR_BIT;

        //copy locations
        for(uint32_t i = warpLane; i < bucketSize; i += WARPSIZE) {
            out[totalLocations + i] = locations[bucketOffset + i];
        }

        totalLocations += bucketSize;
    }

    // add padding for aligned processing
    const uint32_t padding = totalLocations % 2;
    if(padding && (warpLane == 0)) {
        out[totalLocations] = location{~0, ~0};
    }

    return totalLocations + padding;
}



/****************************************************************
 * @brief query the hash table with each feature in sketch;
 *        write the result back into sketch in place of the feature
 */
template<
    int MAX_SKETCH_SIZE,
    class Hashtable>
__device__ __inline__
uint32_t query_bucket_hashtable(
    Hashtable& hashtable,
    feature * sketch,
    uint32_t * sizes,
    uint32_t sketchCounter,
    location * out)
{
    namespace cg = cooperative_groups;

    const auto group =
        cg::tiled_partition<Hashtable::cg_size()>(cg::this_thread_block());

    const int warpLane = threadIdx.x % WARPSIZE;
    const int groupId = warpLane / Hashtable::cg_size();
    constexpr int groupsPerWarp = WARPSIZE / Hashtable::cg_size();

    for(int i = groupId; i < sketchCounter; i += groupsPerWarp)
    {
        uint64_t numValues = 0;
        //if key not found numValues stays 0
        const auto status = hashtable.num_values(
            sketch[i], numValues, group);

        if(group.thread_rank() == 0) {
            // if(status.has_any()) {
                // printf("status %d\n", status.has_any());
                // printf("status %d\n", status.has_key_not_found());
            // }

            sizes[i+1] = numValues;
        }
    }
    if(warpLane == 0) sizes[0] = 0;

    __syncwarp();

    uint32_t numValues = warpLane <= sketchCounter ? sizes[warpLane] : 0;

    for(int i = 1; i <= MAX_SKETCH_SIZE; i *= 2) {
        const uint32_t other = __shfl_up_sync(0xFFFFFFFF, numValues, i);
        if(warpLane > i)
            numValues += other;
    }
    sizes[warpLane] = numValues;

    for(int i = groupId; i < sketchCounter; i += groupsPerWarp)
    {
        const uint32_t offset = sizes[i];
        const uint32_t size  = sizes[i+1] - offset;

        // only retrieve non empty buckets
        if(size) {
            uint64_t dummy;
            const auto status = hashtable.retrieve(
                sketch[i], out + offset, dummy, group);
        }
    }

    // add padding for aligned processing
    const uint32_t totalLocations = sizes[MAX_SKETCH_SIZE];
    const uint32_t padding = totalLocations % 2;
    if(padding && (warpLane == 0)) {
        out[totalLocations] = location{~0, ~0};
    }

    return totalLocations + padding;
}


//---------------------------------------------------------------
template<class K>
__device__ __inline__
void warp_sort_128(K rg_k[4])
{
    const int warpLane = threadIdx.x % WARPSIZE;
    const int bit1 = (warpLane>>0)&0x1;
    const int bit2 = (warpLane>>1)&0x1;
    const int bit3 = (warpLane>>2)&0x1;
    const int bit4 = (warpLane>>3)&0x1;
    const int bit5 = (warpLane>>4)&0x1;

    // sort 128 elements
    // exch_intxn: switch to exch_local()
    CMP_SWP_KEY(K,rg_k[0] ,rg_k[1] );
    CMP_SWP_KEY(K,rg_k[2] ,rg_k[3] );
    // exch_intxn: switch to exch_local()
    CMP_SWP_KEY(K,rg_k[0] ,rg_k[3] );
    CMP_SWP_KEY(K,rg_k[1] ,rg_k[2] );
    // exch_paral: switch to exch_local()
    CMP_SWP_KEY(K,rg_k[0] ,rg_k[1] );
    CMP_SWP_KEY(K,rg_k[2] ,rg_k[3] );
    // exch_intxn: generate exch_intxn_keys()
    exch_intxn_keys(rg_k[0] ,rg_k[1] ,rg_k[2] ,rg_k[3] ,
                    0x1,bit1);
    // exch_paral: switch to exch_local()
    CMP_SWP_KEY(K,rg_k[0] ,rg_k[2] );
    CMP_SWP_KEY(K,rg_k[1] ,rg_k[3] );
    // exch_paral: switch to exch_local()
    CMP_SWP_KEY(K,rg_k[0] ,rg_k[1] );
    CMP_SWP_KEY(K,rg_k[2] ,rg_k[3] );
    // exch_intxn: generate exch_intxn_keys()
    exch_intxn_keys(rg_k[0] ,rg_k[1] ,rg_k[2] ,rg_k[3] ,
                    0x3,bit2);
    // exch_paral: generate exch_paral_keys()
    exch_paral_keys(rg_k[0] ,rg_k[1] ,rg_k[2] ,rg_k[3] ,
                    0x1,bit1);
    // exch_paral: switch to exch_local()
    CMP_SWP_KEY(K,rg_k[0] ,rg_k[2] );
    CMP_SWP_KEY(K,rg_k[1] ,rg_k[3] );
    // exch_paral: switch to exch_local()
    CMP_SWP_KEY(K,rg_k[0] ,rg_k[1] );
    CMP_SWP_KEY(K,rg_k[2] ,rg_k[3] );
    // exch_intxn: generate exch_intxn_keys()
    exch_intxn_keys(rg_k[0] ,rg_k[1] ,rg_k[2] ,rg_k[3] ,
                    0x7,bit3);
    // exch_paral: generate exch_paral_keys()
    exch_paral_keys(rg_k[0] ,rg_k[1] ,rg_k[2] ,rg_k[3] ,
                    0x2,bit2);
    // exch_paral: generate exch_paral_keys()
    exch_paral_keys(rg_k[0] ,rg_k[1] ,rg_k[2] ,rg_k[3] ,
                    0x1,bit1);
    // exch_paral: switch to exch_local()
    CMP_SWP_KEY(K,rg_k[0] ,rg_k[2] );
    CMP_SWP_KEY(K,rg_k[1] ,rg_k[3] );
    // exch_paral: switch to exch_local()
    CMP_SWP_KEY(K,rg_k[0] ,rg_k[1] );
    CMP_SWP_KEY(K,rg_k[2] ,rg_k[3] );
    // exch_intxn: generate exch_intxn_keys()
    exch_intxn_keys(rg_k[0] ,rg_k[1] ,rg_k[2] ,rg_k[3] ,
                    0xf,bit4);
    // exch_paral: generate exch_paral_keys()
    exch_paral_keys(rg_k[0] ,rg_k[1] ,rg_k[2] ,rg_k[3] ,
                    0x4,bit3);
    // exch_paral: generate exch_paral_keys()
    exch_paral_keys(rg_k[0] ,rg_k[1] ,rg_k[2] ,rg_k[3] ,
                    0x2,bit2);
    // exch_paral: generate exch_paral_keys()
    exch_paral_keys(rg_k[0] ,rg_k[1] ,rg_k[2] ,rg_k[3] ,
                    0x1,bit1);
    // exch_paral: switch to exch_local()
    CMP_SWP_KEY(K,rg_k[0] ,rg_k[2] );
    CMP_SWP_KEY(K,rg_k[1] ,rg_k[3] );
    // exch_paral: switch to exch_local()
    CMP_SWP_KEY(K,rg_k[0] ,rg_k[1] );
    CMP_SWP_KEY(K,rg_k[2] ,rg_k[3] );
    // exch_intxn: generate exch_intxn_keys()
    exch_intxn_keys(rg_k[0] ,rg_k[1] ,rg_k[2] ,rg_k[3] ,
                    0x1f,bit5);
    // exch_paral: generate exch_paral_keys()
    exch_paral_keys(rg_k[0] ,rg_k[1] ,rg_k[2] ,rg_k[3] ,
                    0x8,bit4);
    // exch_paral: generate exch_paral_keys()
    exch_paral_keys(rg_k[0] ,rg_k[1] ,rg_k[2] ,rg_k[3] ,
                    0x4,bit3);
    // exch_paral: generate exch_paral_keys()
    exch_paral_keys(rg_k[0] ,rg_k[1] ,rg_k[2] ,rg_k[3] ,
                    0x2,bit2);
    // exch_paral: generate exch_paral_keys()
    exch_paral_keys(rg_k[0] ,rg_k[1] ,rg_k[2] ,rg_k[3] ,
                    0x1,bit1);
    // exch_paral: switch to exch_local()
    CMP_SWP_KEY(K,rg_k[0] ,rg_k[2] );
    CMP_SWP_KEY(K,rg_k[1] ,rg_k[3] );
    // exch_paral: switch to exch_local()
    CMP_SWP_KEY(K,rg_k[0] ,rg_k[1] );
    CMP_SWP_KEY(K,rg_k[2] ,rg_k[3] );

    // transpose blocked to striped
    if(warpLane&0x1 ) SWP_KEY(K, rg_k[0] , rg_k[1] );
    if(warpLane&0x1 ) SWP_KEY(K, rg_k[2] , rg_k[3] );
    rg_k[1]  = __shfl_xor_sync(0xffffffff,rg_k[1] , 0x1 );
    rg_k[3]  = __shfl_xor_sync(0xffffffff,rg_k[3] , 0x1 );
    if(warpLane&0x1 ) SWP_KEY(K, rg_k[0] , rg_k[1] );
    if(warpLane&0x1 ) SWP_KEY(K, rg_k[2] , rg_k[3] );
    if(warpLane&0x2 ) SWP_KEY(K, rg_k[0] , rg_k[2] );
    if(warpLane&0x2 ) SWP_KEY(K, rg_k[1] , rg_k[3] );
    rg_k[2]  = __shfl_xor_sync(0xffffffff,rg_k[2] , 0x2 );
    rg_k[3]  = __shfl_xor_sync(0xffffffff,rg_k[3] , 0x2 );
    if(warpLane&0x2 ) SWP_KEY(K, rg_k[0] , rg_k[2] );
    if(warpLane&0x2 ) SWP_KEY(K, rg_k[1] , rg_k[3] );
    if(warpLane&0x4 ) SWP_KEY(K, rg_k[0] , rg_k[1] );
    if(warpLane&0x4 ) SWP_KEY(K, rg_k[2] , rg_k[3] );
    rg_k[1]  = __shfl_xor_sync(0xffffffff,rg_k[1] , 0x4 );
    rg_k[3]  = __shfl_xor_sync(0xffffffff,rg_k[3] , 0x4 );
    if(warpLane&0x4 ) SWP_KEY(K, rg_k[0] , rg_k[1] );
    if(warpLane&0x4 ) SWP_KEY(K, rg_k[2] , rg_k[3] );
    if(warpLane&0x8 ) SWP_KEY(K, rg_k[0] , rg_k[2] );
    if(warpLane&0x8 ) SWP_KEY(K, rg_k[1] , rg_k[3] );
    rg_k[2]  = __shfl_xor_sync(0xffffffff,rg_k[2] , 0x8 );
    rg_k[3]  = __shfl_xor_sync(0xffffffff,rg_k[3] , 0x8 );
    if(warpLane&0x8 ) SWP_KEY(K, rg_k[0] , rg_k[2] );
    if(warpLane&0x8 ) SWP_KEY(K, rg_k[1] , rg_k[3] );
    if(warpLane&0x10) SWP_KEY(K, rg_k[0] , rg_k[1] );
    if(warpLane&0x10) SWP_KEY(K, rg_k[2] , rg_k[3] );
    rg_k[1]  = __shfl_xor_sync(0xffffffff,rg_k[1] , 0x10);
    rg_k[3]  = __shfl_xor_sync(0xffffffff,rg_k[3] , 0x10);
    if(warpLane&0x10) SWP_KEY(K, rg_k[0] , rg_k[1] );
    if(warpLane&0x10) SWP_KEY(K, rg_k[2] , rg_k[3] );

    SWP_KEY(K, rg_k[1] , rg_k[2] );
}


/****************************************************************
 * @brief extract kmers from sequence
 *
 * @details each thread processes 4 chars and extracts up to 4 kmers
 *
 * @param sequenceLength max length is 128 (32*4)
 */
 template<
    class SizeT>
__device__ __inline__
void warp_kmerize(
    const char * sequenceBegin,
    SizeT sequenceLength,
    numk_t kmerSize,
    feature items[4]
)
{
    const int warpLane = threadIdx.x % WARPSIZE;

    constexpr int subWarpSize = 4;
    const int lane = subWarpSize-1 - warpLane % subWarpSize;

    uint32_t seq = 0;
    uint32_t ambig = 0;

    constexpr int charsPerThread = 4;
    constexpr int charsPerSubWarp = subWarpSize*charsPerThread;
    char4 chars = reinterpret_cast<const char4*>(sequenceBegin)[warpLane];

    char c;
    c = (charsPerThread*warpLane+0 < sequenceLength) ? chars.x : 'N';
    c &= 0b11011111; // to upper case

    seq |= (c == 'A') ? 0 : 0;
    seq |= (c == 'C') ? 1 : 0;
    seq |= (c == 'G') ? 2 : 0;
    seq |= (c == 'T') ? 3 : 0;
    ambig |= (c == 'A' || c == 'C' || c == 'G' || c == 'T') ? 0 : 1;

    seq   <<= 2;
    ambig <<= 1;

    c = (charsPerThread*warpLane+1 < sequenceLength) ? chars.y : 'N';
    c &= 0b11011111; // to upper case

    seq |= (c == 'A') ? 0 : 0;
    seq |= (c == 'C') ? 1 : 0;
    seq |= (c == 'G') ? 2 : 0;
    seq |= (c == 'T') ? 3 : 0;
    ambig |= (c == 'A' || c == 'C' || c == 'G' || c == 'T') ? 0 : 1;

    seq   <<= 2;
    ambig <<= 1;

    c = (charsPerThread*warpLane+2 < sequenceLength) ? chars.z : 'N';
    c &= 0b11011111; // to upper case

    seq |= (c == 'A') ? 0 : 0;
    seq |= (c == 'C') ? 1 : 0;
    seq |= (c == 'G') ? 2 : 0;
    seq |= (c == 'T') ? 3 : 0;
    ambig |= (c == 'A' || c == 'C' || c == 'G' || c == 'T') ? 0 : 1;

    seq   <<= 2;
    ambig <<= 1;

    c = (charsPerThread*warpLane+3 < sequenceLength) ? chars.w : 'N';
    c &= 0b11011111; // to upper case

    seq |= (c == 'A') ? 0 : 0;
    seq |= (c == 'C') ? 1 : 0;
    seq |= (c == 'G') ? 2 : 0;
    seq |= (c == 'T') ? 3 : 0;
    ambig |= (c == 'A' || c == 'C' || c == 'G' || c == 'T') ? 0 : 1;

    seq <<= 2*charsPerThread*lane;
    ambig <<= charsPerThread*lane;

    // combine threads of same subwarp
    #pragma unroll
    for(int stride = 1; stride < subWarpSize; stride *= 2) {
        seq   |= __shfl_xor_sync(0xFFFFFFFF, seq, stride);
        ambig |= __shfl_xor_sync(0xFFFFFFFF, ambig, stride);
    }

    // combine subwarp with subsequent subwarp
    uint64_t seqStore = __shfl_down_sync(0xFFFFFFFF, seq, subWarpSize);
    seqStore |= (uint64_t(seq) << 32) | seqStore;

    uint32_t ambigStore = __shfl_down_sync(0xFFFFFFFF, ambig, subWarpSize);
    // bounds checking: last subwarp adds ~0
    ambigStore = (warpLane < WARPSIZE - subWarpSize) ? ambigStore : uint32_t(~0);
    ambigStore = (ambig << 16) | ambigStore;

    // if(lane == 0) {
        // const int subwarpId = warpLane / subWarpSize;
        // printf("id: %d seq: %u ambig: %u\n", subwarpId, seq, ambig);
    // }

    // letters stored from high to low bits
    // mask get lowest bits
    constexpr uint8_t kmerBits = sizeof(kmer_type)*CHAR_BIT;
    const kmer_type kmerMask = kmer_type(~0) >> (kmerBits - 2*kmerSize);
    // mask get highest bits
    const uint32_t ambigMask = uint32_t(~0) << (32 - kmerSize);

    // each thread extracts one feature
    #pragma unroll
    for(int i = 0; i < 4; ++i) {
        const uint8_t seqSlot  = (2*i + warpLane / charsPerSubWarp);
        const uint8_t seqSrc   = seqSlot * 4;
        const uint8_t kmerSlot = warpLane % charsPerSubWarp;

        uint32_t ambig = __shfl_sync(0xFFFFFFFF, ambigStore, seqSrc);
        // shift to highest bits
        ambig <<= kmerSlot;

        uint64_t seq = __shfl_sync(0xFFFFFFFF, seqStore, seqSrc);
        // shift kmer to lowest bits
        seq >>= 64 - 2*(kmerSlot + kmerSize);

        // get lowest bits
        kmer_type kmer = kmer_type(seq) & kmerMask;

        kmer = make_canonical_2bit(kmer);

        bool pred = !(ambig & ambigMask);
        items[i] = pred ? sketching_hash{}(kmer) : feature(~0);
    }
}





/****************************************************************
 *
 * @brief insert batch of sequences into hashtable
 *
 * @details each block processes one window of at most 128 characters
 *
 */
template<
    int BLOCK_THREADS,
    int MAX_SKETCH_SIZE,
    class Hashtable,
    class IndexT,
    class SizeT>
__global__
__launch_bounds__(BLOCK_THREADS)
void insert_features(
    Hashtable hashtable,
    IndexT numTargets,
    const target_id * targetIds,
    const window_id * winOffsets,
    const char      * sequence,
    const SizeT     * sequenceOffsets,
    numk_t kmerSize,
    uint32_t maxSketchSize,
    uint32_t windowSize,
    uint32_t windowStride
)
{
    using index_type = IndexT;
    using size_type = SizeT;

    constexpr int warpsPerBlock = BLOCK_THREADS / WARPSIZE;
    constexpr int itemsPerThread = 4; // 32*4=128 items max for warp sort

    __shared__ feature sketch[warpsPerBlock][MAX_SKETCH_SIZE];

    for(index_type tgt = blockIdx.y; tgt < numTargets; tgt += gridDim.y) {
        const size_type sequenceBegin = sequenceOffsets[tgt];
        const size_type sequenceEnd   = sequenceOffsets[tgt+1];
        const size_type sequenceLength = sequenceEnd - sequenceBegin;
        const window_id lastWindowSize = sequenceLength % windowStride;
        const window_id numWindows = sequenceLength / windowStride + (lastWindowSize >= kmerSize);

        const int warpId = threadIdx.x / WARPSIZE;

        //each block processes one window
        for(window_id win = blockIdx.x * warpsPerBlock + warpId;
            win < numWindows;
            win += warpsPerBlock * gridDim.x)
        {
            const size_type offsetBegin = sequenceBegin + win*windowStride;
            const size_type offsetEnd = min(sequenceEnd, offsetBegin + windowSize);
            const size_type thisWindowSize = offsetEnd - offsetBegin;
            const char * windowBegin = sequence + offsetBegin;

            feature items[itemsPerThread];
            warp_kmerize(windowBegin, thisWindowSize, kmerSize, items);

            warp_sort_128(items);

            uint32_t sketchCounter =
                unique_sketch<itemsPerThread>(items, sketch[warpId], maxSketchSize);

            __syncwarp();

            const target_id targetId = targetIds[tgt];
            const window_id winOffset = winOffsets[tgt]+win ;

            insert_into_hashtable(
                hashtable, sketch[warpId], sketchCounter, targetId, winOffset);

            __syncwarp();
        }
    }
}



//---------------------------------------------------------------
template<
    int BLOCK_THREADS,
    int MAX_SKETCH_SIZE,
    class Hashtable,
    class SizeT,
    class BucketSizeT,
    class Location>
__global__
__launch_bounds__(BLOCK_THREADS)
void gpu_hahstable_query(
    Hashtable hashtable,
    uint32_t numQueries,
    const SizeT * sequenceOffsets,
    const char * sequences,
    feature * sketches,
    numk_t kmerSize,
    uint32_t maxSketchSize,
    uint32_t windowSize,
    uint32_t windowStride,
    const Location * locations,
    BucketSizeT maxLocationsPerFeature,
    Location * locations_out,
    int      * locationCounts
)
{
    using size_type = SizeT;
    using location = Location;

    constexpr int warpsPerBlock = BLOCK_THREADS / WARPSIZE;
    constexpr int itemsPerThread = 4; // 32*4=128 items max for warp sort

    // using uint64_t enables to store either "feature" or "Hashtable::value_type"
    __shared__ uint64_t sketch[warpsPerBlock][MAX_SKETCH_SIZE];

    const int warpId = threadIdx.x / WARPSIZE;
    const int warpLane = threadIdx.x % WARPSIZE;

    //each block processes one query (= one window)
    for(int queryId = blockIdx.x * warpsPerBlock + warpId;
        queryId < numQueries;
        queryId += warpsPerBlock * gridDim.x)
    {
        const size_type sequenceLength = sequenceOffsets[queryId+1] - sequenceOffsets[queryId];

        //only process non-empty queries
        if(sequenceLength) {
            const char * sequenceBegin = sequences + sequenceOffsets[queryId];

            feature items[itemsPerThread];
            warp_kmerize(sequenceBegin, sequenceLength, kmerSize, items);

            warp_sort_128(items);

            uint32_t realSketchSize = unique_sketch<itemsPerThread>(
                items, sketch[warpId], maxSketchSize);

            __syncwarp();

            // if(warpLane == 0) {
            //     printf("query %d features %llu %llu %llu %llu %llu %llu %llu %llu ...\n", queryId,
            //         sketch[warpId][0],
            //         sketch[warpId][1],
            //         sketch[warpId][2],
            //         sketch[warpId][3],
            //         sketch[warpId][4],
            //         sketch[warpId][5],
            //         sketch[warpId][6],
            //         sketch[warpId][7]);
            // }

            if(sketches != nullptr) {
                if(warpLane < maxSketchSize)
                    sketches[queryId*maxSketchSize + warpLane] =
                        (warpLane < realSketchSize) ?
                        sketch[warpId][warpLane] :
                        feature(~0);
            }

            query_hashtable(hashtable, sketch[warpId], realSketchSize);

            __syncwarp();

            location * out = locations_out + queryId*maxSketchSize*maxLocationsPerFeature;

            uint32_t numLocations = copy_loctions(
                sketch[warpId], realSketchSize, locations, maxLocationsPerFeature, out);

            if(warpLane == 0) {
                //store number of results of this query
                locationCounts[queryId] = numLocations;
            }

            __syncwarp();
        }
        else {
            if(sketches != nullptr) {
                if(warpLane < maxSketchSize)
                    sketches[queryId*maxSketchSize + warpLane] = feature(~0);
            }

            if(warpLane == 0) {
                //store number of results of this query
                locationCounts[queryId] = 0;
            }
        }
    }
}


//---------------------------------------------------------------
template<
    int BLOCK_THREADS,
    int MAX_SKETCH_SIZE,
    class Hashtable,
    class BucketSizeT,
    class Location>
__global__
__launch_bounds__(BLOCK_THREADS)
void gpu_hahstable_query(
    Hashtable hashtable,
    uint32_t numQueries,
    const feature * sketches,
    numk_t kmerSize,
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

    constexpr int warpsPerBlock = BLOCK_THREADS / WARPSIZE;

    // using uint64_t enables to store either "feature" or "Hashtable::value_type"
    __shared__ uint64_t sketch[warpsPerBlock][MAX_SKETCH_SIZE];

    const int warpId = threadIdx.x / WARPSIZE;
    const int warpLane = threadIdx.x % WARPSIZE;

    //each block processes one query (= one window)
    for(int queryId = blockIdx.x * warpsPerBlock + warpId;
        queryId < numQueries;
        queryId += warpsPerBlock * gridDim.x)
    {
        const feature f = (warpLane < maxSketchSize) ?
                          sketches[queryId*maxSketchSize + warpLane] :
                          feature(~0);

        const bool validFeature = f != feature(~0);

        const uint32_t sketchMask = __ballot_sync(0xFFFFFFFF, validFeature);

        const uint32_t realSketchSize = __popc(sketchMask);

        //only process non-empty queries
        if(realSketchSize) {

            if(warpLane < realSketchSize)
                sketch[warpId][warpLane] = f;

            __syncwarp();

            // if(warpLane == 0) {
            //     printf("query %d features %llu %llu %llu %llu %llu %llu %llu %llu ...\n",
            //         queryId,
            //         sketch[warpId][0],
            //         sketch[warpId][1],
            //         sketch[warpId][2],
            //         sketch[warpId][3],
            //         sketch[warpId][4],
            //         sketch[warpId][5],
            //         sketch[warpId][6],
            //         sketch[warpId][7]);
            // }

            query_hashtable(hashtable, sketch[warpId], realSketchSize);

            __syncwarp();

            location * out = locations_out + queryId*maxSketchSize*maxLocationsPerFeature;

            uint32_t numLocations = copy_loctions(
                sketch[warpId], realSketchSize, locations, maxLocationsPerFeature, out);

            if(warpLane == 0) {
                //store number of results of this query
                locationCounts[queryId] = numLocations;
            }

            __syncwarp();
        }
        else {
            if(warpLane == 0) {
                //store number of results of this query
                locationCounts[queryId] = 0;
            }
        }
    }
}


//---------------------------------------------------------------
template<
    int BLOCK_THREADS,
    int MAX_SKETCH_SIZE,
    class Hashtable,
    class SizeT,
    class BucketSizeT,
    class Location>
__global__
__launch_bounds__(BLOCK_THREADS)
void gpu_hahstable_query(
    Hashtable hashtable,
    uint32_t numQueries,
    const SizeT * sequenceOffsets,
    const char * sequences,
    feature * sketches,
    numk_t kmerSize,
    uint32_t maxSketchSize,
    uint32_t windowSize,
    uint32_t windowStride,
    BucketSizeT maxLocationsPerFeature,
    Location * locations_out,
    int      * locationCounts
)
{
    using size_type = SizeT;
    using location = Location;

    constexpr int warpsPerBlock = BLOCK_THREADS / WARPSIZE;
    constexpr int itemsPerThread = 4; // 32*4=128 items max for warp sort

    __shared__ feature sketch[warpsPerBlock][MAX_SKETCH_SIZE];
    __shared__ uint32_t sizes[warpsPerBlock][MAX_SKETCH_SIZE+1];

    const int warpId = threadIdx.x / WARPSIZE;
    const int warpLane = threadIdx.x % WARPSIZE;

    //each block processes one query (= one window)
    for(int queryId = blockIdx.x * warpsPerBlock + warpId;
        queryId < numQueries;
        queryId += warpsPerBlock * gridDim.x)
    {
        const size_type sequenceLength = sequenceOffsets[queryId+1] - sequenceOffsets[queryId];

        //only process non-empty queries
        if(sequenceLength) {
            const char * sequenceBegin = sequences + sequenceOffsets[queryId];

            feature items[itemsPerThread];
            warp_kmerize(sequenceBegin, sequenceLength, kmerSize, items);

            warp_sort_128(items);

            uint32_t realSketchSize = unique_sketch<itemsPerThread>(
                items, sketch[warpId], maxSketchSize);

            __syncwarp();

            // if(warpLane == 0) {
            //     printf("query %d features %u %u %u %u %u %u %u %u ...\n", queryId,
            //         sketch[warpId][0],
            //         sketch[warpId][1],
            //         sketch[warpId][2],
            //         sketch[warpId][3],
            //         sketch[warpId][4],
            //         sketch[warpId][5],
            //         sketch[warpId][6],
            //         sketch[warpId][7]);
            // }

            if(sketches != nullptr) {
                if(warpLane < maxSketchSize)
                    sketches[queryId*maxSketchSize + warpLane] =
                        (warpLane < realSketchSize) ?
                        sketch[warpId][warpLane] :
                        feature(~0);
            }

            location * out = locations_out + queryId*maxSketchSize*maxLocationsPerFeature;

            uint32_t numLocations = query_bucket_hashtable<MAX_SKETCH_SIZE>(
                hashtable, sketch[warpId], sizes[warpId], realSketchSize, out);

            if(warpLane == 0) {
                //store number of results of this query
                locationCounts[queryId] = numLocations;
            }

            __syncwarp();
        }
        else {
            if(sketches != nullptr) {
                if(warpLane < maxSketchSize)
                    sketches[queryId*maxSketchSize + warpLane] = feature(~0);
            }

            if(warpLane == 0) {
                //store number of results of this query
                locationCounts[queryId] = 0;
            }
        }
    }
}


//---------------------------------------------------------------
template<
    int BLOCK_THREADS,
    int MAX_SKETCH_SIZE,
    class Hashtable,
    class BucketSizeT,
    class Location>
__global__
__launch_bounds__(BLOCK_THREADS)
void gpu_hahstable_query(
    Hashtable hashtable,
    uint32_t numQueries,
    const feature * sketches,
    numk_t kmerSize,
    uint32_t maxSketchSize,
    uint32_t windowSize,
    uint32_t windowStride,
    BucketSizeT maxLocationsPerFeature,
    Location * locations_out,
    int      * locationCounts
)
{
    using location = Location;

    constexpr int warpsPerBlock = BLOCK_THREADS / WARPSIZE;

    __shared__ feature sketch[warpsPerBlock][MAX_SKETCH_SIZE];
    __shared__ uint32_t sizes[warpsPerBlock][MAX_SKETCH_SIZE+1];

    const int warpId = threadIdx.x / WARPSIZE;
    const int warpLane = threadIdx.x % WARPSIZE;

    //each block processes one query (= one window)
    for(int queryId = blockIdx.x * warpsPerBlock + warpId;
        queryId < numQueries;
        queryId += warpsPerBlock * gridDim.x)
    {
        const feature f = (warpLane < maxSketchSize) ?
                          sketches[queryId*maxSketchSize + warpLane] :
                          feature(~0);

        const bool validFeature = f != feature(~0);

        const uint32_t sketchMask = __ballot_sync(0xFFFFFFFF, validFeature);

        const uint32_t realSketchSize = __popc(sketchMask);

        //only process non-empty queries
        if(realSketchSize) {

            if(warpLane < realSketchSize)
                sketch[warpId][warpLane] = f;

            __syncwarp();

            // if(warpLane == 0) {
            //     printf("query %d features %u %u %u %u %u %u %u %u ...\n",
            //         queryId,
            //         sketch[warpId][0],
            //         sketch[warpId][1],
            //         sketch[warpId][2],
            //         sketch[warpId][3],
            //         sketch[warpId][4],
            //         sketch[warpId][5],
            //         sketch[warpId][6],
            //         sketch[warpId][7]);
            // }

            location * out = locations_out + queryId*maxSketchSize*maxLocationsPerFeature;

            uint32_t numLocations = query_bucket_hashtable<MAX_SKETCH_SIZE>(
                hashtable, sketch[warpId], sizes[warpId], realSketchSize, out);

            if(warpLane == 0) {
                //store number of results of this query
                locationCounts[queryId] = numLocations;
            }

            __syncwarp();
        }
        else {
            if(warpLane == 0) {
                //store number of results of this query
                locationCounts[queryId] = 0;
            }
        }
    }
}



} // namespace mc

#endif

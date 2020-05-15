#ifndef MC_GPU_HASHMAP_OPERATIONS_H_
#define MC_GPU_HASHMAP_OPERATIONS_H_

#include "config.h"
#include "dna_encoding.h"

#include "../dep/cudahelpers/cuda_helpers.cuh"
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

    for(uint32_t i = (threadIdx.x % WARPSIZE); i < sketchCounter*Hashtable::cg_size(); i += WARPSIZE) {
        const uint8_t gid = i / Hashtable::cg_size();

        const auto status = hashtable.insert(
            sketch[gid], location{winOffset, targetId}, group);

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

    for(uint32_t i = (threadIdx.x % WARPSIZE); i < sketchCounter*Hashtable::cg_size(); i += WARPSIZE) {
        const uint8_t gid = i / Hashtable::cg_size();
        typename Hashtable::value_type valuesOffset = 0;
        //if key not found valuesOffset stays 0
        const auto status = hashtable.retrieve(
            sketch[gid], valuesOffset, group);

        if(group.thread_rank() == 0) {
            // if(status.has_any()) {
                // printf("status %d\n", status.has_any());
                // printf("status %d\n", status.has_key_not_found());
            // }

            sketch[gid] = valuesOffset;
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
    using bucket_size_type = BucketSizeT;

    uint32_t totalLocations = 0;

    for(uint32_t s = 0; s < sketchSize; ++s) {
        auto bucketOffset = sketch[s];
        bucket_size_type bucketSize = min(bucket_size_type(bucketOffset), maxLocationsPerFeature);
        bucketOffset >>= sizeof(bucket_size_type)*CHAR_BIT;

        //copy locations
        for(uint32_t i = (threadIdx.x % WARPSIZE); i < bucketSize; i += WARPSIZE) {
            out[totalLocations + i] = locations[bucketOffset + i];
        }

        totalLocations += bucketSize;
    }

    return totalLocations;
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
// BLOCK_THREADS has to be 32
template<
    int BLOCK_THREADS = 32,
    int ITEMS_PER_THREAD = 4,
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

    __shared__ union {
        //TODO MAX_SKETCH_SIZE
        feature sketch[warpsPerBlock*16];
    } tempStorage;

    for(index_type tgt = blockIdx.y; tgt < numTargets; tgt += gridDim.y) {
        const size_type sequenceLength = sequenceOffsets[tgt+1] - sequenceOffsets[tgt];
        const window_id lastWindowSize = sequenceLength % windowStride;
        const window_id numWindows = sequenceLength / windowStride + (lastWindowSize >= kmerSize);

        const int warpId = threadIdx.x / WARPSIZE;

        //each block processes one window
        for(window_id win = blockIdx.x * warpsPerBlock + warpId; win < numWindows; win += warpsPerBlock * gridDim.x)
        {
            const size_type offsetBegin = sequenceOffsets[tgt] + win*windowStride;
            const size_type offsetEnd = min(sequenceOffsets[tgt+1], offsetBegin + windowSize);
            const size_type thisWindowSize = offsetEnd - offsetBegin;
            const char * windowBegin = sequence + offsetBegin;

            feature items[ITEMS_PER_THREAD];
            warp_kmerize(windowBegin, thisWindowSize, kmerSize, items);

            warp_sort_128(items);

            uint32_t sketchCounter =
                unique_sketch<ITEMS_PER_THREAD>(items, tempStorage.sketch+warpId*16, maxSketchSize);

            __syncwarp();

            const target_id targetId = targetIds[tgt];
            const window_id winOffset = winOffsets[tgt]+win ;

            insert_into_hashtable(hashtable, tempStorage.sketch+warpId*16, sketchCounter, targetId, winOffset);

            __syncwarp();
        }
    }
}



//---------------------------------------------------------------
// BLOCK_THREADS has to be 32
template<
    int BLOCK_THREADS,
    int ITEMS_PER_THREAD,
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

    __shared__ union {
        //TODO MAX_SKETCH_SIZE
        uint64_t sketch[warpsPerBlock*16];
    } tempStorage;

    const int warpId = threadIdx.x / WARPSIZE;
    const int warpLane = threadIdx.x % WARPSIZE;

    //each block processes one query (= one window)
    for(int queryId = blockIdx.x * warpsPerBlock + warpId; queryId < numQueries; queryId += warpsPerBlock * gridDim.x) {
        const size_type sequenceLength = sequenceOffsets[queryId+1] - sequenceOffsets[queryId];

        //only process non-empty queries
        if(sequenceLength) {
            const char * sequenceBegin = sequences + sequenceOffsets[queryId];

            feature items[ITEMS_PER_THREAD];
            warp_kmerize(sequenceBegin, sequenceLength, kmerSize, items);

            warp_sort_128(items);

            uint32_t realSketchSize =
                unique_sketch<ITEMS_PER_THREAD>(items, tempStorage.sketch+warpId*16, maxSketchSize);

            __syncwarp();

            // if(warpLane == 0) {
            //     printf("query %d features %llu %llu %llu %llu %llu %llu %llu %llu ...\n", queryId,
            //         tempStorage.sketch[warpId*16 + 0],
            //         tempStorage.sketch[warpId*16 + 1],
            //         tempStorage.sketch[warpId*16 + 2],
            //         tempStorage.sketch[warpId*16 + 3],
            //         tempStorage.sketch[warpId*16 + 4],
            //         tempStorage.sketch[warpId*16 + 5],
            //         tempStorage.sketch[warpId*16 + 6],
            //         tempStorage.sketch[warpId*16 + 7]);
            // }

            if(sketches != nullptr) {
                if(warpLane < maxSketchSize)
                    sketches[queryId*maxSketchSize + warpLane] =
                        (warpLane < realSketchSize) ?
                        tempStorage.sketch[warpId*16 + warpLane] :
                        feature(~0);
            }

            query_hashtable(hashtable, tempStorage.sketch+warpId*16, realSketchSize);

            __syncwarp();

            location * out = locations_out + queryId*maxSketchSize*maxLocationsPerFeature;

            uint32_t numLocations = copy_loctions(
                tempStorage.sketch+warpId*16, realSketchSize, locations, maxLocationsPerFeature, out);

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
// BLOCK_THREADS has to be 32
template<
    int BLOCK_THREADS,
    int ITEMS_PER_THREAD,
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

    __shared__ union {
        //TODO MAX_SKETCH_SIZE
        uint64_t sketch[warpsPerBlock*16];
    } tempStorage;

    const int warpId = threadIdx.x / WARPSIZE;
    const int warpLane = threadIdx.x % WARPSIZE;

    //each block processes one query (= one window)
    for(int queryId = blockIdx.x * warpsPerBlock + warpId; queryId < numQueries; queryId += warpsPerBlock * gridDim.x)
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
                tempStorage.sketch[warpId*16 + warpLane] = f;

            __syncwarp();

            // if(warpLane == 0) {
            //     printf("query %d features %llu %llu %llu %llu %llu %llu %llu %llu ...\n", queryId,
            //         tempStorage.sketch[warpId*16 + 0],
            //         tempStorage.sketch[warpId*16 + 1],
            //         tempStorage.sketch[warpId*16 + 2],
            //         tempStorage.sketch[warpId*16 + 3],
            //         tempStorage.sketch[warpId*16 + 4],
            //         tempStorage.sketch[warpId*16 + 5],
            //         tempStorage.sketch[warpId*16 + 6],
            //         tempStorage.sketch[warpId*16 + 7]);
            // }

            query_hashtable(hashtable, tempStorage.sketch+warpId*16, realSketchSize);

            __syncwarp();

            location * out = locations_out + queryId*maxSketchSize*maxLocationsPerFeature;

            uint32_t numLocations = copy_loctions(
                tempStorage.sketch+warpId*16, realSketchSize, locations, maxLocationsPerFeature, out);

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

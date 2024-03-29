/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2022 Robin Kobus  (kobus@uni-mainz.de)
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

#ifndef MC_GPU_HASHMAP_OPERATIONS_H_
#define MC_GPU_HASHMAP_OPERATIONS_H_

#include "config.h"
#include "dna_encoding.h"

#include "../dep/hpc_helpers/include/cuda_helpers.cuh"
#include "../dep/bb_segsort/src/bb_exch_keys.cuh"

namespace mc {


/**************************************************************************//**
 *
 * @brief extract kmers from sequence
 *
 * @details each thread processes 4 characters and extracts up to 4 kmers
 *
 * @param sequenceBegin pointer to sequence begin
 * @param sequenceLength max length is 128 (32*4)
 * @param kmerSize size of kemrs to be extracted
 * @param items array to store 4 kmers per thread
 *
 *****************************************************************************/
 template<
    class Offset>
__device__ __inline__
void warp_kmerize(
    const char * sequenceBegin,
    Offset sequenceLength,
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

    // combine threads of same subwarp,
    // so each thread has charsPerThread*subWarpSize characters
    #pragma unroll
    for (int stride = 1; stride < subWarpSize; stride *= 2) {
        seq   |= __shfl_xor_sync(0xFFFFFFFF, seq, stride);
        ambig |= __shfl_xor_sync(0xFFFFFFFF, ambig, stride);
    }

    // combine subwarp with subsequent subwarp,
    // so each thread has 2*charsPerThread*subWarpSize characters
    uint64_t seqStore = __shfl_down_sync(0xFFFFFFFF, seq, subWarpSize);
    seqStore |= (uint64_t(seq) << 32) | seqStore;

    uint32_t ambigStore = __shfl_down_sync(0xFFFFFFFF, ambig, subWarpSize);
    // bounds checking: last subwarp adds ~0
    ambigStore = (warpLane < WARPSIZE - subWarpSize) ? ambigStore : (uint32_t(~0) >> 16);
    ambigStore = (ambig << 16) | ambigStore;

    // letters stored from high to low bits
    // thread i uses characters [charsPerThread*i, charsPerThread*(i+1)+kmerSize-1]
    // shift to low bits,
    // so each thread uses lowest kmerSize+charsPerThread-1 characters
    seqStore >>= 2*charsPerThread*lane + 2 * (17 - kmerSize);
    ambigStore >>= charsPerThread*lane + 1 * (17 - kmerSize);

    // mask get lowest kmerSize characters
    constexpr uint8_t kmerBits = sizeof(kmer_type)*CHAR_BIT;
    const kmer_type kmerMask = kmer_type(~0) >> (kmerBits - 2*kmerSize);
    const uint32_t ambigMask = uint32_t(~0) >> (32 - kmerSize);

    // each thread generates charsPerThread kmers
    #pragma unroll
    for (int i = 0; i < charsPerThread; ++i) {
        // shift by one character each iteration
        ambig = ambigStore >> (i);
        seq = seqStore >> (2*i);

        // get lowest bits
        kmer_type kmer = kmer_type(seq) & kmerMask;

        kmer = make_canonical_2bit(kmer, kmerSize);

        // check for ambiguous base
        bool pred = !(ambig & ambigMask);

        // ignore kmers with ambiguous base
        items[i] = pred ? sketching_hash{}(kmer) : feature(~0);
    }
}


/**************************************************************************//**
 *
 * @brief in-place bitonic sort of 128 (4*32) elements
 *
 * @details see bb_segsort submodule
 *
 * @param rg_k array of 4 items per thread
 *
 *****************************************************************************/
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
    if (warpLane&0x1 ) SWP_KEY(K, rg_k[0] , rg_k[1] );
    if (warpLane&0x1 ) SWP_KEY(K, rg_k[2] , rg_k[3] );
    rg_k[1]  = __shfl_xor_sync(0xffffffff,rg_k[1] , 0x1 );
    rg_k[3]  = __shfl_xor_sync(0xffffffff,rg_k[3] , 0x1 );
    if (warpLane&0x1 ) SWP_KEY(K, rg_k[0] , rg_k[1] );
    if (warpLane&0x1 ) SWP_KEY(K, rg_k[2] , rg_k[3] );
    if (warpLane&0x2 ) SWP_KEY(K, rg_k[0] , rg_k[2] );
    if (warpLane&0x2 ) SWP_KEY(K, rg_k[1] , rg_k[3] );
    rg_k[2]  = __shfl_xor_sync(0xffffffff,rg_k[2] , 0x2 );
    rg_k[3]  = __shfl_xor_sync(0xffffffff,rg_k[3] , 0x2 );
    if (warpLane&0x2 ) SWP_KEY(K, rg_k[0] , rg_k[2] );
    if (warpLane&0x2 ) SWP_KEY(K, rg_k[1] , rg_k[3] );
    if (warpLane&0x4 ) SWP_KEY(K, rg_k[0] , rg_k[1] );
    if (warpLane&0x4 ) SWP_KEY(K, rg_k[2] , rg_k[3] );
    rg_k[1]  = __shfl_xor_sync(0xffffffff,rg_k[1] , 0x4 );
    rg_k[3]  = __shfl_xor_sync(0xffffffff,rg_k[3] , 0x4 );
    if (warpLane&0x4 ) SWP_KEY(K, rg_k[0] , rg_k[1] );
    if (warpLane&0x4 ) SWP_KEY(K, rg_k[2] , rg_k[3] );
    if (warpLane&0x8 ) SWP_KEY(K, rg_k[0] , rg_k[2] );
    if (warpLane&0x8 ) SWP_KEY(K, rg_k[1] , rg_k[3] );
    rg_k[2]  = __shfl_xor_sync(0xffffffff,rg_k[2] , 0x8 );
    rg_k[3]  = __shfl_xor_sync(0xffffffff,rg_k[3] , 0x8 );
    if (warpLane&0x8 ) SWP_KEY(K, rg_k[0] , rg_k[2] );
    if (warpLane&0x8 ) SWP_KEY(K, rg_k[1] , rg_k[3] );
    if (warpLane&0x10) SWP_KEY(K, rg_k[0] , rg_k[1] );
    if (warpLane&0x10) SWP_KEY(K, rg_k[2] , rg_k[3] );
    rg_k[1]  = __shfl_xor_sync(0xffffffff,rg_k[1] , 0x10);
    rg_k[3]  = __shfl_xor_sync(0xffffffff,rg_k[3] , 0x10);
    if (warpLane&0x10) SWP_KEY(K, rg_k[0] , rg_k[1] );
    if (warpLane&0x10) SWP_KEY(K, rg_k[2] , rg_k[3] );

    SWP_KEY(K, rg_k[1] , rg_k[2] );
}


/**************************************************************************//**
 *
 * @brief extract unique features from sorted items
 *
 * @param items         sorted features
 * @param sketch        pointer to sketch array
 * @param maxSketchSize maximum allowed sketch size
 *
 * @return number of unique features stored in sketch
 *
 *****************************************************************************/
template<
    int ITEMS_PER_THREAD,
    class Feature,
    class Feature2>
__device__ __inline__
uint32_t unique_sketch(
    const Feature items[ITEMS_PER_THREAD],
    Feature2 * sketch,
    uint32_t maxSketchSize)
{
    const int warpLane = threadIdx.x % WARPSIZE;

    uint32_t uniqueFeatureCounter = 0;

    #pragma unroll
    for (int i=0; i<ITEMS_PER_THREAD; ++i)
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
        if (predicate && (sketchPos < maxSketchSize)) {
            sketch[sketchPos] = kmer;
            // printf("sketchPos: %2d feature: %d\n", sketchPos, kmer);
        }

        uniqueFeatureCounter += count;
    }
    // cap counter
    uniqueFeatureCounter = min(uniqueFeatureCounter, maxSketchSize);

    return uniqueFeatureCounter;
}


/**************************************************************************//**
 *
 * @brief use one warp to generate min-hasher sketch from sequence;
 *
 * @details smallest *unique* hash values of *one* hash function
 *
 * @param sequenceBegin pointer to sequence begin
 * @param sequenceEnd   pointer to sequence end
 * @param kmerSize      size of kemrs to be extracted
 * @param sketch        pointer to sketch storage
 * @param maxSketchSize maximum allowed sketch size
 *
 * @return number of features in sketch
 *
 *****************************************************************************/
 template<class Offset, class Feature>
__device__ __inline__
uint32_t warp_make_sketch(
    const char * sequenceBegin,
    const char * sequenceEnd,
    const numk_t kmerSize,
    Feature * sketch,
    const uint32_t maxSketchSize)
{
    using offset_type = Offset;
    const int warpLane = threadIdx.x % WARPSIZE;

    uint32_t sketchSize = 0;

    const offset_type sequenceLength = sequenceEnd - sequenceBegin;

    // skip empty queries
    if (!sequenceLength) return sketchSize;

    constexpr int itemsPerThread = 4; // 32*4=128 items
    feature items[itemsPerThread];

    if (sequenceLength <= 128) {
        warp_kmerize(sequenceBegin, sequenceLength, kmerSize, items);

        warp_sort_128(items);

        sketchSize = unique_sketch<itemsPerThread>(items, sketch, maxSketchSize);

        __syncwarp();
    }
    else {
        // set stride so that maxSketchSize items remain free to store previous sketch
        uint32_t paddedSketchSize = (maxSketchSize+3)/4*4;
        offset_type stride = 4*WARPSIZE - paddedSketchSize;
        offset_type maxLength = stride + kmerSize - 1;

        if (warpLane < paddedSketchSize)
            sketch[warpLane] = feature(~0);
        __syncwarp();

        for (; sequenceBegin < sequenceEnd; sequenceBegin += stride) {
            const offset_type sequenceLength = std::min(sequenceEnd - sequenceBegin, long(maxLength));

            // only process non-empty queries
            if (sequenceLength >= kmerSize) {
                // limited length leads to unused items in hindmost threads
                warp_kmerize(sequenceBegin, sequenceLength, kmerSize, items);

                // load previous sketch in hindmost threads
                if (warpLane-int(WARPSIZE-paddedSketchSize/4) >= 0) {
                    int pos = warpLane-(WARPSIZE-paddedSketchSize/4);
                    items[0] = sketch[4*pos];
                    items[1] = sketch[4*pos+1];
                    items[2] = sketch[4*pos+2];
                    items[3] = sketch[4*pos+3];
                }

                warp_sort_128(items);

                sketchSize = unique_sketch<itemsPerThread>(items, sketch, maxSketchSize);
            }

            __syncwarp();
        }
    }

    return sketchSize;
}


/**************************************************************************//**
 *
 * @brief use one warp to load sketch from memory
 *
 * @param sketches       pointer to sketch in memory
 * @param sketch         pointer to sketch storage
 * @param maxSketchSize  maximum allowed sketch size
 *
 * @return number of features in sketch
 *
 *****************************************************************************/
 template<class Feature>
__device__ __inline__
uint32_t warp_load_sketch(
    const feature * sketches,
    Feature * sketch,
    const uint32_t maxSketchSize)
{
    const int warpLane = threadIdx.x % WARPSIZE;

    const feature f =
        (warpLane < maxSketchSize) ? sketches[warpLane] : feature(~0);

    // count valid features
    const bool validFeature = f != feature(~0);
    const uint32_t sketchMask = __ballot_sync(0xFFFFFFFF, validFeature);
    const uint32_t sketchSize = __popc(sketchMask);

    // only process non-empty queries
    if (sketchSize) {
        if (warpLane < sketchSize)
            sketch[warpLane] = f;

        __syncwarp();
    }

    return sketchSize;
}


/**************************************************************************//**
 *
 * @brief for each feature in sketch insert location into hash table
 *
 * @param hashtable  hashtable where sketch is inserted
 * @param sketch     pointer to sketch array (keys to insert)
 * @param sketchSize sketch size
 * @param loc        location (value to insert)
 *
 *****************************************************************************/
template<
    class Hashtable,
    class Feature>
__device__ __inline__
void insert_into_hashtable(
    Hashtable& hashtable,
    const Feature * sketch,
    uint32_t sketchSize,
    typename Hashtable::value_type loc)
{
    namespace cg = cooperative_groups;

    const auto group =
        cg::tiled_partition<Hashtable::cg_size()>(cg::this_thread_block());

    const int warpLane = threadIdx.x % WARPSIZE;
    const int groupId = warpLane / Hashtable::cg_size();
    constexpr int groupsPerWarp = WARPSIZE / Hashtable::cg_size();

    for (int i = groupId; i < sketchSize; i += groupsPerWarp)
    {
        const auto status = hashtable.insert(
            sketch[i], loc, group, 1000000);

        // if (group.thread_rank() == 0 && status.has_any()) {
        //     printf("status %d\n", status.has_any());
        // }
    }
}


/**************************************************************************//**
 *
 * @brief query the hash table with each feature in sketch;
 *        write the result back into sketch in place of the feature
 *
 * @param hashtable  hashtable which is queried
 * @param sketch     pointer to sketch array (keys to query)
 * @param sketchSize sketch size
 *
 *****************************************************************************/
template<
    class Hashtable,
    class Feature>
__device__ __inline__
void query_hashtable(
    Hashtable& hashtable,
    Feature * sketch,
    uint32_t sketchSize)
{
    namespace cg = cooperative_groups;

    const auto group =
        cg::tiled_partition<Hashtable::cg_size()>(cg::this_thread_block());

    const int warpLane = threadIdx.x % WARPSIZE;
    const int groupId = warpLane / Hashtable::cg_size();
    constexpr int groupsPerWarp = WARPSIZE / Hashtable::cg_size();

    for (int i = groupId; i < sketchSize; i += groupsPerWarp)
    {
        typename Hashtable::value_type valuesOffset = 0;
        // if key not found valuesOffset stays 0
        const auto status = hashtable.retrieve(
            sketch[i], valuesOffset, group);

        if (group.thread_rank() == 0) {
            // if (status.has_any()) {
                // printf("status %d\n", status.has_any());
                // printf("status %d\n", status.has_key_not_found());
            // }

            sketch[i] = valuesOffset;
        }
    }
}


/**************************************************************************//**
 *
 * @brief for each feature in sketch copy its loctions to the output
 *
 * @param sketch     pointer to sketch array (offset and size of location buckets)
 * @param sketchSize sketch size
 * @param locations  pointer location bucket storage
 * @param maxLocationsPerFeature
 *                   maximum number of locations per feature to copy
 * @param out        pointer to output array
 * @param outOffset  pointer to output offset
 *
 * @return total number of locations copied for sketch
 *
 *****************************************************************************/
template<
    int MAX_SKETCH_SIZE,
    class Feature,
    class Location,
    class BucketSizeT>
__device__ __inline__
void copy_loctions(
    Feature * sketch,
    uint32_t sketchSize,
    const Location * locations,
    BucketSizeT maxLocationsPerFeature,
    Location * out,
    int * outOffset)
{
    using location = Location;
    using bucket_size_type = BucketSizeT;

    const int warpLane = threadIdx.x % WARPSIZE;

    uint32_t totalLocations = 0;

    for (uint32_t s = warpLane; s < sketchSize; s += WARPSIZE) {
        auto bucketOffset = sketch[s];
        bucket_size_type bucketSize = min(bucket_size_type(bucketOffset), maxLocationsPerFeature);
        totalLocations += bucketSize;
    }

    // reduction into first thread
    for (int i = 1; i < MAX_SKETCH_SIZE; i *= 2) {
        totalLocations += __shfl_down_sync(0xFFFFFFFF, totalLocations, i);
    }

    int offset = 0;
    if (warpLane == 0) {
        // reserve locations for this query
        offset = atomicAdd(outOffset, totalLocations);
    }
    offset = __shfl_sync(0xFFFFFFFF, offset, 0);
    out += offset;

    for (uint32_t s = 0; s < sketchSize; ++s) {
        auto bucketOffset = sketch[s];
        bucket_size_type bucketSize = min(bucket_size_type(bucketOffset), maxLocationsPerFeature);
        bucketOffset >>= sizeof(bucket_size_type)*CHAR_BIT;

        // copy locations
        for (uint32_t i = warpLane; i < bucketSize; i += WARPSIZE) {
            out[i] = locations[bucketOffset + i];
        }

        out += bucketSize;
    }
}


/**************************************************************************//**
 *
 * @brief query the hash table with each feature in sketch;
 *        for each feature copy its loctions to the output
 *
 * @param hashtable  hashtable which is queried
 * @param sketch     pointer to sketch array (keys to query)
 * @param sizes      pointer to sizes array (sizes of location buckets)
 * @param sketchSize sketch size
 * @param out        pointer to output array
 * @param outOffset  pointer to output offset
 *
 *****************************************************************************/
template<
    int MAX_SKETCH_SIZE,
    class Hashtable>
__device__ __inline__
void query_bucket_hashtable(
    Hashtable& hashtable,
    const feature * sketch,
    uint32_t * sizes,
    uint32_t sketchSize,
    location * out,
    int * outOffset)
{
    using slot = typename Hashtable::value_type;

    namespace cg = cooperative_groups;

    const auto group =
        cg::tiled_partition<Hashtable::cg_size()>(cg::this_thread_block());

    const int warpLane = threadIdx.x % WARPSIZE;
    const int groupId = warpLane / Hashtable::cg_size();
    constexpr int groupsPerWarp = WARPSIZE / Hashtable::cg_size();

    for (int i = groupId; i < sketchSize; i += groupsPerWarp)
    {
        uint64_t numValues = 0;
        // if key not found numValues stays 0
        const auto status = hashtable.num_values(
            sketch[i], numValues, group);

        if (group.thread_rank() == 0) {
            // if (status.has_any()) {
                // printf("status %d\n", status.has_any());
                // printf("status %d\n", status.has_key_not_found());
            // }

            sizes[i] = numValues;
        }
    }

    __syncwarp();

    // prefix sum to get output offsets
    uint32_t numValues = (warpLane < sketchSize) ? sizes[warpLane] : 0;
    for (int i = 1; i < MAX_SKETCH_SIZE; i *= 2) {
        const uint32_t other = __shfl_up_sync(0xFFFFFFFF, numValues, i);
        if (warpLane >= i)
            numValues += other;
    }
    if (warpLane < MAX_SKETCH_SIZE)
        sizes[warpLane] = numValues;

    int globalOffset = 0;
    if (warpLane == MAX_SKETCH_SIZE-1) {
        // reserve locations for this query
        globalOffset = atomicAdd(outOffset, numValues);
    }
    globalOffset = __shfl_sync(0xFFFFFFFF, globalOffset, MAX_SKETCH_SIZE-1);
    out += globalOffset;

    for (int i = groupId; i < sketchSize; i += groupsPerWarp)
    {
        const uint32_t offset = (i != 0) ? sizes[i-1] : 0;
        const uint32_t size  = sizes[i] - offset;

        // only retrieve non empty buckets
        if (size) {
            uint64_t dummy;
            const auto status = hashtable.retrieve(
                sketch[i], reinterpret_cast<slot*>(out + offset), dummy, group);
        }
    }
}




/**************************************************************************//**
 *
 * @brief insert batch of sequences into bucket list hashtable
 *
 * @details each warp processes one sequence (window)
 *
 *****************************************************************************/
template<
    int BLOCK_THREADS,
    int MAX_SKETCH_SIZE,
    class Hashtable,
    class IndexT,
    class Offset>
__global__
__launch_bounds__(BLOCK_THREADS)
void insert_features(
    Hashtable hashtable,
    IndexT numTargets,
    const target_id * targetIds,
    const window_id * winOffsets,
    const char      * sequence,
    const Offset    * sequenceOffsets,
    const sketching_opt targetSketching
)
{
    using index_type = IndexT;
    using offset_type = Offset;

    const auto& kmerSize = targetSketching.kmerlen;
    const auto& maxSketchSize = targetSketching.sketchlen;
    const auto& windowSize = targetSketching.winlen;
    const auto& windowStride = targetSketching.winstride;

    constexpr int warpsPerBlock = BLOCK_THREADS / WARPSIZE;

    __shared__ feature sketch[warpsPerBlock][MAX_SKETCH_SIZE];

    for (index_type tgt = blockIdx.y; tgt < numTargets; tgt += gridDim.y) {
        const offset_type sequenceBegin = sequenceOffsets[tgt];
        const offset_type sequenceEnd   = sequenceOffsets[tgt+1];
        const offset_type sequenceLength = sequenceEnd - sequenceBegin;
        const window_id lastWindowSize = sequenceLength % windowStride;
        const window_id numWindows = sequenceLength / windowStride + (lastWindowSize >= kmerSize);

        const int warpId = threadIdx.x / WARPSIZE;

        // each warp processes one window
        for (window_id win = blockIdx.x * warpsPerBlock + warpId;
            win < numWindows;
            win += warpsPerBlock * gridDim.x)
        {
            const offset_type offsetBegin = sequenceBegin + win*windowStride;
            const offset_type offsetEnd = min(sequenceEnd, offsetBegin + windowSize);
            const char * windowBegin = sequence + offsetBegin;
            const char * windowEnd = sequence + offsetEnd;

            uint32_t sketchSize = warp_make_sketch<offset_type>(
                windowBegin, windowEnd, kmerSize, sketch[warpId], maxSketchSize);

            const target_id targetId = targetIds[tgt];
            const window_id winOffset = winOffsets[tgt]+win ;
            const uint64_t loc = (uint64_t(targetId) << 32) + winOffset;

            // const int warpLane = threadIdx.x % WARPSIZE;
            // if (warpLane == 0)
            //     printf("tgt %d win %d features %u %u %u %u %u %u %u %u "
            //                                   "%u %u %u %u %u %u %u %u\n",
            //         tgt, win,
            //         sketch[warpId][0],
            //         sketch[warpId][1],
            //         sketch[warpId][2],
            //         sketch[warpId][3],
            //         sketch[warpId][4],
            //         sketch[warpId][5],
            //         sketch[warpId][6],
            //         sketch[warpId][7],
            //         sketch[warpId][8],
            //         sketch[warpId][9],
            //         sketch[warpId][10],
            //         sketch[warpId][11],
            //         sketch[warpId][12],
            //         sketch[warpId][13],
            //         sketch[warpId][14],
            //         sketch[warpId][15]);

            insert_into_hashtable(
                hashtable, sketch[warpId], sketchSize, loc);

            __syncwarp();
        }
    }
}


/**************************************************************************//**
 *
 * @brief query batch of sequences against hashtable
 *
 * @details each warp processes one sequence (window)
 *
 *****************************************************************************/
template<
    int BLOCK_THREADS,
    int MAX_SKETCH_SIZE,
    class Hashtable,
    class Offset,
    class BucketSizeT,
    class Location>
__global__
__launch_bounds__(BLOCK_THREADS)
void gpu_hahstable_query(
    Hashtable hashtable,
    const uint32_t * readIds,
    uint32_t numQueries,
    const Offset * sequenceOffsets,
    const char * sequences,
    feature * sketches,
    numk_t kmerSize,
    uint32_t maxSketchSize,
    const Location * locations,
    BucketSizeT maxLocationsPerFeature,
    Location * locationsOut,
    int      * outOffset
)
{
    using location = Location;

    constexpr int warpsPerBlock = BLOCK_THREADS / WARPSIZE;

    // using uint64_t enables to store either "feature" or "Hashtable::value_type"
    __shared__ uint64_t sketch[warpsPerBlock][MAX_SKETCH_SIZE];

    const int warpId = threadIdx.x / WARPSIZE;
    const int warpLane = threadIdx.x % WARPSIZE;

    // each warp processes one query (= one window)
    for (int queryId = blockIdx.x * warpsPerBlock + warpId;
        queryId < numQueries;
        queryId += warpsPerBlock * gridDim.x)
    {
        uint32_t sketchSize = 0;

        if (sequences != nullptr)
        {
            const char * sequenceBegin = sequences + sequenceOffsets[queryId];
            const char * sequenceEnd   = sequences + sequenceOffsets[queryId+1];

            sketchSize = warp_make_sketch<Offset>(
                sequenceBegin, sequenceEnd, kmerSize, sketch[warpId], maxSketchSize);

            if (sketches != nullptr) {
                // store sketch for multi-gpu processing
                if (warpLane < maxSketchSize)
                    sketches[queryId*maxSketchSize + warpLane] =
                        (warpLane < sketchSize) ? sketch[warpId][warpLane] : feature(~0);
            }
        }
        else
        {
            sketchSize = warp_load_sketch(sketches+queryId*maxSketchSize, sketch[warpId], maxSketchSize);
        }

        // only process non-empty queries
        if (sketchSize) {
            // if (warpLane == 0)
            //     printf("query %d features %llu %llu %llu %llu %llu %llu %llu %llu "
            //                              "%llu %llu %llu %llu %llu %llu %llu %llu\n",
            //         queryId,
            //         sketch[warpId][0],
            //         sketch[warpId][1],
            //         sketch[warpId][2],
            //         sketch[warpId][3],
            //         sketch[warpId][4],
            //         sketch[warpId][5],
            //         sketch[warpId][6],
            //         sketch[warpId][7],
            //         sketch[warpId][8],
            //         sketch[warpId][9],
            //         sketch[warpId][10],
            //         sketch[warpId][11],
            //         sketch[warpId][12],
            //         sketch[warpId][13],
            //         sketch[warpId][14],
            //         sketch[warpId][15]);

            query_hashtable(hashtable, sketch[warpId], sketchSize);

            __syncwarp();

            copy_loctions<MAX_SKETCH_SIZE>(
                sketch[warpId], sketchSize, locations, maxLocationsPerFeature,
                locationsOut, outOffset + readIds[queryId]);

            __syncwarp();
        }
    }
}


/**************************************************************************//**
 *
 * @brief query batch of sequences against bucket list hashtable
 *
 * @details each warp processes one sequence (window)
 *
 *****************************************************************************/
template<
    int BLOCK_THREADS,
    int MAX_SKETCH_SIZE,
    class Hashtable,
    class Offset,
    class BucketSizeT,
    class Location>
__global__
__launch_bounds__(BLOCK_THREADS)
void gpu_hahstable_query(
    Hashtable hashtable,
    const uint32_t * readIds,
    uint32_t numQueries,
    const Offset * sequenceOffsets,
    const char * sequences,
    feature * sketches,
    numk_t kmerSize,
    uint32_t maxSketchSize,
    BucketSizeT maxLocationsPerFeature,
    Location * locationsOut,
    int      * outOffset
)
{
    using location = Location;

    constexpr int warpsPerBlock = BLOCK_THREADS / WARPSIZE;

    __shared__ feature sketch[warpsPerBlock][MAX_SKETCH_SIZE];
    __shared__ uint32_t sizes[warpsPerBlock][MAX_SKETCH_SIZE];

    const int warpId = threadIdx.x / WARPSIZE;
    const int warpLane = threadIdx.x % WARPSIZE;

    // each warp processes one query (= one window)
    for (int queryId = blockIdx.x * warpsPerBlock + warpId;
        queryId < numQueries;
        queryId += warpsPerBlock * gridDim.x)
    {
        uint32_t sketchSize = 0;

        if (sequences != nullptr)
        {
            const char * sequenceBegin = sequences + sequenceOffsets[queryId];
            const char * sequenceEnd   = sequences + sequenceOffsets[queryId+1];

            sketchSize = warp_make_sketch<Offset>(
                sequenceBegin, sequenceEnd, kmerSize, sketch[warpId], maxSketchSize);

            if (sketches != nullptr) {
                // store sketch for multi-gpu processing
                if (warpLane < maxSketchSize)
                    sketches[queryId*maxSketchSize + warpLane] =
                        (warpLane < sketchSize) ? sketch[warpId][warpLane] : feature(~0);
            }
        }
        else
        {
            sketchSize = warp_load_sketch(sketches+queryId*maxSketchSize, sketch[warpId], maxSketchSize);
        }

        // only process non-empty queries
        if (sketchSize) {
            // if (warpLane == 0)
            //     printf("query %d features %llu %llu %llu %llu %llu %llu %llu %llu "
            //                              "%llu %llu %llu %llu %llu %llu %llu %llu\n",
            //         queryId,
            //         sketch[warpId][0],
            //         sketch[warpId][1],
            //         sketch[warpId][2],
            //         sketch[warpId][3],
            //         sketch[warpId][4],
            //         sketch[warpId][5],
            //         sketch[warpId][6],
            //         sketch[warpId][7],
            //         sketch[warpId][8],
            //         sketch[warpId][9],
            //         sketch[warpId][10],
            //         sketch[warpId][11],
            //         sketch[warpId][12],
            //         sketch[warpId][13],
            //         sketch[warpId][14],
            //         sketch[warpId][15]);

            query_bucket_hashtable<MAX_SKETCH_SIZE>(
                hashtable, sketch[warpId], sizes[warpId], sketchSize,
                locationsOut, outOffset + readIds[queryId]);

            __syncwarp();
        }
    }
}


} // namespace mc

#endif

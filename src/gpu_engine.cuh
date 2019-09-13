#ifndef MC_GPU_ENGINE_H_
#define MC_GPU_ENGINE_H_

#include "config.h"
#include "dna_encoding.h"
#include "sketch_database.h"
#include "../dep/cudahelpers/cuda_helpers.cuh"
#include "../dep/cub/cub/block/block_radix_sort.cuh"
#include "../dep/cub/cub/block/block_store.cuh"

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
    target_id * targetIds,
    window_id * winOffsets,
    uint32_t  * encodedSeq,
    encodedseq_t* seq,
    encodedambig_t* ambig,
    numk_t k,
    size_t windowStride,
    uint32_t sketchSize,
    kmer_type* features_out,
    location* locations_out,
    uint64_t* size_out)
{
    const target_id tgt = targetIds[blockIdx.y];
    const window_id winOffset = winOffsets[blockIdx.y];
    const uint32_t encodedLength = encodedSeq[blockIdx.y+1] - encodedSeq[blockIdx.y];

    typedef cub::BlockRadixSort<kmer_type, BLOCK_THREADS, ITEMS_PER_THREAD> BlockRadixSortT;

    __shared__ typename BlockRadixSortT::TempStorage temp_storage;

    constexpr uint8_t encBits   = sizeof(encodedseq_t)*CHAR_BIT;
    constexpr uint8_t ambigBits = sizeof(encodedambig_t)*CHAR_BIT;

    //letters stored from high to low bits
    //masks get highest bits
    const encodedseq_t   kmerMask  = encodedseq_t(~0) << (encBits - 2*k);
    const encodedambig_t ambigMask = encodedambig_t(~0) << (ambigBits - k);

    const size_t seqLength = encodedLength * ambigBits;
    const window_id numWindows = (seqLength-k + windowStride) / windowStride;

    //each block processes one window
    for(window_id bid = blockIdx.x; bid < numWindows; bid += gridDim.x)
    {
        kmer_type items[ITEMS_PER_THREAD];
        uint8_t numItems = 0;
        for(int i=0; i<ITEMS_PER_THREAD; ++i)
        {
            items[i] = kmer_type(~0);
        }

        //each thread extracts one feature
        const size_t offset = bid * windowStride;
        for(size_t tid  = offset + threadIdx.x;
                   tid  < min(offset + windowStride, seqLength-k+1);
                   tid += blockDim.x)
        {
            const std::uint32_t  seq_slot    = tid / (ambigBits);
            const std::uint8_t   kmer_slot   = tid & (ambigBits-1);
            const encodedambig_t ambig_left  = ambig[seq_slot] << kmer_slot;
            const encodedambig_t ambig_right = ambig[seq_slot+1] >> (ambigBits-kmer_slot);

            //continue only if no bases ambiguous
            if(!((ambig_left+ambig_right) & ambigMask))
            {
                const encodedseq_t   seq_left    = seq[seq_slot] << (2*kmer_slot);
                const encodedseq_t   seq_right   = seq[seq_slot+1] >> (encBits-(2*kmer_slot));

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
            const kmer_type kmer = items[i];
            kmer_type nextKmer = threadIdx.x ? items[i] : items[i+1];
            nextKmer = __shfl_sync(0xFFFFFFFF, nextKmer, threadIdx.x+1);

            //find valid candidates
            const bool predicate = (kmer != kmer_type(~0)) && (kmer != nextKmer);

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
                features_out[pos]  = kmer;
                locations_out[pos] = location{winOffset+bid, tgt};
            }

            kmerCounter += capped;
        }
    }
}



} // namespace mc


#endif
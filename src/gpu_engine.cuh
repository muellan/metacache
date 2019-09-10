#ifndef MC_GPU_ENGINE_H_
#define MC_GPU_ENGINE_H_

#include "config.h"
#include "dna_encoding.h"
#include "../dep/cudahelpers/cuda_helpers.cuh"

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


template<bool CanonicalKmers = true>
__global__
void extract_kmers(
    encodedseq_t* seq_in,
    encodedambig_t* ambig_in,
    size_t size_in,
    numk_t k,
    kmer_type* kmers_out,
    uint64_t* size_out)
{
    //bases stored from high to low bits
    //masks get highest bits
    const encodedseq_t   kmer_mask  = encodedseq_t(~0) << (sizeof(encodedseq_t)*CHAR_BIT - 2*k);
    const encodedambig_t ambig_mask = encodedambig_t(~0) << (sizeof(encodedambig_t)*CHAR_BIT - k);

    for (size_t tid  = blockIdx.x * blockDim.x + threadIdx.x;
                tid  < ((size_in*sizeof(encodedambig_t)*CHAR_BIT)-k+1);
                tid += blockDim.x * gridDim.x)
    {
        const std::uint32_t  seq_slot    = tid / sizeof(encodedambig_t)*CHAR_BIT;
        const std::uint8_t   kmer_slot   = tid & (sizeof(encodedambig_t)*CHAR_BIT-1);
        const encodedseq_t   seq_left    = seq_in[seq_slot] << (2*kmer_slot);
        const encodedseq_t   seq_right   = seq_in[seq_slot+1] >> (sizeof(encodedseq_t)*CHAR_BIT-(2*kmer_slot));
        const encodedambig_t ambig_left  = ambig_in[seq_slot] << kmer_slot;
        const encodedambig_t ambig_right = ambig_in[seq_slot+1] >> (sizeof(encodedambig_t)*CHAR_BIT-kmer_slot);

        //continue only if no bases ambiguous
        if(!((ambig_left+ambig_right) & ambig_mask))
        {
            //get highest bits
            kmer_type kmer = (seq_left+seq_right) & kmer_mask;
            //shift kmer to lowest bits
            kmer >>= (sizeof(encodedseq_t)*CHAR_BIT - 2*k);

            if(CanonicalKmers)
            {
                kmer = make_canonical_2bit(kmer);
            }

            kmers_out[atomicAggInc(size_out)] = kmer;
        }
    }
}



} // namespace mc


#endif
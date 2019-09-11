#include <limits>

#include "gpu_hashmap.cuh"
#include "hash_dna.h"
#include "hash_int.h"
#include "sketch_database.h"
#include "gpu_engine.cuh"

namespace mc {

    __global__
    void insert_kernel(target_id tgt, encodedseq_t* encodedSeq, encodedambig_t* encodedAmbig) {

    }


    //-----------------------------------------------------
    template<
        class Key,
        class ValueT,
        class Hash,
        class KeyEqual,
        class BucketSizeT
    >
    void gpu_hashmap<Key,ValueT,Hash,KeyEqual,BucketSizeT>::init()
    {
        cudaSetDevice(0);
    }


    //-----------------------------------------------------
    template<
        class Key,
        class ValueT,
        class Hash,
        class KeyEqual,
        class BucketSizeT
    >
    void gpu_hashmap<Key,ValueT,Hash,KeyEqual,BucketSizeT>::insert(
        target_id tgt,
        std::vector<encodedseq_t> encodedSeq,
        std::vector<encodedambig_t> encodedAmbig,
        size_t seqLength,
        numk_t k, size_t windowStride)
    {
        encodedseq_t * d_encodedSeq;
        encodedambig_t * d_encodedAmbig;
        size_t encodedLength = encodedSeq.size();

        cudaMalloc(&d_encodedSeq, encodedLength*sizeof(encodedseq_t));
        cudaMalloc(&d_encodedAmbig, encodedLength*sizeof(encodedambig_t));
        CUERR

        // kmer_type * h_kmers;
        kmer_type * d_kmers;
        size_t numKmers = seqLength-k+1;
        uint64_t * d_kmerCounter;

        std::vector<kmer_type> h_kmers(numKmers);

        // cudaMallocHost(&h_kmers, numKmers*sizeof(kmer_type));
        cudaMalloc(&d_kmers, numKmers*sizeof(kmer_type));
        CUERR
        cudaMalloc(&d_kmerCounter, sizeof(uint64_t));
        cudaMemset(d_kmerCounter, 0, sizeof(uint64_t));
        CUERR

        cudaMemcpy(d_encodedSeq, encodedSeq.data(),
                   encodedLength*sizeof(encodedseq_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_encodedAmbig, encodedAmbig.data(),
                   encodedLength*sizeof(encodedambig_t), cudaMemcpyHostToDevice);
        CUERR

        // insert_kernel<<<1,1>>>(tgt, d_encodedSeq, d_encodedAmbig);


        extract_features<32,4><<<1,32>>>(
            d_encodedSeq, d_encodedAmbig, seqLength,
            k, windowStride,
            d_kmers, d_kmerCounter);
        cudaDeviceSynchronize();
        CUERR

        cudaMemcpy(h_kmers.data(), d_kmers, numKmers*sizeof(kmer_type), cudaMemcpyDeviceToHost);
        CUERR

        std::cout << "Target ID: " << tgt << '\n';
        //sort kmers
        // std::sort(h_kmers.begin(), h_kmers.end());
        //print kmers
        size_t kmerCounter = 0;
        for(const auto& kmer : h_kmers) {
            std:: cout << kmer << ' ';
            if(++kmerCounter % windowStride == 0)
                std::cout << '\n';
        }
        std::cout << std::endl;
    }





    //-----------------------------------------------------
    template class gpu_hashmap<
            unsigned int,
            // sketch_database<
            //     std::__cxx11::basic_string<
            //         char,
            //         std::char_traits<char>,
            //         std::allocator<char> >,
            //     single_function_unique_min_hasher<
            //         unsigned int, same_size_hash<unsigned int> >,
            //     same_size_hash<unsigned int>,
            //     unsigned short,
            //     unsigned int,
            //     unsigned char
            //     >::target_location,
            uint64_t,
            same_size_hash<unsigned int>,
            std::equal_to<unsigned int>,
            unsigned char
            >;

} // namespace mc

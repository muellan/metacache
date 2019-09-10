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
        numk_t k)
    {
        encodedseq_t * d_encodedSeq;
        encodedambig_t * d_encodedAmbig;
        size_t encodedLength = encodedSeq.size();

        cudaMalloc(&d_encodedSeq, encodedLength*sizeof(encodedseq_t));
        cudaMalloc(&d_encodedAmbig, encodedLength*sizeof(encodedambig_t));

        kmer_type * h_kmers, * d_kmers;
        size_t numKmers = encodedLength*sizeof(encodedambig_t)*CHAR_BIT-k+1;
        uint64_t * d_kmerCounter;

        cudaMallocHost(&h_kmers, numKmers*sizeof(kmer_type));
        cudaMalloc(&d_kmers, numKmers*sizeof(kmer_type));
        cudaMalloc(&d_kmerCounter, sizeof(uint64_t));

        cudaMemcpy(d_encodedSeq, encodedSeq.data(),
                   encodedLength*sizeof(encodedseq_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_encodedAmbig, encodedAmbig.data(),
                   encodedLength*sizeof(encodedambig_t), cudaMemcpyHostToDevice);

        // insert_kernel<<<1,1>>>(tgt, d_encodedSeq, d_encodedAmbig);

        extract_kmers<<<1024,1024>>>(d_encodedSeq, d_encodedAmbig, encodedLength,
                                     k, d_kmers, d_kmerCounter);

        cudaMemcpy(h_kmers, d_kmers, numKmers*sizeof(kmer_type), cudaMemcpyDeviceToHost);

        std::cout << "Target ID: " << tgt << '\n';

        //print kmers
        for(size_t i=0; i<numKmers; ++i) {
            std:: cout << h_kmers[i] << ' ';
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

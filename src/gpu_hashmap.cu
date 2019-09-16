#include <limits>

#include "gpu_hashmap.cuh"
#include "hash_dna.h"
#include "hash_int.h"
#include "sketch_database.h"
#include "gpu_engine.cuh"

namespace mc {

    //-----------------------------------------------------
    template<
        class Key,
        class ValueT,
        class Hash,
        class KeyEqual,
        class BucketSizeT
    >
    void gpu_hashmap<Key,ValueT,Hash,KeyEqual,BucketSizeT>::init(
        numk_t kmerLength,
        sketch_size_type sketchSize,
        size_t windowStride,
        size_t windowSize
    )
    {
        kmerLength_ = kmerLength;
        sketchSize_ = sketchSize;
        windowStride_ = windowStride;
        windowSize_ = windowSize;

        cudaSetDevice(0);
        size_t maxTargets = 1;
        size_t maxEncodeLength = 100;
        size_t maxSeqLength = maxEncodeLength*sizeof(encodedambig_t)*CHAR_BIT;
        size_t maxWindows = (maxSeqLength-kmerLength_ + windowStride_) / windowStride_;
        size_t maxFeatures = maxWindows * sketchSize_;

        cudaMalloc(&seqBatch_.targetIds, maxTargets*sizeof(target_id));
        cudaMalloc(&seqBatch_.windowOffsets, maxTargets*sizeof(window_id));
        cudaMalloc(&seqBatch_.encodeOffsets, (maxTargets+1)*sizeof(uint32_t));

        cudaMalloc(&seqBatch_.encodedSeq, maxEncodeLength*sizeof(encodedseq_t));
        cudaMalloc(&seqBatch_.encodedAmbig, maxEncodeLength*sizeof(encodedambig_t));

        cudaMalloc(&seqBatch_.features, maxFeatures*sizeof(Key));
        cudaMalloc(&seqBatch_.values, maxFeatures*sizeof(ValueT));
        cudaMalloc(&seqBatch_.featureCounter, sizeof(size_t));
        CUERR
    }


    //-----------------------------------------------------
    template<
        class Key,
        class ValueT,
        class Hash,
        class KeyEqual,
        class BucketSizeT
    >
    std::vector<Key> gpu_hashmap<Key,ValueT,Hash,KeyEqual,BucketSizeT>::insert(
        target_id tgt,
        std::vector<encodedseq_t> encodedSeq,
        std::vector<encodedambig_t> encodedAmbig)
    {

        const uint32_t encodedLength = encodedSeq.size();
        //not real seq length but upper bound
        const size_t seqLength = encodedLength*sizeof(encodedambig_t)*CHAR_BIT;
        const size_t numWindows = (seqLength-kmerLength_ + windowStride_) / windowStride_;
        const size_t numFeatures = numWindows * sketchSize_;
        std::cout << "Target ID: " << tgt
                  << " Length: " << seqLength
                  << " Windows: " << numWindows
                  << " Features: " << numFeatures
                  << '\n';

        cudaStream_t stream = 0;

        cudaMemsetAsync(seqBatch_.featureCounter, 0, sizeof(uint64_t), stream);

        cudaMemcpyAsync(seqBatch_.targetIds, &tgt,
                   sizeof(target_id), cudaMemcpyHostToDevice, stream);
        const window_id winOffset = 0;
        cudaMemcpyAsync(seqBatch_.windowOffsets, &winOffset,
                   sizeof(window_id), cudaMemcpyHostToDevice, stream);
        const uint32_t encOffsets[2] = {0,encodedLength};
        cudaMemcpyAsync(seqBatch_.encodeOffsets, encOffsets,
                   2*sizeof(uint32_t), cudaMemcpyHostToDevice, stream);
        cudaMemcpyAsync(seqBatch_.encodedSeq, encodedSeq.data(),
                   encodedLength*sizeof(encodedseq_t), cudaMemcpyHostToDevice, stream);
        cudaMemcpyAsync(seqBatch_.encodedAmbig, encodedAmbig.data(),
                   encodedLength*sizeof(encodedambig_t), cudaMemcpyHostToDevice, stream);
        // cudaStreamSynchronize(stream);
        // CUERR

        #define BLOCK_THREADS 32

        extract_features<BLOCK_THREADS,4><<<1,BLOCK_THREADS,0,stream>>>(
            seqBatch_.targetIds,
            seqBatch_.windowOffsets,
            seqBatch_.encodeOffsets,
            seqBatch_.encodedSeq,
            seqBatch_.encodedAmbig,
            kmerLength_, sketchSize_, windowStride_, windowSize_,
            seqBatch_.features,
            seqBatch_.values,
            seqBatch_.featureCounter);
        // cudaStreamSynchronize(stream);
        // CUERR

        uint64_t h_featureCounter = 0;
        cudaMemcpyAsync(&h_featureCounter, seqBatch_.featureCounter,
                   sizeof(uint64_t), cudaMemcpyDeviceToHost, stream);
        // cudaStreamSynchronize(stream);
        // std::cout << "Counter: " << h_featureCounter << '\n';

        // kmer_type * h_features;
        // cudaMallocHost(&h_features, numFeatures*sizeof(Key));
        std::vector<Key> h_features(h_featureCounter);
        std::vector<ValueT> h_values(h_featureCounter);
        cudaMemcpyAsync(h_features.data(), seqBatch_.features,
                   h_featureCounter*sizeof(Key), cudaMemcpyDeviceToHost, stream);
        cudaMemcpyAsync(h_values.data(), seqBatch_.values,
                   h_featureCounter*sizeof(ValueT), cudaMemcpyDeviceToHost, stream);
        // cudaStreamSynchronize(stream);
        // CUERR

        //print features
        // for(size_t i=0; i<h_featureCounter; ++i) {
        //     std:: cout << h_features[i] << ' ';
        //     if((i+1) % sketchSize_ == 0)
        //         std::cout << '\n';
        // }
        // std::cout << std::endl;

        // for(size_t i=0; i<h_featureCounter; ++i) {
        //     std:: cout << '(' << h_values[i].tgt << ',' << h_values[i].win << ") " ;
        //     if((i+1) % sketchSize_ == 0)
        //         std::cout << '\n';
        // }
        std::cout << std::endl;

        cudaStreamSynchronize(stream);
        CUERR

        return h_features;
    }





    //-----------------------------------------------------
    template class gpu_hashmap<
            kmer_type,
            location,
            // uint64_t,
            same_size_hash<kmer_type>,
            std::equal_to<kmer_type>,
            unsigned char
            >;

} // namespace mc

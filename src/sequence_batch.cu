
#include "sequence_batch.cuh"

namespace mc {


//---------------------------------------------------------------
template<>
sequence_batch<policy::Host>::sequence_batch(uint32_t maxTargets, size_t maxEncodeLength) :
    maxTargets_{maxTargets}, numTargets_{0}, maxEncodeLength_{maxEncodeLength}
{
    if(maxTargets_) {
        cudaMallocHost(&targetIds_, maxTargets_*sizeof(target_id));
        cudaMallocHost(&windowOffsets_, maxTargets_*sizeof(window_id));
        cudaMallocHost(&encodeOffsets_, (maxTargets_+1)*sizeof(encodinglen_t));
        encodeOffsets_[0] = 0;
    }
    if(maxEncodeLength_) {
        cudaMallocHost(&encodedSeq_, maxEncodeLength_*sizeof(encodedseq_t));
        cudaMallocHost(&encodedAmbig_, maxEncodeLength_*sizeof(encodedambig_t));
    }
    CUERR
}
//---------------------------------------------------------------
template<>
sequence_batch<policy::Host>::~sequence_batch() {
    if(maxTargets_) {
        cudaFreeHost(targetIds_);
        cudaFreeHost(windowOffsets_);
        cudaFreeHost(encodeOffsets_);
    }
    if(maxEncodeLength_) {
        cudaFreeHost(encodedSeq_);
        cudaFreeHost(encodedAmbig_);
    }
    CUERR
}

//---------------------------------------------------------------
template<>
sequence_batch<policy::Device>::sequence_batch(uint32_t maxTargets, size_t maxEncodeLength) :
    maxTargets_{maxTargets}, numTargets_{0}, maxEncodeLength_{maxEncodeLength}
{
    if(maxTargets_) {
        cudaMalloc(&targetIds_, maxTargets_*sizeof(target_id));
        cudaMalloc(&windowOffsets_, maxTargets_*sizeof(window_id));
        cudaMalloc(&encodeOffsets_, (maxTargets_+1)*sizeof(encodinglen_t));
    }
    if(maxEncodeLength_) {
        cudaMalloc(&encodedSeq_, maxEncodeLength_*sizeof(encodedseq_t));
        cudaMalloc(&encodedAmbig_, maxEncodeLength_*sizeof(encodedambig_t));
    }
    CUERR

    size_t totalSize = maxTargets_*(sizeof(target_id) + sizeof(window_id)) +
                       (maxTargets_+1)*sizeof(encodinglen_t) +
                       maxEncodeLength_*(sizeof(encodedseq_t) + sizeof(encodedambig_t));
    std::cerr << "total batch size: " << (totalSize >> 10) << " KB\n";
}
//---------------------------------------------------------------
template<>
sequence_batch<policy::Device>::~sequence_batch() {
    if(maxTargets_) {
        cudaFree(targetIds_);
        cudaFree(windowOffsets_);
        cudaFree(encodeOffsets_);
    }
    if(maxEncodeLength_) {
        cudaFree(encodedSeq_);
        cudaFree(encodedAmbig_);
    }
    CUERR
}


} // namespace mc

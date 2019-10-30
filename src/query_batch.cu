
#include "query_batch.cuh"
#include "gpu_engine.cuh"

namespace mc {


//---------------------------------------------------------------
template<class result_type>
query_batch<result_type>::query_batch(size_t maxQueries, size_t maxEncodeLength, size_t maxResultsPerQuery) :
    numQueries_{0},
    maxQueries_{maxQueries},
    maxEncodeLength_{maxEncodeLength},
    maxResultsPerQuery_{maxResultsPerQuery}
{
    if(maxQueries_) {
        cudaMallocHost(&queryIds_, maxQueries_*sizeof(query_id));

        cudaMallocHost(&h_encodeOffsets_, (maxQueries_+1)*sizeof(encodinglen_t));
        h_encodeOffsets_[0] = 0;
        cudaMallocHost(&h_queryResults_, maxQueries_*maxResultsPerQuery_*sizeof(result_type));

        cudaMalloc(&d_encodeOffsets_, (maxQueries_+1)*sizeof(encodinglen_t));
        cudaMalloc(&d_queryResults_, maxQueries_*maxResultsPerQuery_*sizeof(result_type));
    }
    if(maxEncodeLength_) {
        cudaMallocHost(&h_encodedSeq_, maxEncodeLength_*sizeof(encodedseq_t));
        cudaMallocHost(&h_encodedAmbig_, maxEncodeLength_*sizeof(encodedambig_t));

        cudaMalloc(&d_encodedSeq_, maxEncodeLength_*sizeof(encodedseq_t));
        cudaMalloc(&d_encodedAmbig_, maxEncodeLength_*sizeof(encodedambig_t));
    }
    CUERR
}
//---------------------------------------------------------------
template<class result_type>
query_batch<result_type>::~query_batch() {
    if(maxQueries_) {
        cudaFreeHost(queryIds_);

        cudaFreeHost(h_encodeOffsets_);
        cudaFreeHost(h_queryResults_);

        cudaFree(d_encodeOffsets_);
        cudaFree(d_queryResults_);
    }
    if(maxEncodeLength_) {
        cudaFreeHost(h_encodedSeq_);
        cudaFreeHost(h_encodedAmbig_);

        cudaFree(d_encodedSeq_);
        cudaFree(d_encodedAmbig_);
    }
    CUERR
}


//---------------------------------------------------------------
template<class result_type>
void query_batch<result_type>::copy_queries_to_device(cudaStream_t stream) {
    cudaMemcpyAsync(d_encodeOffsets_, h_encodeOffsets_,
               (numQueries_+1)*sizeof(encodinglen_t), cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_encodedSeq_, h_encodedSeq_,
               h_encodeOffsets_[numQueries_]*sizeof(encodedseq_t), cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_encodedAmbig_, h_encodedAmbig_,
               h_encodeOffsets_[numQueries_]*sizeof(encodedambig_t), cudaMemcpyHostToDevice, stream);
}


//---------------------------------------------------------------
template<class result_type>
void query_batch<result_type>::copy_results_to_host(cudaStream_t stream) {
    cudaMemcpyAsync(h_queryResults_, d_queryResults_,
               numQueries_*maxResultsPerQuery_*sizeof(result_type), cudaMemcpyHostToDevice, stream);
}



//---------------------------------------------------------------
template class query_batch<location>;

} // namespace mc

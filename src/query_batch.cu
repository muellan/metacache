
#include "query_batch.cuh"
#include "sketch_database.h"

#include "../dep/cub/cub/device/device_segmented_radix_sort.cuh"

namespace mc {


//---------------------------------------------------------------
template<class result_type>
query_batch<result_type>::query_batch(
    size_t maxQueries,
    size_t maxEncodeLength,
    size_t maxResultsPerQuery) :
    numSegments_{0},
    numQueries_{0},
    maxQueries_{maxQueries},
    maxEncodeLength_{maxEncodeLength},
    maxResultsPerQuery_{maxResultsPerQuery}
{
    //TODO reuse/combine device arrays:
    // d_encodeOffsets_ + d_resultOffsets_

    if(maxQueries_) {
        cudaMallocHost(&queryIds_, maxQueries_*sizeof(query_id));

        cudaMallocHost(&h_encodeOffsets_, (maxQueries_+1)*sizeof(encodinglen_t));
        cudaMalloc    (&d_encodeOffsets_, (maxQueries_+1)*sizeof(encodinglen_t));
        h_encodeOffsets_[0] = 0;

        cudaMallocHost(&h_queryResults_, maxQueries_*maxResultsPerQuery_*sizeof(result_type));
        cudaMalloc    (&d_queryResults_, maxQueries_*maxResultsPerQuery_*sizeof(result_type));

        cudaMallocHost(&h_queryResultsTmp_, maxQueries_*maxResultsPerQuery_*sizeof(result_type));
        cudaMalloc    (&d_queryResultsTmp_, maxQueries_*maxResultsPerQuery_*sizeof(result_type));

        cudaMallocHost(&h_resultOffsets_, (maxQueries_+1)*sizeof(int));
        cudaMalloc    (&d_resultOffsets_, (maxQueries_+1)*sizeof(int));
        h_resultOffsets_[0] = 0;
    }
    CUERR

    if(maxEncodeLength_) {
        cudaMallocHost(&h_encodedSeq_, maxEncodeLength_*sizeof(encodedseq_t));
        cudaMalloc    (&d_encodedSeq_, maxEncodeLength_*sizeof(encodedseq_t));

        cudaMallocHost(&h_encodedAmbig_, maxEncodeLength_*sizeof(encodedambig_t));
        cudaMalloc    (&d_encodedAmbig_, maxEncodeLength_*sizeof(encodedambig_t));
    }
    CUERR

    cudaStreamCreate(&stream_);
    CUERR
}
//---------------------------------------------------------------
template<class result_type>
query_batch<result_type>::~query_batch() {
    if(maxQueries_) {
        cudaFreeHost(queryIds_);

        cudaFreeHost(h_encodeOffsets_);
        cudaFree    (d_encodeOffsets_);

        cudaFreeHost(h_queryResults_);
        cudaFree    (d_queryResults_);

        cudaFreeHost(h_queryResultsTmp_);
        cudaFree    (d_queryResultsTmp_);

        cudaFreeHost(h_resultOffsets_);
        cudaFree    (d_resultOffsets_);
    }
    CUERR

    if(maxEncodeLength_) {
        cudaFreeHost(h_encodedSeq_);
        cudaFree    (d_encodedSeq_);

        cudaFreeHost(h_encodedAmbig_);
        cudaFree    (d_encodedAmbig_);
    }
    CUERR

    cudaStreamDestroy(stream_);
    CUERR
}


//---------------------------------------------------------------
template<class result_type>
void query_batch<result_type>::copy_queries_to_device_async() {
    cudaMemcpyAsync(d_encodeOffsets_, h_encodeOffsets_,
                    (numQueries_+1)*sizeof(encodinglen_t),
                    cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_encodedSeq_, h_encodedSeq_,
                    h_encodeOffsets_[numQueries_]*sizeof(encodedseq_t),
                    cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_encodedAmbig_, h_encodedAmbig_,
                    h_encodeOffsets_[numQueries_]*sizeof(encodedambig_t),
                    cudaMemcpyHostToDevice, stream_);
}


//---------------------------------------------------------------
template<class result_type>
void query_batch<result_type>::copy_results_to_host_async() {
    cudaMemcpyAsync(h_queryResults_, d_queryResults_,
                    numQueries_*maxResultsPerQuery_*sizeof(result_type),
                    cudaMemcpyHostToDevice, stream_);
}

//---------------------------------------------------------------
template<class result_type>
void query_batch<result_type>::copy_results_to_host_tmp_async() {
    cudaMemcpyAsync(h_queryResultsTmp_, d_queryResults_,
                    numQueries_*maxResultsPerQuery_*sizeof(result_type),
                    cudaMemcpyHostToDevice, stream_);
}


//---------------------------------------------------------------
template<class result_type>
void query_batch<result_type>::sync_stream() {
    cudaStreamSynchronize(stream_);
}


//---------------------------------------------------------------
template<class result_type>
void query_batch<result_type>::sort_results() {

    using result_type_equivalent = uint64_t;

    static_assert(sizeof(result_type) == sizeof(result_type_equivalent), "result_type must be 64 bit");

    size_t segment = 0;
    h_resultOffsets_[segment+1] = maxResultsPerQuery_;

    for(size_t i = 1; i < numQueries_; ++i) {
        if(queryIds_[i] != queryIds_[i-1]) {
            ++segment;
        }
        h_resultOffsets_[segment+1] = (i+1) * maxResultsPerQuery_;
    }

    int numItems = numQueries_*maxResultsPerQuery_;
    int numSegments = segment+1;

    cudaMemcpyAsync(d_resultOffsets_, h_resultOffsets_,
                    (numSegments+1)*sizeof(int),
                    cudaMemcpyHostToDevice, stream_);


    cub::DoubleBuffer<result_type_equivalent> d_keys(
        (result_type_equivalent*)(d_queryResults_),
        (result_type_equivalent*)(d_queryResultsTmp_));

    // size_t tempStorageBytes = 255;
    size_t tempStorageBytes = maxEncodeLength_*sizeof(encodedseq_t);

    cudaError_t err = cub::DeviceSegmentedRadixSort::SortKeys(
        (void*)(d_encodedSeq_), tempStorageBytes,
        d_keys,
        numItems, numSegments,
        d_resultOffsets_, d_resultOffsets_ + 1,
        0, sizeof(result_type_equivalent) * CHAR_BIT,
        stream_
    );

    cudaStreamSynchronize(stream_);

    if (err != cudaSuccess) {                       \
        std::cout << "CUDA error: " << cudaGetErrorString(err) << " : "    \
                  << __FILE__ << ", line " << __LINE__ << std::endl;       \
        exit(1);                                                           \
    }

    d_queryResults_    = (result_type*)d_keys.Current();
    d_queryResultsTmp_ = (result_type*)d_keys.Alternate();
}


//---------------------------------------------------------------
template class query_batch<location>;

} // namespace mc

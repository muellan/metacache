
#include "query_batch.cuh"
#include "sketch_database.h"
#include "gpu_engine.cuh"

#include "../dep/cub/cub/device/device_segmented_radix_sort.cuh"
#include "../dep/cub/cub/device/device_scan.cuh"

namespace mc {


//---------------------------------------------------------------
template<class result_type>
query_batch<result_type>::query_batch(
    id_type maxQueries,
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

    if(maxQueries_ && maxEncodeLength_ && maxResultsPerQuery_) {
        cudaMallocHost(&h_queryIds_, maxQueries_*sizeof(id_type));
        cudaMalloc    (&d_queryIds_, maxQueries_*sizeof(id_type));

        cudaMallocHost(&h_encodeOffsets_, (maxQueries_+1)*sizeof(encodinglen_t));
        cudaMalloc    (&d_encodeOffsets_, (maxQueries_+1)*sizeof(encodinglen_t));
        h_encodeOffsets_[0] = 0;

        cudaMallocHost(&h_encodedSeq_, maxEncodeLength_*sizeof(encodedseq_t));
        cudaMalloc    (&d_encodedSeq_, maxEncodeLength_*sizeof(encodedseq_t));

        cudaMallocHost(&h_encodedAmbig_, maxEncodeLength_*sizeof(encodedambig_t));
        cudaMalloc    (&d_encodedAmbig_, maxEncodeLength_*sizeof(encodedambig_t));

        cudaMallocHost(&h_queryResults_, maxQueries_*maxResultsPerQuery_*sizeof(result_type));
        cudaMalloc    (&d_queryResults_, maxQueries_*maxResultsPerQuery_*sizeof(result_type));

        cudaMalloc    (&d_queryResultsTmp_, maxQueries_*maxResultsPerQuery_*sizeof(result_type));

        cudaMallocHost(&h_resultOffsets_, (maxQueries_+1)*sizeof(int));
        cudaMalloc    (&d_resultOffsets_, (maxQueries_+1)*sizeof(int));
        h_resultOffsets_[0] = 0;
        cudaMemcpy(d_resultOffsets_, h_resultOffsets_, sizeof(int), cudaMemcpyHostToDevice);

        cudaMalloc    (&d_resultCounts_, maxQueries_*sizeof(int));

        cudaMallocHost(&h_topCandidates_, maxQueries_*sizeof(candidate_target));
        cudaMalloc    (&d_topCandidates_, maxQueries_*sizeof(candidate_target));

        cudaMallocHost(&h_maxWindowsInRange_, maxQueries_*sizeof(window_id));
        cudaMalloc    (&d_maxWindowsInRange_, maxQueries_*sizeof(window_id));
    }
    CUERR

    cudaStreamCreate(&stream_);
    cudaEventCreate(&event_);
    CUERR
}
//---------------------------------------------------------------
template<class result_type>
query_batch<result_type>::~query_batch() {
    CUERR

    if(maxQueries_ && maxEncodeLength_ && maxResultsPerQuery_) {
        cudaFreeHost(h_queryIds_);
        cudaFree    (d_queryIds_);

        cudaFreeHost(h_encodeOffsets_);
        cudaFree    (d_encodeOffsets_);

        cudaFreeHost(h_encodedSeq_);
        cudaFree    (d_encodedSeq_);

        cudaFreeHost(h_encodedAmbig_);
        cudaFree    (d_encodedAmbig_);

        cudaFreeHost(h_queryResults_);
        cudaFree    (d_queryResults_);

        cudaFree    (d_queryResultsTmp_);

        cudaFreeHost(h_resultOffsets_);
        cudaFree    (d_resultOffsets_);

        cudaFree    (d_resultCounts_);

        cudaFreeHost(h_topCandidates_);
        cudaFree    (d_topCandidates_);

        cudaFreeHost(h_maxWindowsInRange_);
        cudaFree    (d_maxWindowsInRange_);
    }
    CUERR

    cudaStreamDestroy(stream_);
    cudaEventDestroy(event_);
    CUERR
}


//---------------------------------------------------------------
template<class result_type>
void query_batch<result_type>::copy_queries_to_device_async() {
    cudaMemcpyAsync(d_queryIds_, h_queryIds_,
                    numQueries_*sizeof(id_type),
                    cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_encodeOffsets_, h_encodeOffsets_,
                    (numQueries_+1)*sizeof(encodinglen_t),
                    cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_encodedSeq_, h_encodedSeq_,
                    h_encodeOffsets_[numQueries_]*sizeof(encodedseq_t),
                    cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_encodedAmbig_, h_encodedAmbig_,
                    h_encodeOffsets_[numQueries_]*sizeof(encodedambig_t),
                    cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_maxWindowsInRange_, h_maxWindowsInRange_,
                    numSegments_*sizeof(window_id),
                    cudaMemcpyHostToDevice, stream_);
    // cudaStreamSynchronize(stream_);
    // CUERR
}



//---------------------------------------------------------------
template<class result_type>
void query_batch<result_type>::sync_stream() {
    cudaStreamSynchronize(stream_);
}


//---------------------------------------------------------------
template<class result_type>
void query_batch<result_type>::compact_sort_and_copy_results_async() {

    {
        size_t tempStorageBytes = maxQueries_*maxResultsPerQuery_*sizeof(result_type);
        void * d_tempStorage = (void*)(d_queryResultsTmp_);

        // std::cout << "temp size: " << tempStorageBytes << std::endl;

        int numItems = numQueries_;
        // size_t tempStorageBytes = 767;
        cudaError_t err = cub::DeviceScan::InclusiveSum(
            d_tempStorage, tempStorageBytes,
            d_resultCounts_, d_resultCounts_,
            numItems,
            stream_
        );

        // cudaStreamSynchronize(stream_);
        // CUERR

        if (err != cudaSuccess) {                       \
            std::cout << "CUDA error: " << cudaGetErrorString(err) << " : "    \
                    << __FILE__ << ", line " << __LINE__ << std::endl;       \
            exit(1);                                                           \
        }

        compact_kernel<<<numQueries_,32,0,stream_>>>(
            numQueries_,
            d_resultCounts_,
            maxResultsPerQuery_,
            d_queryResults_,
            d_queryResultsTmp_,
            d_queryIds_,
            d_resultOffsets_
        );

        cudaMemcpyAsync(h_resultOffsets_, d_resultOffsets_,
                        (numSegments_+1)*sizeof(int),
                        cudaMemcpyDeviceToHost, stream_);

        cudaEventRecord(event_, stream_);

        // cudaStreamSynchronize(stream_);
        // CUERR
    }

    using result_type_equivalent = uint64_t;

    static_assert(sizeof(result_type) == sizeof(result_type_equivalent), "result_type must be 64 bit");

    cub::DoubleBuffer<result_type_equivalent> d_keys(
        (result_type_equivalent*)(d_queryResultsTmp_),
        (result_type_equivalent*)(d_queryResults_));

    size_t tempStorageBytes = maxEncodeLength_*sizeof(encodedseq_t);
    void * d_tempStorage = (void*)(d_encodedSeq_);

    int numItems = numQueries_*maxResultsPerQuery_;
    // int numItems = h_resultOffsets_[numSegments_+1];
    // size_t tempStorageBytes = 255;
    cudaError_t err = cub::DeviceSegmentedRadixSort::SortKeys(
        d_tempStorage, tempStorageBytes,
        d_keys,
        numItems, numSegments_,
        d_resultOffsets_, d_resultOffsets_ + 1,
        0, sizeof(result_type_equivalent) * CHAR_BIT,
        stream_
    );

    // cudaStreamSynchronize(stream_);
    // CUERR

    if (err != cudaSuccess) {                       \
        std::cout << "CUDA error: " << cudaGetErrorString(err) << " : "    \
                  << __FILE__ << ", line " << __LINE__ << std::endl;       \
        exit(1);                                                           \
    }

    d_queryResults_    = (result_type*)d_keys.Current();
    d_queryResultsTmp_ = (result_type*)d_keys.Alternate();


    cudaEventSynchronize(event_);

    cudaMemcpyAsync(h_queryResults_, d_queryResults_,
        // numQueries_*maxResultsPerQuery_*sizeof(result_type),
        h_resultOffsets_[numSegments_]*sizeof(result_type),
        cudaMemcpyDeviceToHost, stream_);


    //TODO max cand as template?
    #define MAX_CANDIDATES 8

    const size_t numBlocks = SDIV(numSegments_, 32);

    generate_top_candidates<MAX_CANDIDATES><<<numBlocks,32,0,stream_>>>(
        numSegments_,
        d_resultOffsets_,
        d_queryResults_,
        d_maxWindowsInRange_,
        d_topCandidates_);

    cudaStreamSynchronize(stream_);
}


//---------------------------------------------------------------
template class query_batch<location>;

} // namespace mc

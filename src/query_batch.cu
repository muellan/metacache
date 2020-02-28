
#include "query_batch.cuh"
#include "sketch_database.h"
#include "gpu_result_processing.cuh"

#include "../dep/cub/cub/device/device_scan.cuh"

#include "../dep/bb_segsort/src/bb_segsort_keys.cuh"

namespace mc {


//---------------------------------------------------------------
template<class result_type>
class query_batch<result_type>::segmented_sort
{
    using result_type_equivalent = uint64_t;

    static_assert(sizeof(result_type) == sizeof(result_type_equivalent), "result_type must be 64 bit");

public:
    segmented_sort(
        result_type_equivalent *d_keys, result_type_equivalent *d_keysB,
        const int *d_segs,
        int *d_bin_segs_id, int *d_bin_counter,
        cudaStream_t stream)
    :
        sorter_(d_keys, d_keysB,
            d_segs, d_bin_segs_id, d_bin_counter,
            stream)
    {}

    void run(int num_segs, cudaStream_t stream) const {
        sorter_.run(num_segs, stream);
    }

private:
    bb_segsort_keys<result_type_equivalent> sorter_;
};



//---------------------------------------------------------------
template<class result_type>
query_batch<result_type>::query_batch(
    id_type maxQueries,
    size_t maxSequenceLength,
    size_t maxResultsPerQuery,
    uint32_t maxCandidatesPerQuery
) :
    h_numSegments_{0},
    d_numSegments_{0},
    h_numQueries_{0},
    d_numQueries_{0},
    maxQueries_{maxQueries},
    maxSequenceLength_{maxSequenceLength},
    maxResultsPerQuery_{maxResultsPerQuery},
    maxCandidatesPerQuery_{maxCandidatesPerQuery}
{
    //TODO reuse/combine device arrays:
    // d_sequenceOffsets_ + d_resultOffsets_

    if(maxQueries_ && maxSequenceLength_ && maxResultsPerQuery_) {
        cudaMallocHost(&h_queryIds_, maxQueries_*sizeof(id_type));
        cudaMalloc    (&d_queryIds_, maxQueries_*sizeof(id_type));

        cudaMallocHost(&h_sequenceOffsets_, (maxQueries_+1)*sizeof(encodinglen_t));
        cudaMalloc    (&d_sequenceOffsets_, (maxQueries_+1)*sizeof(encodinglen_t));
        h_sequenceOffsets_[0] = 0;

        cudaMallocHost(&h_sequences_, maxSequenceLength_*sizeof(char));
        cudaMalloc    (&d_sequences_, maxSequenceLength_*sizeof(char));

        cudaMallocHost(&h_queryResults_, maxQueries_*maxResultsPerQuery_*sizeof(result_type));
        cudaMalloc    (&d_queryResults_, maxQueries_*maxResultsPerQuery_*sizeof(result_type));

        cudaMalloc    (&d_queryResultsTmp_, maxQueries_*maxResultsPerQuery_*sizeof(result_type));

        cudaMallocHost(&h_resultOffsets_, (maxQueries_+1)*sizeof(int));
        cudaMalloc    (&d_resultOffsets_, (maxQueries_+1)*sizeof(int));
        h_resultOffsets_[0] = 0;
        cudaMemcpy(d_resultOffsets_, h_resultOffsets_, sizeof(int), cudaMemcpyHostToDevice);

        cudaMalloc    (&d_resultCounts_, maxQueries_*sizeof(int));
        d_binnedSegIds_ = d_resultCounts_;

        cudaMallocHost(&h_segBinCounters_, (SEGBIN_NUM+1)*sizeof(int));
        cudaMalloc    (&d_segBinCounters_, (SEGBIN_NUM+1)*sizeof(int));

        cudaMallocHost(&h_topCandidates_, maxQueries_*maxCandidatesPerQuery_*sizeof(match_candidate));
        cudaMalloc    (&d_topCandidates_, maxQueries_*maxCandidatesPerQuery_*sizeof(match_candidate));

        cudaMallocHost(&h_maxWindowsInRange_, maxQueries_*sizeof(window_id));
        cudaMalloc    (&d_maxWindowsInRange_, maxQueries_*sizeof(window_id));
    }
    CUERR

    cudaStreamCreate(&stream_);
    cudaStreamCreate(&resultCopyStream_);

    cudaEventCreate(&queriesCopiedEvent_);
    cudaEventCreate(&offsetsCopiedEvent_);
    cudaEventCreate(&resultReadyEvent_);
    CUERR

    using result_type_equivalent = uint64_t;

    static_assert(sizeof(result_type) == sizeof(result_type_equivalent), "result_type must be 64 bit");

    sorter_ = std::make_unique<segmented_sort>(
        (result_type_equivalent*)(d_queryResultsTmp_),
        (result_type_equivalent*)(d_queryResults_),
        d_resultOffsets_,
        d_binnedSegIds_, d_segBinCounters_,
        stream_);
}
//---------------------------------------------------------------
template<class result_type>
query_batch<result_type>::~query_batch() {
    CUERR

    if(maxQueries_ && maxSequenceLength_ && maxResultsPerQuery_) {
        cudaFreeHost(h_queryIds_);
        cudaFree    (d_queryIds_);

        cudaFreeHost(h_sequenceOffsets_);
        cudaFree    (d_sequenceOffsets_);

        cudaFreeHost(h_sequences_);
        cudaFree    (d_sequences_);

        cudaFreeHost(h_queryResults_);
        cudaFree    (d_queryResults_);

        cudaFree    (d_queryResultsTmp_);

        cudaFreeHost(h_resultOffsets_);
        cudaFree    (d_resultOffsets_);

        cudaFree    (d_resultCounts_);
        d_binnedSegIds_ = nullptr;

        cudaFreeHost(h_segBinCounters_);
        cudaFree    (d_segBinCounters_);

        cudaFreeHost(h_topCandidates_);
        cudaFree    (d_topCandidates_);

        cudaFreeHost(h_maxWindowsInRange_);
        cudaFree    (d_maxWindowsInRange_);
    }
    CUERR

    cudaStreamDestroy(stream_);
    cudaStreamDestroy(resultCopyStream_);

    cudaEventDestroy(queriesCopiedEvent_);
    cudaEventDestroy(offsetsCopiedEvent_);
    cudaEventDestroy(resultReadyEvent_);
    CUERR
}


//---------------------------------------------------------------
template<class result_type>
void query_batch<result_type>::copy_queries_to_device_async() {
    d_numQueries_ = h_numQueries_;
    d_numSegments_ = h_numSegments_;

    cudaMemcpyAsync(d_queryIds_, h_queryIds_,
                    d_numQueries_*sizeof(id_type),
                    cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_sequenceOffsets_, h_sequenceOffsets_,
                    (d_numQueries_+1)*sizeof(encodinglen_t),
                    cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_sequences_, h_sequences_,
                    h_sequenceOffsets_[d_numQueries_]*sizeof(char),
                    cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(d_maxWindowsInRange_, h_maxWindowsInRange_,
                    d_numSegments_*sizeof(window_id),
                    cudaMemcpyHostToDevice, stream_);

    cudaEventRecord(queriesCopiedEvent_, stream_);

    // cudaStreamSynchronize(stream_);
    // CUERR
}

//---------------------------------------------------------------
template<class result_type>
void query_batch<result_type>::wait_for_queries_copied() {
    cudaEventSynchronize(queriesCopiedEvent_);
}



//---------------------------------------------------------------
template<class result_type>
void query_batch<result_type>::sync_streams() {
    cudaStreamSynchronize(stream_);
    cudaStreamSynchronize(resultCopyStream_);
}

//---------------------------------------------------------------
template<class result_type>
void query_batch<result_type>::sync_result_stream() {
    cudaStreamSynchronize(resultCopyStream_);
}


//---------------------------------------------------------------
template<class result_type>
void query_batch<result_type>::compact_results_async() {

    size_t tempStorageBytes = maxQueries_*maxResultsPerQuery_*sizeof(result_type);
    void * d_tempStorage = (void*)(d_queryResultsTmp_);

    cudaError_t err = cub::DeviceScan::InclusiveSum(
        d_tempStorage, tempStorageBytes,
        d_resultCounts_, d_resultCounts_,
        d_numQueries_,
        stream_
    );
    // cudaStreamSynchronize(stream_);
    // CUERR

    if (err != cudaSuccess) {                       \
        std::cout << "CUDA error: " << cudaGetErrorString(err) << " : "    \
        << __FILE__ << ", line " << __LINE__ << std::endl;       \
        exit(1);                                                           \
    }

    compact_kernel<<<d_numQueries_,32,0,stream_>>>(
        d_numQueries_,
        d_resultCounts_,
        maxResultsPerQuery_,
        d_queryResults_,
        d_queryResultsTmp_,
        d_queryIds_,
        d_resultOffsets_);
    // cudaStreamSynchronize(stream_);
    // CUERR
}



//---------------------------------------------------------------
template<class result_type>
void query_batch<result_type>::compact_sort_and_copy_results_async(bool copyAllHits)
{
    compact_results_async();

    if(copyAllHits) {
        cudaEventRecord(resultReadyEvent_, stream_);
        cudaStreamWaitEvent(resultCopyStream_, resultReadyEvent_, 0);

        cudaMemcpyAsync(h_resultOffsets_, d_resultOffsets_,
                        (d_numSegments_+1)*sizeof(int),
                        cudaMemcpyDeviceToHost, resultCopyStream_);
    }
    // cudaStreamSynchronize(resultCopyStream_);
    // CUERR

    sorter_->run(d_numSegments_, stream_);
    // cudaStreamSynchronize(stream_);
    // CUERR

    if(copyAllHits) {
        cudaEventRecord(offsetsCopiedEvent_, resultCopyStream_);

        cudaEventRecord(resultReadyEvent_, stream_);
        cudaStreamWaitEvent(resultCopyStream_, resultReadyEvent_, 0);

        cudaEventSynchronize(offsetsCopiedEvent_);

        cudaMemcpyAsync(h_queryResults_, d_queryResults_,
                        h_resultOffsets_[d_numSegments_]*sizeof(result_type),
                        cudaMemcpyDeviceToHost, resultCopyStream_);
    }
    // cudaStreamSynchronize(resultCopyStream_);
    // CUERR
}


//---------------------------------------------------------------
template<class result_type>
void query_batch<result_type>::generate_and_copy_top_candidates_async(
    const ranked_lineage * lineages,
    taxon_rank lowestRank)
{
    const size_t numBlocks = d_numSegments_;

    //TODO different max cand cases
    if(maxCandidatesPerQuery_ <= 2) {
        constexpr int maxCandidates = 2;

        generate_top_candidates<maxCandidates><<<numBlocks,32,0,stream_>>>(
            d_numSegments_,
            d_resultOffsets_,
            d_queryResults_,
            d_maxWindowsInRange_,
            lineages,
            lowestRank,
            maxCandidatesPerQuery_,
            d_topCandidates_);

        // cudaStreamSynchronize(stream_);
        // CUERR
    }

    cudaEventRecord(resultReadyEvent_, stream_);
    cudaStreamWaitEvent(resultCopyStream_, resultReadyEvent_, 0);

    // copy candidates to host
    cudaMemcpyAsync(h_topCandidates_, d_topCandidates_,
                    d_numSegments_*maxCandidatesPerQuery_*sizeof(match_candidate),
                    cudaMemcpyDeviceToHost, resultCopyStream_);

    // cudaStreamSynchronize(resultCopyStream_);
    // CUERR
}


//---------------------------------------------------------------
template class query_batch<location>;

} // namespace mc


#include "query_batch.cuh"
#include "sketch_database.h"
#include "gpu_result_processing.cuh"

#include "cub/device/device_scan.cuh"

#include "../dep/bb_segsort/src/bb_segsort_keys.cuh"

namespace mc {


//---------------------------------------------------------------
template<class Location>
class query_batch<Location>::segmented_sort
{
    using location_type_equivalent = uint64_t;

    static_assert(sizeof(location_type) == sizeof(location_type_equivalent), "location_type must be 64 bit");

public:
    segmented_sort(
        location_type_equivalent *d_keys, location_type_equivalent *d_keysB,
        const int *d_segs,
        int *d_binnedSegIds,
        int *d_segBinCounters_,
        cudaStream_t stream)
    :
        sorter_{d_keys, d_keysB,
            d_segs, d_binnedSegIds, d_segBinCounters_,
            stream}
    {}

    void run(int numSegs, cudaStream_t stream) const {
        sorter_.run(numSegs, stream);
    }

private:
    bb_segsort_keys<location_type_equivalent> sorter_;
};



//---------------------------------------------------------------
template<class Location>
query_batch<Location>::query_batch(
    index_type maxQueries,
    size_type maxSequenceLength,
    size_type maxResultsPerQuery,
    size_type maxCandidatesPerQuery
) :
    hostInput_{},
    hostOutput_{},
    gpuData_{},
    maxQueries_{maxQueries},
    maxSequenceLength_{maxSequenceLength},
    maxResultsPerQuery_{maxResultsPerQuery},
    maxCandidatesPerQuery_{maxCandidatesPerQuery}
{
    hostInput_.numSegments_ = 0;
    hostOutput_.numSegments_ = 0;
    gpuData_.numSegments_ = 0;

    hostInput_.numQueries_ = 0;
    gpuData_.numQueries_ = 0;

    size_t allocatedGpuMem = 0;

    if(maxQueries_ && maxSequenceLength_ && maxResultsPerQuery_) {
        cudaMallocHost(&hostInput_.queryIds_, maxQueries_*sizeof(index_type));
        cudaMalloc    (&gpuData_.queryIds_, maxQueries_*sizeof(index_type));
        allocatedGpuMem += maxQueries_*sizeof(index_type);

        cudaMallocHost(&hostInput_.sequenceOffsets_, (maxQueries_+1)*sizeof(size_type));
        cudaMalloc    (&gpuData_.sequenceOffsets_, (maxQueries_+1)*sizeof(size_type));
        allocatedGpuMem += (maxQueries_+1)*sizeof(size_type);
        hostInput_.sequenceOffsets_[0] = 0;

        cudaMallocHost(&hostInput_.sequences_, maxSequenceLength_*sizeof(char));
        cudaMalloc    (&gpuData_.sequences_, maxSequenceLength_*sizeof(char));
        allocatedGpuMem += maxSequenceLength_*sizeof(char);

        cudaMallocHost(&hostOutput_.queryResults_, maxQueries_*maxResultsPerQuery_*sizeof(location_type));
        cudaMalloc    (&gpuData_.queryResults_, maxQueries_*maxResultsPerQuery_*sizeof(location_type));
        cudaMalloc    (&gpuData_.queryResultsTmp_, maxQueries_*maxResultsPerQuery_*sizeof(location_type));
        allocatedGpuMem += 2*maxQueries_*maxResultsPerQuery_*sizeof(location_type);

        cudaMallocHost(&hostOutput_.resultOffsets_, (maxQueries_+1)*sizeof(int));
        cudaMalloc    (&gpuData_.resultOffsets_, (maxQueries_+1)*sizeof(int));
        allocatedGpuMem += (maxQueries_+1)*sizeof(int);
        hostOutput_.resultOffsets_[0] = 0;
        cudaMemcpy(gpuData_.resultOffsets_, hostOutput_.resultOffsets_, sizeof(int), cudaMemcpyHostToDevice);

        cudaMalloc    (&gpuData_.resultCounts_, maxQueries_*sizeof(int));
        allocatedGpuMem += maxQueries_*sizeof(int);

        cudaMalloc    (&gpuData_.segBinCounters_, (SEGBIN_NUM+1)*sizeof(int));
        allocatedGpuMem += (SEGBIN_NUM+1)*sizeof(int);

        cudaMallocHost(&hostOutput_.topCandidates_, maxQueries_*maxCandidatesPerQuery_*sizeof(match_candidate));
        cudaMalloc    (&gpuData_.topCandidates_, maxQueries_*maxCandidatesPerQuery_*sizeof(match_candidate));
        allocatedGpuMem += maxQueries_*maxCandidatesPerQuery_*sizeof(match_candidate);

        cudaMallocHost(&hostInput_.maxWindowsInRange_, maxQueries_*sizeof(window_id));
        cudaMalloc    (&gpuData_.maxWindowsInRange_, maxQueries_*sizeof(window_id));
        allocatedGpuMem += maxQueries_*sizeof(window_id);
    }
    CUERR

    // std::cerr << "query batch size on gpu: " << (allocatedGpuMem >> 20) << " MB\n";

    cudaStreamCreate(&stream_);
    cudaStreamCreate(&resultCopyStream_);

    cudaEventCreate(&queriesCopiedEvent_);
    cudaEventCreate(&offsetsCopiedEvent_);
    cudaEventCreate(&resultReadyEvent_);
    CUERR

    using location_type_equivalent = uint64_t;

    static_assert(sizeof(location_type) == sizeof(location_type_equivalent), "location_type must be 64 bit");

    sorter_ = std::make_unique<segmented_sort>(
        (location_type_equivalent*)(gpuData_.queryResultsTmp_),
        (location_type_equivalent*)(gpuData_.queryResults_),
        gpuData_.resultOffsets_,
        gpuData_.resultCounts_, // reuse for binning
        gpuData_.segBinCounters_,
        stream_);
    CUERR
}
//---------------------------------------------------------------
template<class Location>
query_batch<Location>::~query_batch() {
    CUERR

    if(maxQueries_ && maxSequenceLength_ && maxResultsPerQuery_) {
        cudaFreeHost(hostInput_.queryIds_);
        cudaFree    (gpuData_.queryIds_);

        cudaFreeHost(hostInput_.sequenceOffsets_);
        cudaFree    (gpuData_.sequenceOffsets_);

        cudaFreeHost(hostInput_.sequences_);
        cudaFree    (gpuData_.sequences_);

        cudaFreeHost(hostOutput_.queryResults_);
        cudaFree    (gpuData_.queryResults_);

        cudaFree    (gpuData_.queryResultsTmp_);

        cudaFreeHost(hostOutput_.resultOffsets_);
        cudaFree    (gpuData_.resultOffsets_);

        cudaFree    (gpuData_.resultCounts_);

        cudaFree    (gpuData_.segBinCounters_);

        cudaFreeHost(hostOutput_.topCandidates_);
        cudaFree    (gpuData_.topCandidates_);

        cudaFreeHost(hostInput_.maxWindowsInRange_);
        cudaFree    (gpuData_.maxWindowsInRange_);
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
template<class Location>
void query_batch<Location>::copy_queries_to_device_async() {
    gpuData_.numQueries_ = hostInput_.numQueries_;
    gpuData_.numSegments_ = hostInput_.numSegments_;
    hostOutput_.numSegments_ = hostInput_.numSegments_;

    cudaMemcpyAsync(gpuData_.queryIds_, hostInput_.queryIds_,
                    gpuData_.numQueries_*sizeof(index_type),
                    cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(gpuData_.sequenceOffsets_, hostInput_.sequenceOffsets_,
                    (gpuData_.numQueries_+1)*sizeof(size_type),
                    cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(gpuData_.sequences_, hostInput_.sequences_,
                    hostInput_.sequenceOffsets_[gpuData_.numQueries_]*sizeof(char),
                    cudaMemcpyHostToDevice, stream_);
    cudaMemcpyAsync(gpuData_.maxWindowsInRange_, hostInput_.maxWindowsInRange_,
                    gpuData_.numSegments_*sizeof(window_id),
                    cudaMemcpyHostToDevice, stream_);

    cudaEventRecord(queriesCopiedEvent_, stream_);

    // cudaStreamSynchronize(stream_);
    // CUERR
}

//---------------------------------------------------------------
template<class Location>
void query_batch<Location>::wait_for_queries_copied() {
    cudaEventSynchronize(queriesCopiedEvent_);
}



//---------------------------------------------------------------
template<class Location>
void query_batch<Location>::sync_streams() {
    cudaStreamSynchronize(stream_);
    cudaStreamSynchronize(resultCopyStream_);
}

//---------------------------------------------------------------
template<class Location>
void query_batch<Location>::sync_result_stream() {
    cudaStreamSynchronize(resultCopyStream_);
}


//---------------------------------------------------------------
template<class Location>
void query_batch<Location>::compact_results_async() {

    size_t tempStorageBytes = maxQueries_*maxResultsPerQuery_*sizeof(location_type);
    void * d_tempStorage = (void*)(gpuData_.queryResultsTmp_);

    cudaError_t err = cub::DeviceScan::InclusiveSum(
        d_tempStorage, tempStorageBytes,
        gpuData_.resultCounts_, gpuData_.resultCounts_,
        gpuData_.numQueries_,
        stream_
    );
    // cudaStreamSynchronize(stream_);
    // CUERR

    if (err != cudaSuccess) {                       \
        std::cout << "CUDA error: " << cudaGetErrorString(err) << " : "    \
        << __FILE__ << ", line " << __LINE__ << std::endl;       \
        exit(1);                                                           \
    }

    compact_kernel<<<gpuData_.numQueries_,128,0,stream_>>>(
        gpuData_.numQueries_,
        gpuData_.resultCounts_,
        maxResultsPerQuery_,
        gpuData_.queryResults_,
        gpuData_.queryResultsTmp_,
        gpuData_.queryIds_,
        gpuData_.resultOffsets_);
    // cudaStreamSynchronize(stream_);
    // CUERR
}



//---------------------------------------------------------------
template<class Location>
void query_batch<Location>::compact_sort_and_copy_allhits_async(bool copyAllHits)
{
    compact_results_async();

    if(copyAllHits) {
        cudaEventRecord(resultReadyEvent_, stream_);
        cudaStreamWaitEvent(resultCopyStream_, resultReadyEvent_, 0);

        cudaMemcpyAsync(hostOutput_.resultOffsets_, gpuData_.resultOffsets_,
                        (gpuData_.numSegments_+1)*sizeof(int),
                        cudaMemcpyDeviceToHost, resultCopyStream_);
    }
    // cudaStreamSynchronize(resultCopyStream_);
    // CUERR

    sorter_->run(gpuData_.numSegments_, stream_);
    // cudaStreamSynchronize(stream_);
    // CUERR

    if(copyAllHits) {
        cudaEventRecord(offsetsCopiedEvent_, resultCopyStream_);

        cudaEventRecord(resultReadyEvent_, stream_);
        cudaStreamWaitEvent(resultCopyStream_, resultReadyEvent_, 0);

        cudaEventSynchronize(offsetsCopiedEvent_);

        cudaMemcpyAsync(hostOutput_.queryResults_, gpuData_.queryResults_,
                        hostOutput_.resultOffsets_[gpuData_.numSegments_]*sizeof(location_type),
                        cudaMemcpyDeviceToHost, resultCopyStream_);
    }
    // cudaStreamSynchronize(resultCopyStream_);
    // CUERR
}


//---------------------------------------------------------------
template<class Location>
void query_batch<Location>::generate_and_copy_top_candidates_async(
    const ranked_lineage * lineages,
    taxon_rank lowestRank)
{
    const index_type numBlocks = gpuData_.numSegments_;

    //TODO different max cand cases
    if(maxCandidatesPerQuery_ <= 2) {
        constexpr int maxCandidates = 2;

        generate_top_candidates<maxCandidates><<<numBlocks,32,0,stream_>>>(
            gpuData_.numSegments_,
            gpuData_.resultOffsets_,
            gpuData_.queryResults_,
            gpuData_.maxWindowsInRange_,
            lineages,
            lowestRank,
            maxCandidatesPerQuery_,
            gpuData_.topCandidates_);

        // cudaStreamSynchronize(stream_);
        // CUERR
    }

    cudaEventRecord(resultReadyEvent_, stream_);
    cudaStreamWaitEvent(resultCopyStream_, resultReadyEvent_, 0);

    // copy candidates to host
    cudaMemcpyAsync(hostOutput_.topCandidates_, gpuData_.topCandidates_,
                    gpuData_.numSegments_*maxCandidatesPerQuery_*sizeof(match_candidate),
                    cudaMemcpyDeviceToHost, resultCopyStream_);

    // cudaStreamSynchronize(resultCopyStream_);
    // CUERR
}


//---------------------------------------------------------------
template class query_batch<location>;

} // namespace mc

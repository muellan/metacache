/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2020 Robin Kobus  (kobus@uni-mainz.de)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *****************************************************************************/

#include "query_batch.cuh"
#include "database.h"
#include "gpu_result_processing.cuh"

#include <cub/device/device_scan.cuh>

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
        int *d_segBinCounters,
        cudaStream_t stream)
    :
        sorter_{d_keys, d_keysB,
            d_segs, d_binnedSegIds, d_segBinCounters,
            stream}
    {}

    void run(int numSegs, int maxSegmentSize, cudaStream_t stream) const {
        sorter_.run(numSegs, maxSegmentSize, stream);
    }

private:
    bb_segsort_keys<location_type_equivalent, int> sorter_;
};


//---------------------------------------------------------------
template<class Location>
query_batch<Location>::query_host_data::query_host_data(
    index_type maxQueries,
    size_type maxSequenceLength,
    size_type maxResultsPerWindow,
    size_type maxCandidatesPerQuery,
    bool copyAllHits
) :
    numQueries_{0},
    numWindows_{0},
    largestSegmentSize_{0},
    maxCandidatesPerQuery_{maxCandidatesPerQuery}
{
    cudaMallocHost(&queryIds_, maxQueries*sizeof(index_type));
    cudaMallocHost(&sequenceOffsets_, (maxQueries+1)*sizeof(size_type));
    sequenceOffsets_[0] = 0;
    cudaMallocHost(&sequences_, maxSequenceLength*sizeof(char));
    cudaMallocHost(&maxWindowsInRange_, maxQueries*sizeof(window_id));
    CUERR
    if(copyAllHits)
        cudaMallocHost(&queryResults_, maxQueries*maxResultsPerWindow*sizeof(location_type));
    else
        queryResults_ = nullptr;
    cudaMallocHost(&resultOffsets_, (maxQueries+1)*sizeof(int));
    resultOffsets_[0] = 0;
    cudaMallocHost(&topCandidates_, maxQueries*maxCandidatesPerQuery_*sizeof(match_candidate));
    CUERR

    cudaEventCreate(&resultsCopiedEvent_);
    CUERR
}
//---------------------------------------------------------------
template<class Location>
query_batch<Location>::query_host_data::~query_host_data()
{
    if(queryIds_)          cudaFreeHost(queryIds_);
    if(sequenceOffsets_)   cudaFreeHost(sequenceOffsets_);
    if(sequences_)         cudaFreeHost(sequences_);
    if(maxWindowsInRange_) cudaFreeHost(maxWindowsInRange_);
    CUERR
    if(queryResults_)  cudaFreeHost(queryResults_);
    if(resultOffsets_) cudaFreeHost(resultOffsets_);
    if(topCandidates_) cudaFreeHost(topCandidates_);
    CUERR

    if(resultsCopiedEvent_)  cudaEventDestroy(resultsCopiedEvent_);
    CUERR
}
//---------------------------------------------------------------
template<class Location>
query_batch<Location>::query_host_data::query_host_data(query_host_data&& other)
{
    numQueries_ = other.numQueries_;
    numWindows_  = other.numWindows_;
    largestSegmentSize_ = other.largestSegmentSize_;
    maxCandidatesPerQuery_  = other.maxCandidatesPerQuery_;

    queryIds_       = other.queryIds_;
    other.queryIds_ = nullptr;
    sequenceOffsets_       = other.sequenceOffsets_;
    other.sequenceOffsets_ = nullptr;
    sequences_       = other.sequences_;
    other.sequences_ = nullptr;
    maxWindowsInRange_       = other.maxWindowsInRange_;
    other.maxWindowsInRange_ = nullptr;

    queryResults_       = other.queryResults_;
    other.queryResults_ = nullptr;
    resultOffsets_       = other.resultOffsets_;
    other.resultOffsets_ = nullptr;
    topCandidates_       = other.topCandidates_;
    other.topCandidates_ = nullptr;

    resultsCopiedEvent_ = other.resultsCopiedEvent_;
    other.resultsCopiedEvent_ = 0;
}
//---------------------------------------------------------------
template<class Location>
void query_batch<Location>::query_host_data::wait_for_results()
{
    cudaEventSynchronize(resultsCopiedEvent_);
    CUERR
}



//---------------------------------------------------------------
template<class Location>
query_batch<Location>::query_gpu_data::query_gpu_data(
    index_type maxQueries,
    size_type maxSequenceLength,
    size_type maxSketchSize,
    size_type maxResultsPerWindow,
    size_type maxCandidatesPerQuery,
    bool multiGPU,
    part_id gpuId
)
{
    size_t allocatedGpuMem = 0;

    cudaSetDevice(gpuId); CUERR

    cudaMalloc    (&queryIds_, maxQueries*sizeof(index_type));
    allocatedGpuMem += maxQueries*sizeof(index_type);
    if(gpuId == 0) {
        // only first gpu holds sequences
        // other gpus use sketches generated by first gpu
        cudaMalloc    (&sequenceOffsets_, (maxQueries+1)*sizeof(size_type));
        allocatedGpuMem += (maxQueries+1)*sizeof(size_type);
        cudaMalloc    (&sequences_, maxSequenceLength*sizeof(char));
        allocatedGpuMem += maxSequenceLength*sizeof(char);
    }
    else {
        sequenceOffsets_ = nullptr;
        sequences_ = nullptr;
    }

    if(multiGPU) {
        cudaMalloc    (&sketches_, maxQueries*maxSketchSize*sizeof(feature_type));
        allocatedGpuMem += maxQueries*maxSketchSize*sizeof(feature_type);
    }
    else {
        sketches_ = nullptr;
    }
    cudaMalloc    (&queryResults_, maxQueries*maxResultsPerWindow*sizeof(location_type));
    cudaMalloc    (&queryResultsTmp_, maxQueries*maxResultsPerWindow*sizeof(location_type));
    allocatedGpuMem += 2*maxQueries*maxResultsPerWindow*sizeof(location_type);
    cudaMalloc    (&resultOffsets_, (maxQueries+1)*sizeof(int));
    allocatedGpuMem += (maxQueries+1)*sizeof(int);
    cudaMemset(resultOffsets_, 0, sizeof(int));
    cudaMalloc    (&resultCounts_, maxQueries*sizeof(int));
    allocatedGpuMem += maxQueries*sizeof(int);
    cudaMalloc    (&segBinCounters_, (SEGBIN_NUM+1)*sizeof(int));
    allocatedGpuMem += (SEGBIN_NUM+1)*sizeof(int);
    cudaMalloc    (&topCandidates_, maxQueries*maxCandidatesPerQuery*sizeof(match_candidate));
    allocatedGpuMem += maxQueries*maxCandidatesPerQuery*sizeof(match_candidate);
    cudaMalloc    (&maxWindowsInRange_, maxQueries*sizeof(window_id));
    allocatedGpuMem += maxQueries*sizeof(window_id);
    CUERR

    // std::cerr << "query batch size on gpu: " << (allocatedGpuMem >> 20) << " MB\n";

    cudaStreamCreate(&workStream_);
    cudaStreamCreate(&copyStream_);
    CUERR
    cudaEventCreate(&queryFinishedEvent_);
    cudaEventCreate(&sketchesCopiedEvent_);
    cudaEventCreate(&queryIdsCopiedEvent_);
    cudaEventCreate(&queryIdsFinishedEvent_);
    cudaEventCreate(&offsetsReadyEvent_);
    cudaEventCreate(&offsetsCopiedEvent_);
    cudaEventCreate(&allhitsReadyEvent_);
    cudaEventCreate(&allhitsCopiedEvent_);
    cudaEventCreate(&tophitsReadyEvent_);
    cudaEventCreate(&tophitsCopiedEvent_);
    CUERR
}
//---------------------------------------------------------------
template<class Location>
query_batch<Location>::query_gpu_data::query_gpu_data(query_gpu_data&& other)
{
    queryIds_ = other.queryIds_;
    other.queryIds_ = nullptr;
    sequenceOffsets_ = other.sequenceOffsets_;
    other.sequenceOffsets_ = nullptr;
    sequences_ = other.sequences_;
    other.sequences_ = nullptr;
    sketches_ = other.sketches_;
    other.sketches_ = nullptr;
    queryResults_ = other.queryResults_;
    other.queryResults_ = nullptr;
    queryResultsTmp_ = other.queryResultsTmp_;
    other.queryResultsTmp_ = nullptr;
    resultOffsets_ = other.resultOffsets_;
    other.resultOffsets_ = nullptr;
    resultCounts_ = other.resultCounts_;
    other.resultCounts_ = nullptr;
    segBinCounters_ = other.segBinCounters_;
    other.segBinCounters_ = nullptr;
    topCandidates_ = other.topCandidates_;
    other.topCandidates_ = nullptr;
    maxWindowsInRange_ = other.maxWindowsInRange_;
    other.maxWindowsInRange_ = nullptr;

    workStream_ = other.workStream_;
    other.workStream_ = 0;
    copyStream_ = other.copyStream_;
    other.copyStream_ = 0;

    queryFinishedEvent_ = other.queryFinishedEvent_;
    other.queryFinishedEvent_ = 0;
    sketchesCopiedEvent_ = other.sketchesCopiedEvent_;
    other.sketchesCopiedEvent_ = 0;
    queryIdsCopiedEvent_ = other.queryIdsCopiedEvent_;
    other.queryIdsCopiedEvent_ = 0;
    queryIdsFinishedEvent_ = other.queryIdsFinishedEvent_;
    other.queryIdsFinishedEvent_ = 0;
    offsetsReadyEvent_ = other.offsetsReadyEvent_;
    other.offsetsReadyEvent_ = 0;
    offsetsCopiedEvent_ = other.offsetsCopiedEvent_;
    other.offsetsCopiedEvent_ = 0;
    allhitsReadyEvent_ = other.allhitsReadyEvent_;
    other.allhitsReadyEvent_ = 0;
    allhitsCopiedEvent_ = other.allhitsCopiedEvent_;
    other.allhitsCopiedEvent_ = 0;
    tophitsReadyEvent_ = other.tophitsReadyEvent_;
    other.tophitsReadyEvent_ = 0;
    tophitsCopiedEvent_ = other.tophitsCopiedEvent_;
    other.tophitsCopiedEvent_ = 0;
}
//---------------------------------------------------------------
template<class Location>
query_batch<Location>::query_gpu_data::~query_gpu_data()
{
    if(queryIds_)          cudaFree    (queryIds_);
    if(sequenceOffsets_)   cudaFree    (sequenceOffsets_);
    if(sequences_)         cudaFree    (sequences_);
    if(sketches_)          cudaFree    (sketches_);
    if(queryResults_)      cudaFree    (queryResults_);
    if(queryResultsTmp_)   cudaFree    (queryResultsTmp_);
    if(resultOffsets_)     cudaFree    (resultOffsets_);
    if(resultCounts_)      cudaFree    (resultCounts_);
    if(segBinCounters_)    cudaFree    (segBinCounters_);
    if(topCandidates_)     cudaFree    (topCandidates_);
    if(maxWindowsInRange_) cudaFree    (maxWindowsInRange_);
    CUERR
    if(workStream_) cudaStreamDestroy(workStream_);
    if(copyStream_) cudaStreamDestroy(copyStream_);
    CUERR
    if(queryFinishedEvent_)  cudaEventDestroy(queryFinishedEvent_);
    if(sketchesCopiedEvent_) cudaEventDestroy(sketchesCopiedEvent_);
    if(queryIdsCopiedEvent_)  cudaEventDestroy(queryIdsCopiedEvent_);
    if(queryIdsFinishedEvent_) cudaEventDestroy(queryIdsFinishedEvent_);
    if(offsetsReadyEvent_)  cudaEventDestroy(offsetsReadyEvent_);
    if(offsetsCopiedEvent_) cudaEventDestroy(offsetsCopiedEvent_);
    if(allhitsReadyEvent_)   cudaEventDestroy(allhitsReadyEvent_);
    if(allhitsCopiedEvent_)  cudaEventDestroy(allhitsCopiedEvent_);
    if(tophitsReadyEvent_)   cudaEventDestroy(tophitsReadyEvent_);
    if(tophitsCopiedEvent_)  cudaEventDestroy(tophitsCopiedEvent_);
    CUERR
}


//---------------------------------------------------------------
template<class Location>
query_batch<Location>::query_batch(
    index_type maxQueries,
    size_type maxSequenceLength,
    size_type maxSketchSize,
    size_type maxResultsPerWindow,
    size_type maxCandidatesPerQuery,
    bool copyAllHits,
    part_id numHostThreads,
    part_id numGPUs
) :
    maxWindows_{maxQueries},
    maxSequenceLength_{maxSequenceLength},
    maxSketchSize_{maxSketchSize},
    maxResultsPerWindow_{maxResultsPerWindow},
    maxCandidatesPerQuery_{maxCandidatesPerQuery},
    hostData_{},
    gpuData_{},
    sorters_{},
    numGPUs_{numGPUs}
{
    cudaStreamCreate(&h2dCopyStream_);

    cudaEventCreate(&queriesCopiedEvent_);
    cudaEventCreate(&queryIdsCopiedEvent_);
    cudaEventCreate(&maxWinCopiedEvent_);
    CUERR

    cudaSetDevice(numGPUs_-1); CUERR

    for(part_id hostId = 0; hostId < numHostThreads; ++hostId) {
        hostData_.emplace_back(maxQueries, maxSequenceLength, maxResultsPerWindow, maxCandidatesPerQuery, copyAllHits);
    }

    using location_type_equivalent = uint64_t;
    static_assert(sizeof(location_type) == sizeof(location_type_equivalent), "location_type must be 64 bit");

    for(part_id gpuId = 0; gpuId < numGPUs_; ++gpuId) {
        cudaSetDevice(gpuId); CUERR

        gpuData_.emplace_back(maxQueries, maxSequenceLength, maxSketchSize, maxResultsPerWindow, maxCandidatesPerQuery, numGPUs > 1, gpuId);

        sorters_.emplace_back(
            (location_type_equivalent*)(gpuData_[gpuId].queryResultsTmp_),
            (location_type_equivalent*)(gpuData_[gpuId].queryResults_),
            gpuData_[gpuId].resultOffsets_,
            gpuData_[gpuId].resultCounts_, // reuse for binning
            gpuData_[gpuId].segBinCounters_,
            gpuData_[gpuId].workStream_);
        CUERR
    }
}
//---------------------------------------------------------------
template<class Location>
query_batch<Location>::~query_batch()
{
    cudaStreamDestroy(h2dCopyStream_);
    CUERR
    cudaEventDestroy(queriesCopiedEvent_);
    cudaEventDestroy(queryIdsCopiedEvent_);
    cudaEventDestroy(maxWinCopiedEvent_);
    CUERR
}


//---------------------------------------------------------------
template<class Location>
void query_batch<Location>::copy_queries_to_device_async(part_id hostId)
{
    // copy from host to device 0
    part_id gpuId = 0;

    cudaSetDevice(gpuId);

    cudaStreamWaitEvent(h2dCopyStream_, gpuData_[gpuId].queryFinishedEvent_, 0);

    cudaMemcpyAsync(gpuData_[gpuId].sequenceOffsets_, hostData_[hostId].sequence_offsets(),
                    (hostData_[hostId].num_windows()+1)*sizeof(size_type),
                    cudaMemcpyHostToDevice, h2dCopyStream_);
    cudaMemcpyAsync(gpuData_[gpuId].sequences_, hostData_[hostId].sequences(),
                    hostData_[hostId].sequence_offsets()[hostData_[hostId].num_windows()]*sizeof(char),
                    cudaMemcpyHostToDevice, h2dCopyStream_);

    cudaEventRecord(queriesCopiedEvent_, h2dCopyStream_);

    cudaStreamWaitEvent(gpuData_[gpuId].workStream_, queriesCopiedEvent_, 0);
    cudaStreamWaitEvent(gpuData_[gpuId].workStream_, gpuData_[gpuId].sketchesCopiedEvent_, 0);

    // cudaEventSynchronize(queriesCopiedEvent_);
    // CUERR


    cudaStreamWaitEvent(h2dCopyStream_, gpuData_[gpuId].queryIdsCopiedEvent_, 0);
    cudaStreamWaitEvent(h2dCopyStream_, gpuData_[gpuId].queryIdsFinishedEvent_, 0);

    cudaMemcpyAsync(gpuData_[gpuId].queryIds_, hostData_[hostId].query_ids(),
                    hostData_[hostId].num_windows()*sizeof(index_type),
                    cudaMemcpyHostToDevice, h2dCopyStream_);

    cudaEventRecord(queryIdsCopiedEvent_, h2dCopyStream_);

    // cudaEventSynchronize(queryIdsCopiedEvent_);
    // CUERR


    cudaStreamWaitEvent(h2dCopyStream_, gpuData_[gpuId].tophitsCopiedEvent_, 0);

    cudaMemcpyAsync(gpuData_[gpuId].maxWindowsInRange_, hostData_[hostId].max_windows_in_range(),
                    hostData_[hostId].num_queries()*sizeof(window_id),
                    cudaMemcpyHostToDevice, h2dCopyStream_);

    cudaEventRecord(maxWinCopiedEvent_, h2dCopyStream_);

    // cudaEventSynchronize(maxWinCopiedEvent_);
    // CUERR
}


//---------------------------------------------------------------
template<class Location>
void query_batch<Location>::copy_queries_to_next_device_async(part_id hostId, part_id gpuId)
{
    // copy from device gpuId to device gpuId+1
    cudaSetDevice(gpuId);

    if(gpuId == 0)
        cudaStreamWaitEvent(gpuData_[gpuId].copyStream_, gpuData_[gpuId].queryFinishedEvent_, 0);

    cudaStreamWaitEvent(gpuData_[gpuId].copyStream_, gpuData_[gpuId+1].queryFinishedEvent_, 0);
    cudaStreamWaitEvent(gpuData_[gpuId].copyStream_, gpuData_[gpuId+1].sketchesCopiedEvent_, 0);


    cudaMemcpyPeerAsync(gpuData_[gpuId+1].sketches_, gpuId+1,
                        gpuData_[gpuId].sketches_, gpuId,
                        hostData_[hostId].num_windows()*maxSketchSize_*sizeof(feature_type),
                        gpuData_[gpuId].copyStream_);

    cudaEventRecord(gpuData_[gpuId].sketchesCopiedEvent_, gpuData_[gpuId].copyStream_);

    // cudaEventSynchronize(gpuData_[gpuId].sketchesCopiedEvent_);
    // CUERR

    cudaSetDevice(gpuId+1);

    cudaStreamWaitEvent(gpuData_[gpuId+1].copyStream_, gpuData_[gpuId].sketchesCopiedEvent_, 0);
    cudaStreamWaitEvent(gpuData_[gpuId+1].workStream_, gpuData_[gpuId].sketchesCopiedEvent_, 0);


    cudaSetDevice(gpuId);

    if(gpuId == 0)
        cudaStreamWaitEvent(gpuData_[gpuId].copyStream_, queryIdsCopiedEvent_, 0);
    else
        cudaStreamWaitEvent(gpuData_[gpuId].copyStream_, gpuData_[gpuId-1].queryIdsCopiedEvent_, 0);

    cudaStreamWaitEvent(gpuData_[gpuId].copyStream_, gpuData_[gpuId+1].queryIdsCopiedEvent_, 0);
    cudaStreamWaitEvent(gpuData_[gpuId].copyStream_, gpuData_[gpuId+1].queryIdsFinishedEvent_, 0);

    cudaMemcpyPeerAsync(gpuData_[gpuId+1].queryIds_, gpuId+1,
                        gpuData_[gpuId].queryIds_, gpuId,
                        hostData_[hostId].num_windows()*sizeof(index_type),
                        gpuData_[gpuId].copyStream_);

    cudaEventRecord(gpuData_[gpuId].queryIdsCopiedEvent_, gpuData_[gpuId].copyStream_);

    // cudaEventSynchronize(gpuData_[gpuId].queryIdsCopiedEvent_);
    // CUERR
}


//---------------------------------------------------------------
template<class Location>
void query_batch<Location>::wait_for_queries_copied()
{
    // cudaEventSynchronize(queriesCopiedEvent_);
    cudaEventSynchronize(maxWinCopiedEvent_);
    CUERR
}


//---------------------------------------------------------------
template<class Location>
void query_batch<Location>::mark_query_finished(part_id gpuId)
{
    cudaEventRecord(gpuData_[gpuId].queryFinishedEvent_, gpuData_[gpuId].workStream_);
}


//---------------------------------------------------------------
template<class Location>
void query_batch<Location>::compact_results_async(part_id hostId, part_id gpuId)
{
    cudaSetDevice(gpuId);

    if(gpuId == 0)
        cudaStreamWaitEvent(gpuData_[gpuId].workStream_, queryIdsCopiedEvent_, 0);
    else
        cudaStreamWaitEvent(gpuData_[gpuId].workStream_, gpuData_[gpuId-1].queryIdsCopiedEvent_, 0);

    size_t tempStorageBytes = maxWindows_*maxResultsPerWindow_*sizeof(location_type);
    void * d_tempStorage = (void*)(gpuData_[gpuId].queryResultsTmp_);

    cudaError_t err = cub::DeviceScan::InclusiveSum(
        d_tempStorage, tempStorageBytes,
        gpuData_[gpuId].resultCounts_, gpuData_[gpuId].resultCounts_,
        hostData_[hostId].num_windows(),
        gpuData_[gpuId].workStream_
    );
    // cudaStreamSynchronize(gpuData_[gpuId].workStream_);
    // CUERR

    if (err != cudaSuccess) {                       \
        std::cout << "CUDA error: " << cudaGetErrorString(err) << " : "    \
        << __FILE__ << ", line " << __LINE__ << std::endl;       \
        exit(1);                                                           \
    }

    compact_results<<<hostData_[hostId].num_windows(),128,0,gpuData_[gpuId].workStream_>>>(
        hostData_[hostId].num_windows(),
        gpuData_[gpuId].resultCounts_,
        maxResultsPerWindow_,
        gpuData_[gpuId].queryResults_,
        gpuData_[gpuId].queryResultsTmp_,
        gpuData_[gpuId].queryIds_,
        gpuData_[gpuId].resultOffsets_);

    cudaEventRecord(gpuData_[gpuId].queryIdsFinishedEvent_, gpuData_[gpuId].workStream_);

    // cudaEventSynchronize(gpuData_[gpuId].queryIdsFinishedEvent_);
    // CUERR
}


//---------------------------------------------------------------
template<class Location>
void query_batch<Location>::compact_sort_and_copy_allhits_async(
    part_id hostId,
    part_id gpuId,
    bool copyAllHits)
{
    cudaSetDevice(gpuId);

    compact_results_async(hostId, gpuId);
    // cudaStreamSynchronize(gpuData_[gpuId].workStream_);
    // CUERR

    if(copyAllHits) {
        cudaEventRecord(gpuData_[gpuId].offsetsReadyEvent_, gpuData_[gpuId].workStream_);
        cudaStreamWaitEvent(gpuData_[gpuId].copyStream_, gpuData_[gpuId].offsetsReadyEvent_, 0);

        cudaMemcpyAsync(hostData_[hostId].result_offsets(), gpuData_[gpuId].resultOffsets_,
                        (hostData_[hostId].num_queries()+1)*sizeof(int),
                        cudaMemcpyDeviceToHost, gpuData_[gpuId].copyStream_);

        cudaEventRecord(gpuData_[gpuId].offsetsCopiedEvent_, gpuData_[gpuId].copyStream_);

        // cudaStreamSynchronize(gpuData_[gpuId].copyStream_);
        // CUERR
    }

    sorters_[gpuId].run(hostData_[hostId].num_queries(), hostData_[hostId].largest_segment_size(), gpuData_[gpuId].workStream_);
    // cudaStreamSynchronize(gpuData_[gpuId].workStream_);
    // CUERR

    if(copyAllHits) {
        cudaEventRecord(gpuData_[gpuId].allhitsReadyEvent_, gpuData_[gpuId].workStream_);
        cudaStreamWaitEvent(gpuData_[gpuId].copyStream_, gpuData_[gpuId].allhitsReadyEvent_, 0);

        cudaEventSynchronize(gpuData_[gpuId].offsetsCopiedEvent_);

        cudaMemcpyAsync(hostData_[hostId].query_results(), gpuData_[gpuId].queryResults_,
                        hostData_[hostId].result_offsets()[hostData_[hostId].num_queries()]*sizeof(location_type),
                        cudaMemcpyDeviceToHost, gpuData_[gpuId].copyStream_);

        cudaEventRecord(gpuData_[gpuId].allhitsCopiedEvent_, gpuData_[gpuId].copyStream_);

        // cudaStreamSynchronize(gpuData_[gpuId].copyStream_);
        // CUERR
    }
}


//---------------------------------------------------------------
template<class Location>
void query_batch<Location>::generate_and_copy_top_candidates_async(
    part_id hostId,
    part_id gpuId,
    const ranked_lineage * lineages,
    taxon_rank lowestRank)
{
    cudaSetDevice(gpuId);

    cudaEvent_t& event = (gpuId == 0) ? maxWinCopiedEvent_ : gpuData_[gpuId-1].tophitsCopiedEvent_;

    cudaStreamWaitEvent(gpuData_[gpuId].workStream_, event, 0);
    if(gpuId+1 < numGPUs_) {
        // copy maxWindowsInRange to next device
        cudaStreamWaitEvent(gpuData_[gpuId].copyStream_, event, 0);
        cudaStreamWaitEvent(gpuData_[gpuId].copyStream_, gpuData_[gpuId+1].tophitsCopiedEvent_, 0);

        cudaMemcpyPeerAsync(gpuData_[gpuId+1].maxWindowsInRange_, gpuId+1,
            gpuData_[gpuId].maxWindowsInRange_, gpuId,
            hostData_[hostId].num_queries()*sizeof(window_id),
            gpuData_[gpuId].copyStream_);
    }

    const index_type numBlocks = hostData_[hostId].num_queries();

    //TODO different max cand cases
    if(maxCandidatesPerQuery_ <= 2) {
        constexpr int maxCandidates = 2;

        bool update = gpuId > 0;

        generate_top_candidates<maxCandidates><<<numBlocks,32,0,gpuData_[gpuId].workStream_>>>(
            hostData_[hostId].num_queries(),
            gpuData_[gpuId].resultOffsets_,
            gpuData_[gpuId].queryResults_,
            gpuData_[gpuId].maxWindowsInRange_,
            lineages,
            lowestRank,
            maxCandidatesPerQuery_,
            gpuData_[gpuId].topCandidates_,
            update);

        // cudaStreamSynchronize(gpuData_[gpuId].workStream_);
        // CUERR
    }
    else {
        std::cerr << "At most 2 candidates per query allowed!\n";
    }

    cudaEventRecord(gpuData_[gpuId].tophitsReadyEvent_, gpuData_[gpuId].workStream_);
    cudaStreamWaitEvent(gpuData_[gpuId].copyStream_, gpuData_[gpuId].tophitsReadyEvent_, 0);

    if(gpuId == numGPUs_-1) {
        // copy candidates from last device to host
        cudaMemcpyAsync(hostData_[hostId].top_candidates(), gpuData_[gpuId].topCandidates_,
                        hostData_[hostId].num_queries()*maxCandidatesPerQuery_*sizeof(match_candidate),
                        cudaMemcpyDeviceToHost, gpuData_[gpuId].copyStream_);

        cudaEventRecord(hostData_[hostId].results_copied_event(), gpuData_[gpuId].copyStream_);
    }
    else {
        // copy candidates to next device
        cudaMemcpyPeerAsync(gpuData_[gpuId+1].topCandidates_, gpuId+1,
                            gpuData_[gpuId].topCandidates_, gpuId,
                            hostData_[hostId].num_queries()*maxCandidatesPerQuery_*sizeof(match_candidate),
                            gpuData_[gpuId].copyStream_);
    }

    cudaEventRecord(gpuData_[gpuId].tophitsCopiedEvent_, gpuData_[gpuId].copyStream_);

    // cudaStreamSynchronize(gpuData_[gpuId].copyStream_);
    // CUERR
}


//---------------------------------------------------------------
template<class Location>
void query_batch<Location>::wait_for_allhits_copied(part_id gpuId)
{
    cudaStreamWaitEvent(gpuData_[gpuId].workStream_, gpuData_[gpuId].allhitsCopiedEvent_, 0);
}


//---------------------------------------------------------------
template<class Location>
void query_batch<Location>::sync_work_stream(part_id gpuId)
{
    cudaStreamSynchronize(gpuData_[gpuId].workStream_);
}

//---------------------------------------------------------------
template<class Location>
void query_batch<Location>::sync_copy_stream(part_id gpuId)
{
    cudaStreamSynchronize(gpuData_[gpuId].copyStream_);
}


//---------------------------------------------------------------
template class query_batch<location>;

} // namespace mc

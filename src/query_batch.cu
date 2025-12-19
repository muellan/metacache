/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2022 Robin Kobus  (github.com/funatiq)
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


#include "database.hpp"
#include "gpu_result_processing.cuh"
#include "query_batch.cuh"

#include "../dep/bb_segsort/src/bb_segsort_keys.cuh"


namespace mc {


//---------------------------------------------------------------
template <class Location>
class query_batch<Location>::segmented_sort
{
    using location_type_equivalent = uint64_t;

    static_assert(sizeof(location_type) == sizeof(location_type_equivalent), "location_type must be 64 bit");

public:
    segmented_sort(
        location_type_equivalent *d_keys, location_type_equivalent *d_keysB,
        const int *d_seg_begins, const int *d_seg_ends,
        int *d_binnedSegIds,
        int *d_segBinCounters,
        cudaStream_t stream)
    :
        sorter_{d_keys, d_keysB,
            d_seg_begins, d_seg_ends,
            d_binnedSegIds, d_segBinCounters,
            stream}
    {}

    void run(int numSegs, int maxSegmentSize, cudaStream_t stream) const {
        sorter_.run(numSegs, maxSegmentSize, stream);
    }

private:
    bb_segsort_keys<location_type_equivalent, int> sorter_;
};


//---------------------------------------------------------------
template <class Location>
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
    maxCandidatesPerQuery_{maxCandidatesPerQuery},
    copyAllHits_{copyAllHits}
{
    cudaMallocHost(&queryIds_, maxQueries*sizeof(index_type));
    cudaMallocHost(&sequenceOffsets_, (maxQueries+1)*sizeof(size_type));
    sequenceOffsets_[0] = 0;
    cudaMallocHost(&sequences_, maxSequenceLength*sizeof(char));
    cudaMallocHost(&maxWindowsInRange_, maxQueries*sizeof(window_id));
    CUERR
    if (copyAllHits_)
        cudaMallocHost(&queryResults_, maxQueries*maxResultsPerWindow*sizeof(location_type));
    else
        queryResults_ = nullptr;
    cudaMallocHost(&resultBeginOffsets_, maxQueries*sizeof(int));
    cudaMallocHost(&resultEndOffsets_, maxQueries*sizeof(int));
    cudaMallocHost(&topCandidates_, maxQueries*maxCandidatesPerQuery_*sizeof(match_candidate));
    CUERR

    cudaEventCreate(&resultsCopiedEvent_);
    CUERR
}
//---------------------------------------------------------------
template <class Location>
query_batch<Location>::query_host_data::~query_host_data()
{
    if (queryIds_)           cudaFreeHost(queryIds_);
    if (sequenceOffsets_)    cudaFreeHost(sequenceOffsets_);
    if (sequences_)          cudaFreeHost(sequences_);
    if (maxWindowsInRange_)  cudaFreeHost(maxWindowsInRange_);
    CUERR
    if (queryResults_)       cudaFreeHost(queryResults_);
    if (resultBeginOffsets_) cudaFreeHost(resultBeginOffsets_);
    if (resultEndOffsets_)   cudaFreeHost(resultEndOffsets_);
    if (topCandidates_)      cudaFreeHost(topCandidates_);
    CUERR

    if (resultsCopiedEvent_) cudaEventDestroy(resultsCopiedEvent_);
    CUERR
}
//---------------------------------------------------------------
template <class Location>
query_batch<Location>::query_host_data::query_host_data(query_host_data&& other)
{
    numQueries_ = other.numQueries_;
    numWindows_  = other.numWindows_;
    largestSegmentSize_ = other.largestSegmentSize_;
    maxCandidatesPerQuery_  = other.maxCandidatesPerQuery_;
    copyAllHits_ = other.copyAllHits_;

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
    resultBeginOffsets_       = other.resultBeginOffsets_;
    other.resultBeginOffsets_ = nullptr;
    resultEndOffsets_       = other.resultEndOffsets_;
    other.resultEndOffsets_ = nullptr;
    topCandidates_       = other.topCandidates_;
    other.topCandidates_ = nullptr;

    resultsCopiedEvent_ = other.resultsCopiedEvent_;
    other.resultsCopiedEvent_ = 0;
}
//---------------------------------------------------------------
template <class Location>
void query_batch<Location>::query_host_data::wait_for_results()
{
    cudaEventSynchronize(resultsCopiedEvent_);
    CUERR
}



//---------------------------------------------------------------
template <class Location>
query_batch<Location>::query_gpu_data::query_gpu_data(
    index_type maxQueries,
    size_type maxSequenceLength,
    size_type maxSketchSize,
    size_type maxResultsPerWindow,
    size_type maxCandidatesPerQuery,
    const std::vector<part_id>& gpus,
    part_id gpuId
) :
    offsetSwitcher_{false},
    sorters_{}
{
    size_t allocatedGpuMem = 0;

    cudaSetDevice(gpus[gpuId]); CUERR

    cudaMalloc    (&queryIds_, maxQueries*sizeof(index_type));
    allocatedGpuMem += maxQueries*sizeof(index_type);
    if (gpuId == 0) {
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

    if (gpus.size() > 1) {
        // multi gpu case use sketches buffer
        cudaMalloc    (&sketches_, maxQueries*maxSketchSize*sizeof(feature_type));
        allocatedGpuMem += maxQueries*maxSketchSize*sizeof(feature_type);
    }
    else {
        sketches_ = nullptr;
    }
    cudaMalloc    (&queryResults_, maxQueries*maxResultsPerWindow*sizeof(location_type));
    cudaMalloc    (&queryResultsSorted_, maxQueries*maxResultsPerWindow*sizeof(location_type));
    allocatedGpuMem += 2*maxQueries*maxResultsPerWindow*sizeof(location_type);
    cudaMalloc    (&resultBeginOffsets_[0], maxQueries*sizeof(int));
    cudaMalloc    (&resultBeginOffsets_[1], maxQueries*sizeof(int));
    cudaMalloc    (&resultEndOffsets_[0], maxQueries*sizeof(int));
    cudaMalloc    (&resultEndOffsets_[1], maxQueries*sizeof(int));
    allocatedGpuMem += 4*maxQueries*sizeof(int);
    cudaMalloc    (&binnedSegIds_, maxQueries*sizeof(int));
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
    cudaEventCreate(&d2dCopyEvent_);
    cudaEventCreate(&offsetsCopiedEvent_);
    cudaEventCreate(&allhitsReadyEvent_);
    cudaEventCreate(&allhitsCopiedEvent_);
    cudaEventCreate(&tophitsReadyEvent_);
    cudaEventCreate(&tophitsCopiedEvent_);
    CUERR

    using location_type_equivalent = uint64_t;
    static_assert(sizeof(location_type) == sizeof(location_type_equivalent), "location_type must be 64 bit");

    for (int i = 0; i < 2; ++i) {
        sorters_.emplace_back(
            (location_type_equivalent*)(queryResults_),
            (location_type_equivalent*)(queryResultsSorted_),
            resultBeginOffsets_[i],
            resultEndOffsets_[i],
            binnedSegIds_,
            segBinCounters_,
            workStream_);
        CUERR
    }
}
//---------------------------------------------------------------
template <class Location>
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
    queryResultsSorted_ = other.queryResultsSorted_;
    other.queryResultsSorted_ = nullptr;
    resultBeginOffsets_[0] = other.resultBeginOffsets_[0];
    other.resultBeginOffsets_[0] = nullptr;
    resultBeginOffsets_[1] = other.resultBeginOffsets_[1];
    other.resultBeginOffsets_[1] = nullptr;
    resultEndOffsets_[0] = other.resultEndOffsets_[0];
    other.resultEndOffsets_[0] = nullptr;
    resultEndOffsets_[1] = other.resultEndOffsets_[1];
    other.resultEndOffsets_[1] = nullptr;
    binnedSegIds_ = other.binnedSegIds_;
    other.binnedSegIds_ = nullptr;
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
    d2dCopyEvent_ = other.d2dCopyEvent_;
    other.d2dCopyEvent_ = 0;
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

    sorters_ = std::move(other.sorters_);

    offsetSwitcher_ = other.offsetSwitcher_;
    other.offsetSwitcher_ = false;
}
//---------------------------------------------------------------
template <class Location>
query_batch<Location>::query_gpu_data::~query_gpu_data()
{
    if (queryIds_)           cudaFree(queryIds_);
    if (sequenceOffsets_)    cudaFree(sequenceOffsets_);
    if (sequences_)          cudaFree(sequences_);
    if (sketches_)           cudaFree(sketches_);
    if (queryResults_)       cudaFree(queryResults_);
    if (queryResultsSorted_) cudaFree(queryResultsSorted_);
    if (resultBeginOffsets_[0]) cudaFree(resultBeginOffsets_[0]);
    if (resultBeginOffsets_[1]) cudaFree(resultBeginOffsets_[1]);
    if (resultEndOffsets_[0])   cudaFree(resultEndOffsets_[0]);
    if (resultEndOffsets_[1])   cudaFree(resultEndOffsets_[1]);
    if (binnedSegIds_)       cudaFree(binnedSegIds_);
    if (segBinCounters_)     cudaFree(segBinCounters_);
    if (topCandidates_)      cudaFree(topCandidates_);
    if (maxWindowsInRange_)  cudaFree(maxWindowsInRange_);
    CUERR
    if (workStream_) cudaStreamDestroy(workStream_);
    if (copyStream_) cudaStreamDestroy(copyStream_);
    CUERR
    if (queryFinishedEvent_)    cudaEventDestroy(queryFinishedEvent_);
    if (d2dCopyEvent_)          cudaEventDestroy(d2dCopyEvent_);
    if (offsetsCopiedEvent_)    cudaEventDestroy(offsetsCopiedEvent_);
    if (allhitsReadyEvent_)     cudaEventDestroy(allhitsReadyEvent_);
    if (allhitsCopiedEvent_)    cudaEventDestroy(allhitsCopiedEvent_);
    if (tophitsReadyEvent_)     cudaEventDestroy(tophitsReadyEvent_);
    if (tophitsCopiedEvent_)    cudaEventDestroy(tophitsCopiedEvent_);
    CUERR
}


//---------------------------------------------------------------
template <class Location>
query_batch<Location>::query_batch(
    index_type maxQueries,
    size_type maxSequenceLength,
    size_type maxSketchSize,
    size_type maxResultsPerWindow,
    size_type maxCandidatesPerQuery,
    bool copyAllHits,
    part_id numHostThreads,
    part_id numGPUs,
    unsigned replica
) :
    numGPUs_{numGPUs},
    maxWindows_{maxQueries},
    maxSequenceLength_{maxSequenceLength},
    maxSketchSize_{maxSketchSize},
    maxResultsPerWindow_{maxResultsPerWindow},
    maxCandidatesPerQuery_{maxCandidatesPerQuery},
    gpus_(numGPUs),
    hostData_{},
    gpuData_{}
{
    cudaStreamCreate(&h2dCopyStream_);

    cudaEventCreate(&queriesCopiedEvent_);
    cudaEventCreate(&maxWinCopiedEvent_);
    CUERR

    for (part_id gpuId = 0; gpuId < numGPUs_; ++gpuId)
        gpus_[gpuId] = replica*numGPUs_ + gpuId;

    cudaSetDevice(gpus_.back()); CUERR

    for (part_id hostId = 0; hostId < numHostThreads; ++hostId) {
        hostData_.emplace_back(maxQueries, maxSequenceLength, maxResultsPerWindow, maxCandidatesPerQuery, copyAllHits);
    }

    for (part_id gpuId = 0; gpuId < numGPUs_; ++gpuId) {
        gpuData_.emplace_back(maxQueries, maxSequenceLength, maxSketchSize, maxResultsPerWindow, maxCandidatesPerQuery,
                              gpus_, gpuId);
    }
}
//---------------------------------------------------------------
template <class Location>
query_batch<Location>::query_batch(query_batch&& other) {
    numGPUs_ = other.numGPUs_;

    maxWindows_ = other.maxWindows_;
    maxSequenceLength_ = other.maxSequenceLength_;
    maxSketchSize_ = other.maxSketchSize_;
    maxResultsPerWindow_ = other.maxResultsPerWindow_;
    maxCandidatesPerQuery_ = other.maxCandidatesPerQuery_;

    gpus_ = std::move(other.gpus_);
    hostData_ = std::move(other.hostData_);
    gpuData_ = std::move(other.gpuData_);

    h2dCopyStream_ = other.h2dCopyStream_;
    other.h2dCopyStream_ = 0;

    queriesCopiedEvent_ = other.queriesCopiedEvent_;
    other.queriesCopiedEvent_ = 0;
    maxWinCopiedEvent_ = other.maxWinCopiedEvent_;
    other.maxWinCopiedEvent_ = 0;
};

//---------------------------------------------------------------
template <class Location>
query_batch<Location>::~query_batch()
{
    if (h2dCopyStream_) cudaStreamDestroy(h2dCopyStream_);
    CUERR
    if (queriesCopiedEvent_) cudaEventDestroy(queriesCopiedEvent_);
    if (maxWinCopiedEvent_) cudaEventDestroy(maxWinCopiedEvent_);
    CUERR
}


//---------------------------------------------------------------
template <class Location>
void query_batch<Location>::copy_queries_to_device_async(part_id hostId)
{
    // copy from host to device 0
    part_id gpuId = 0;

    cudaSetDevice(gpus_[gpuId]);

    // wait for previous batch
    cudaStreamWaitEvent(h2dCopyStream_, gpuData_[gpuId].d2dCopyEvent_, 0);

    cudaMemcpyAsync(gpuData_[gpuId].queryIds_, hostData_[hostId].query_ids(),
                    hostData_[hostId].num_windows()*sizeof(index_type),
                    cudaMemcpyHostToDevice, h2dCopyStream_);
    cudaMemcpyAsync(gpuData_[gpuId].sequenceOffsets_, hostData_[hostId].sequence_offsets(),
                    (hostData_[hostId].num_windows()+1)*sizeof(size_type),
                    cudaMemcpyHostToDevice, h2dCopyStream_);
    cudaMemcpyAsync(gpuData_[gpuId].sequences_, hostData_[hostId].sequences(),
                    hostData_[hostId].sequence_offsets()[hostData_[hostId].num_windows()]*sizeof(char),
                    cudaMemcpyHostToDevice, h2dCopyStream_);
    cudaMemcpyAsync(gpuData_[gpuId].result_begin_offsets(), hostData_[hostId].result_begin_offsets(),
                    hostData_[hostId].num_queries()*sizeof(int),
                    cudaMemcpyHostToDevice, h2dCopyStream_);
    cudaMemcpyAsync(gpuData_[gpuId].result_end_offsets(), gpuData_[gpuId].result_begin_offsets(),
                    hostData_[hostId].num_queries()*sizeof(int),
                    cudaMemcpyDeviceToDevice, h2dCopyStream_);

    cudaEventRecord(queriesCopiedEvent_, h2dCopyStream_);
    // work has to wait until copies finished
    cudaStreamWaitEvent(gpuData_[gpuId].workStream_, queriesCopiedEvent_, 0);

    // cudaEventSynchronize(queriesCopiedEvent_);
    // CUERR

    // wait for previous batch
    cudaStreamWaitEvent(h2dCopyStream_, gpuData_[gpuId].tophitsCopiedEvent_, 0);

    cudaMemcpyAsync(gpuData_[gpuId].maxWindowsInRange_, hostData_[hostId].max_windows_in_range(),
                    hostData_[hostId].num_queries()*sizeof(window_id),
                    cudaMemcpyHostToDevice, h2dCopyStream_);

    cudaEventRecord(maxWinCopiedEvent_, h2dCopyStream_);

    // cudaEventSynchronize(maxWinCopiedEvent_);
    // CUERR
}


//---------------------------------------------------------------
template <class Location>
void query_batch<Location>::copy_queries_to_next_device_async(part_id hostId, part_id gpuId)
{
    // copy from device gpuId to device gpuId+1
    cudaSetDevice(gpus_[gpuId]);

    // wait for previous batch
    cudaStreamWaitEvent(gpuData_[gpuId].copyStream_, gpuData_[gpuId+1].queryFinishedEvent_, 0);
    cudaStreamWaitEvent(gpuData_[gpuId].copyStream_, gpuData_[gpuId+1].d2dCopyEvent_, 0);

    // wait until data available
    if (gpuId == 0)
        cudaStreamWaitEvent(gpuData_[gpuId].copyStream_, queriesCopiedEvent_, 0);
    else
        cudaStreamWaitEvent(gpuData_[gpuId].copyStream_, gpuData_[gpuId-1].d2dCopyEvent_, 0);

    cudaMemcpyPeerAsync(gpuData_[gpuId+1].queryIds_, gpus_[gpuId+1],
                        gpuData_[gpuId].queryIds_, gpus_[gpuId],
                        hostData_[hostId].num_windows()*sizeof(index_type),
                        gpuData_[gpuId].copyStream_);

    cudaMemcpyPeerAsync(gpuData_[gpuId+1].result_begin_offsets(), gpus_[gpuId+1],
                        gpuData_[gpuId].result_begin_offsets(), gpus_[gpuId],
                        hostData_[hostId].num_windows()*sizeof(index_type),
                        gpuData_[gpuId].copyStream_);

    // wait until query produced sketches
    if (gpuId == 0)
        cudaStreamWaitEvent(gpuData_[gpuId].copyStream_, gpuData_[gpuId].queryFinishedEvent_, 0);

    cudaMemcpyPeerAsync(gpuData_[gpuId+1].sketches_, gpus_[gpuId+1],
                        gpuData_[gpuId].sketches_, gpus_[gpuId],
                        hostData_[hostId].num_windows()*maxSketchSize_*sizeof(feature_type),
                        gpuData_[gpuId].copyStream_);

    cudaEventRecord(gpuData_[gpuId].d2dCopyEvent_, gpuData_[gpuId].copyStream_);

    // cudaEventSynchronize(gpuData_[gpuId].d2dCopyEvent_);
    // CUERR

    // wait until data has arrived
    if (gpuId == 0)
        cudaStreamWaitEvent(gpuData_[gpuId].copyStream_, maxWinCopiedEvent_, 0);
    else
        cudaStreamWaitEvent(gpuData_[gpuId].copyStream_, gpuData_[gpuId-1].tophitsCopiedEvent_, 0);

    // wait until previous batch finished
    cudaStreamWaitEvent(gpuData_[gpuId].copyStream_, gpuData_[gpuId+1].tophitsCopiedEvent_, 0);

    cudaMemcpyPeerAsync(gpuData_[gpuId+1].maxWindowsInRange_, gpus_[gpuId+1],
                        gpuData_[gpuId].maxWindowsInRange_, gpus_[gpuId],
                        hostData_[hostId].num_queries()*sizeof(window_id),
                        gpuData_[gpuId].copyStream_);


    cudaSetDevice(gpus_[gpuId+1]);

    // work has to wait until copies finished
    cudaStreamWaitEvent(gpuData_[gpuId+1].workStream_, gpuData_[gpuId].d2dCopyEvent_, 0);

    cudaMemcpyAsync(gpuData_[gpuId+1].result_end_offsets(), gpuData_[gpuId+1].result_begin_offsets(),
                    hostData_[hostId].num_queries()*sizeof(int),
                    cudaMemcpyDeviceToDevice, gpuData_[gpuId+1].workStream_);
}


//---------------------------------------------------------------
template <class Location>
void query_batch<Location>::mark_query_finished(part_id gpuId)
{
    cudaEventRecord(gpuData_[gpuId].queryFinishedEvent_, gpuData_[gpuId].workStream_);

    // cudaStreamSynchronize(gpuData_[gpuId].workStream_);
    // CUERR
}


//---------------------------------------------------------------
template <class Location>
void query_batch<Location>::sort_and_copy_allhits_async(part_id hostId, part_id gpuId)
{
    cudaSetDevice(gpus_[gpuId]);

    if (hostData_[hostId].allhits_requested()) {
        cudaStreamWaitEvent(gpuData_[gpuId].copyStream_, gpuData_[gpuId].queryFinishedEvent_, 0);

        cudaMemcpyAsync(hostData_[hostId].result_end_offsets(), gpuData_[gpuId].result_end_offsets(),
                        hostData_[hostId].num_queries()*sizeof(int),
                        cudaMemcpyDeviceToHost, gpuData_[gpuId].copyStream_);

        cudaEventRecord(gpuData_[gpuId].offsetsCopiedEvent_, gpuData_[gpuId].copyStream_);

        // cudaStreamSynchronize(gpuData_[gpuId].copyStream_);
        // CUERR

        // wait until previous batch finished copying
        cudaStreamWaitEvent(gpuData_[gpuId].workStream_, gpuData_[gpuId].allhitsCopiedEvent_, 0);
    }

    gpuData_[gpuId].sort(
        hostData_[hostId].num_queries(), hostData_[hostId].largest_segment_size(), gpuData_[gpuId].workStream_);
    // cudaStreamSynchronize(gpuData_[gpuId].workStream_);
    // CUERR

    if (hostData_[hostId].allhits_requested()) {
        cudaEventRecord(gpuData_[gpuId].allhitsReadyEvent_, gpuData_[gpuId].workStream_);
        cudaStreamWaitEvent(gpuData_[gpuId].copyStream_, gpuData_[gpuId].allhitsReadyEvent_, 0);

        cudaEventSynchronize(gpuData_[gpuId].offsetsCopiedEvent_);

        for (index_type i = 0; i < hostData_[hostId].num_queries(); ++i) {
            int length = hostData_[hostId].result_end_offsets()[i] - hostData_[hostId].result_begin_offsets()[i];

            cudaMemcpyAsync(
                hostData_[hostId].query_results()+hostData_[hostId].result_begin_offsets()[i],
                gpuData_[gpuId].queryResultsSorted_+hostData_[hostId].result_begin_offsets()[i],
                length*sizeof(location_type),
                cudaMemcpyDeviceToHost, gpuData_[gpuId].copyStream_);
        }

        cudaEventRecord(gpuData_[gpuId].allhitsCopiedEvent_, gpuData_[gpuId].copyStream_);

        // cudaStreamSynchronize(gpuData_[gpuId].copyStream_);
        // CUERR
    }
}


//---------------------------------------------------------------
template <class Location>
void query_batch<Location>::generate_and_copy_top_candidates_async(
    part_id hostId,
    part_id gpuId,
    const ranked_lineage * lineages,
    taxon_rank lowestRank)
{
    cudaSetDevice(gpus_[gpuId]);

    // wait until data has arrived
    if (gpuId == 0)
        cudaStreamWaitEvent(gpuData_[gpuId].workStream_, maxWinCopiedEvent_, 0);
    else
        cudaStreamWaitEvent(gpuData_[gpuId].workStream_, gpuData_[gpuId-1].tophitsCopiedEvent_, 0);

    const index_type numBlocks = hostData_[hostId].num_queries();

    // TODO different max cand cases
    if (maxCandidatesPerQuery_ <= 2) {
        constexpr int maxCandidates = 2;

        bool update = gpuId > 0;

        generate_top_candidates<maxCandidates><<<numBlocks,32,0,gpuData_[gpuId].workStream_>>>(
            hostData_[hostId].num_queries(),
            gpuData_[gpuId].result_begin_offsets(),
            gpuData_[gpuId].result_end_offsets(),
            gpuData_[gpuId].queryResultsSorted_,
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

    if (gpuId == numGPUs_-1) {
        // copy candidates from last device to host
        cudaMemcpyAsync(hostData_[hostId].top_candidates(), gpuData_[gpuId].topCandidates_,
                        hostData_[hostId].num_queries()*maxCandidatesPerQuery_*sizeof(match_candidate),
                        cudaMemcpyDeviceToHost, gpuData_[gpuId].copyStream_);

        cudaEventRecord(hostData_[hostId].results_copied_event(), gpuData_[gpuId].copyStream_);
    }
    else {
        // copy candidates to next device
        cudaMemcpyPeerAsync(gpuData_[gpuId+1].topCandidates_, gpus_[gpuId+1],
                            gpuData_[gpuId].topCandidates_, gpus_[gpuId],
                            hostData_[hostId].num_queries()*maxCandidatesPerQuery_*sizeof(match_candidate),
                            gpuData_[gpuId].copyStream_);
    }

    cudaEventRecord(gpuData_[gpuId].tophitsCopiedEvent_, gpuData_[gpuId].copyStream_);

    // cudaStreamSynchronize(gpuData_[gpuId].copyStream_);
    // CUERR
}


//---------------------------------------------------------------
template <class Location>
void query_batch<Location>::sync_work_stream(part_id gpuId)
{
    cudaStreamSynchronize(gpuData_[gpuId].workStream_);
}

//---------------------------------------------------------------
template <class Location>
void query_batch<Location>::sync_copy_stream(part_id gpuId)
{
    cudaStreamSynchronize(gpuData_[gpuId].copyStream_);
}


//---------------------------------------------------------------
template class query_batch<location>;

} // namespace mc

/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2021 Robin Kobus  (kobus@uni-mainz.de)
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

#ifndef MC_QUERY_BATCH_H_
#define MC_QUERY_BATCH_H_


#include "candidate_structs.h"
#include "config.h"
#include "hash_dna.h"
#include "taxonomy.h"

#include "cuda_runtime.h"

#include <functional>
#include <memory>


namespace mc {


/*************************************************************************//**
 *
 * @brief batch contains sequence data & query results of multiple reads,
 *        manages allocated memory on host and device,
 *        moves data between host & device,
 *        uses its own stream
 *
 *****************************************************************************/
template<class Location>
class query_batch
{
    class segmented_sort;

public:
    using index_type      = uint32_t;
    using size_type       = uint32_t;
    using location_type   = Location;
    using match_locations = typename std::vector<location_type>;
    using feature_type    = typename sketcher::feature_type;

    using taxon_rank     = taxonomy::rank;
    using ranked_lineage = taxonomy::ranked_lineage;

    //---------------------------------------------------------------
    class query_host_data
    {
    public:
        query_host_data(index_type maxQueries,
                        size_type maxSequenceLength,
                        size_type maxResultsPerWindow,
                        size_type maxCandidatesPerQuery,
                        bool copyAllHits);
        query_host_data(const query_host_data&) = delete;
        query_host_data(query_host_data&&);
        ~query_host_data();

        /*************************************************************************//**
        *
        * @brief add sequence pair to batch as windows of characters
        *        if sequence pair does not fit into batch, don't add it
        *
        * @details each window of a sequence is added as a separate query
        *
        * @return true if added, false otherwise
        *
        *****************************************************************************/
        template<class InputIterator>
        bool add_paired_read(
            InputIterator first1, InputIterator last1,
            InputIterator first2, InputIterator last2,
            const sketching_opt& querySketching,
            const candidate_generation_rules& rules,
            index_type maxWindows,
            uint32_t maxSequenceLength,
            size_type maxResultsPerWindow
        ) {
            using std::distance;

            const numk_t kmerSize = querySketching.kmerlen;
            const size_t windowSize = querySketching.winlen;
            const size_t windowStride = querySketching.winstride;

            const size_t seqLength1 = distance(first1, last1);
            const size_t seqLength2 = distance(first2, last2);

            // set offset depending on previous number of windows in batch
            resultBeginOffsets_[numQueries_] = numWindows_*maxResultsPerWindow;

            // no kmers in sequence
            if(seqLength1 < kmerSize && seqLength2 < kmerSize) {
                // batch full, nothing processed
                if(numWindows_ + 1 > maxWindows) return false;

                // insert empty query
                queryIds_[numWindows_] = numQueries_;
                sequenceOffsets_[numWindows_+1] = sequenceOffsets_[numWindows_];

                ++numWindows_;
                ++numQueries_;
                return true;
            }

            const window_id windowsInSeq1 = (seqLength1 >= kmerSize) ?
                                            (seqLength1-kmerSize + windowStride) / windowStride : 0;
            const window_id windowsInSeq2 = (seqLength2 >= kmerSize) ?
                                            (seqLength2-kmerSize + windowStride) / windowStride : 0;

            // batch full, nothing processed
            if(numWindows_ + windowsInSeq1 + windowsInSeq2 > maxWindows) return false;

            const auto availableSize = maxSequenceLength - sequenceOffsets_[numWindows_];
            const auto windowSizePadded = (windowSize + 3) / 4 * 4;
            // batch full, nothing processed
            if((windowsInSeq1 + windowsInSeq2)*windowSizePadded > availableSize) return false;

            index_type windowsInQuery = 0;

            // insert first sequence into batch as separate windows
            for_each_window(first1, last1, windowSize, windowStride,
                [&] (InputIterator first, InputIterator last) {
                    auto length = distance(first, last);
                    if(length >= kmerSize) {
                        // store window affiliation
                        queryIds_[numWindows_] = numQueries_;
                        // copy characters and pad for vectorized access
                        std::copy(first, last, sequences_ + sequenceOffsets_[numWindows_]);
                        auto lengthPadded = (length + 3) / 4 * 4;
                        std::fill(sequences_ + sequenceOffsets_[numWindows_] + length,
                                sequences_ + sequenceOffsets_[numWindows_] + lengthPadded,
                                'N');
                        sequenceOffsets_[numWindows_+1] = sequenceOffsets_[numWindows_] + lengthPadded;

                        ++numWindows_;
                        ++windowsInQuery;
                    }
                }
            );

            // insert second sequence into batch as separate windows
            for_each_window(first2, last2, windowSize, windowStride,
                [&] (InputIterator first, InputIterator last) {
                    auto length = distance(first, last);
                    if(length >= kmerSize) {
                        // store window affiliation
                        queryIds_[numWindows_] = numQueries_;
                        // copy characters and pad for vectorized access
                        std::copy(first, last, sequences_ + sequenceOffsets_[numWindows_]);
                        auto lengthPadded = (length + 3) / 4 * 4;
                        std::fill(sequences_ + sequenceOffsets_[numWindows_] + length,
                                sequences_ + sequenceOffsets_[numWindows_] + lengthPadded,
                                'N');
                        sequenceOffsets_[numWindows_+1] = sequenceOffsets_[numWindows_] + lengthPadded;

                        ++numWindows_;
                        ++windowsInQuery;
                    }
                }
            );

            maxWindowsInRange_[numQueries_] = rules.maxWindowsInRange;

            ++numQueries_;

            size_type segmentSize = windowsInQuery*maxResultsPerWindow;
            if(largestSegmentSize_ < segmentSize) largestSegmentSize_ = segmentSize;

            return true;
        }

        //-----------------------------------------------------
        template<class Sequence>
        bool add_paired_read(
            Sequence seq1, Sequence seq2,
            const sketching_opt& querySketching,
            const candidate_generation_rules& rules,
            index_type maxWindows,
            uint32_t maxSequenceLength,
            size_type maxResultsPerWindow
        ) {
            using std::begin;
            using std::end;

            return add_paired_read(
                begin(seq1), end(seq1),
                begin(seq2), end(seq2),
                querySketching,
                rules,
                maxWindows,
                maxSequenceLength,
                maxResultsPerWindow);
        }

        //---------------------------------------------------------------
        span<const location_type> allhits(index_type id) const noexcept {
            if(copyAllHits_&& id < num_queries()) {
                location_type * begin = query_results()+result_begin_offsets()[id];
                location_type * end = query_results()+result_end_offsets()[id];

                return span<const location_type>{begin, size_t(end-begin)};
            }
            else
                return span<const location_type>{};
        }
        //---------------------------------------------------------------
        span<const match_candidate> top_candidates(index_type id) const noexcept {
            if(id < num_queries())
                return span<const match_candidate>{
                    top_candidates()+id*maxCandidatesPerQuery_,
                    maxCandidatesPerQuery_
                };
            else
                return span<const match_candidate>{};
        }

        //---------------------------------------------------------------
        index_type num_queries() const noexcept { return numQueries_; }
        index_type num_windows() const noexcept { return numWindows_; }
        size_type  largest_segment_size() const noexcept { return largestSegmentSize_; }
        bool       allhitsCopyNeeded() const noexcept { return copyAllHits_; }

        index_type * query_ids() const noexcept { return queryIds_; };
        size_type  * sequence_offsets() const noexcept { return sequenceOffsets_; };
        char       * sequences() const noexcept { return sequences_; };
        window_id  * max_windows_in_range() const noexcept { return maxWindowsInRange_; };

        location_type   * query_results() const noexcept { return queryResults_; }
        int             * result_begin_offsets() const noexcept { return resultBeginOffsets_; }
        int             * result_end_offsets() const noexcept { return resultEndOffsets_; }
        match_candidate * top_candidates() const noexcept { return topCandidates_; }

        cudaEvent_t& results_copied_event() noexcept { return resultsCopiedEvent_; }

        //---------------------------------------------------------------
        void wait_for_results();

        //---------------------------------------------------------------
        void clear() noexcept {
            numQueries_ = 0;
            numWindows_ = 0;
            largestSegmentSize_ = 0;
        }

    private:
        index_type numQueries_;
        index_type numWindows_;
        size_type largestSegmentSize_;

        size_type  maxCandidatesPerQuery_;
        bool       copyAllHits_;

        // input
        index_type * queryIds_;
        size_type  * sequenceOffsets_;
        char       * sequences_;
        window_id  * maxWindowsInRange_;
        // output
        location_type   * queryResults_;
        int             * resultBeginOffsets_;
        int             * resultEndOffsets_;
        match_candidate * topCandidates_;

        cudaEvent_t resultsCopiedEvent_;
    };

    //---------------------------------------------------------------
    struct query_gpu_data
    {
        query_gpu_data(index_type maxQueries,
                       size_type maxSequenceLength,
                       size_type maxSketchSize,
                       size_type maxResultsPerWindow,
                       size_type maxCandidatesPerQuery,
                       const std::vector<part_id>& gpus,
                       part_id gpuId);
        query_gpu_data(const query_gpu_data&) = delete;
        query_gpu_data(query_gpu_data&&);
        ~query_gpu_data();

        //---------------------------------------------------------------
        index_type      * queryIds_;
        size_type       * sequenceOffsets_;
        char            * sequences_;
        feature_type    * sketches_;
        window_id       * maxWindowsInRange_;

        location_type   * queryResults_;
        location_type   * queryResultsSorted_;
        int             * resultBeginOffsets_;
        int             * resultEndOffsets_;

        int             * binnedSegIds_;
        int             * segBinCounters_;

        match_candidate * topCandidates_;

        cudaStream_t workStream_;
        cudaStream_t copyStream_;

        cudaEvent_t queryFinishedEvent_;
        cudaEvent_t sketchesCopiedEvent_;
        cudaEvent_t offsetsCopiedEvent_;
        cudaEvent_t allhitsReadyEvent_;
        cudaEvent_t allhitsCopiedEvent_;
        cudaEvent_t tophitsReadyEvent_;
        cudaEvent_t tophitsCopiedEvent_;
    };


    //---------------------------------------------------------------
    /** @brief allocate memory on host and device */
    query_batch(index_type maxQueries,
                size_type maxEncodeLength,
                size_type maxSketchSize,
                size_type maxResultsPerWindow,
                size_type maxCandidatesPerQuery,
                bool copyAllHits,
                part_id numHostThreads,
                part_id numGPUs,
                unsigned replica);
    query_batch(const query_batch&) = delete;
    query_batch(query_batch&&);
    /** @brief free memory allocation */
    ~query_batch();

    //---------------------------------------------------------------
    part_id num_gpus() const noexcept {
        return numGPUs_;
    }
    //---------------------------------------------------------------
    part_id gpu(part_id gpuId) const noexcept {
        return gpus_[gpuId];
    }
    //---------------------------------------------------------------
    query_host_data& host_data(part_id hostId) noexcept {
        return hostData_[hostId];
    }
    //---------------------------------------------------------------
    const query_gpu_data& gpu_data(part_id gpuId) const noexcept {
        return gpuData_[gpuId];
    }

    //---------------------------------------------------------------
    void sync_work_stream(part_id part_id);
    void sync_copy_stream(part_id part_id);

    //---------------------------------------------------------------
    /** @brief add read to host batch */
    template<class... Args>
    bool add_paired_read(part_id hostId, Args&&... args)
    {
        return hostData_[hostId].add_paired_read(
            std::forward<Args>(args)...,
            maxWindows_,
            maxSequenceLength_,
            maxResultsPerWindow_);
    }

    //---------------------------------------------------------------
    /** @brief asynchronously copy queries to device */
    void copy_queries_to_device_async(part_id hostId);
    void copy_queries_to_next_device_async(part_id hostId, part_id gpuId);
    //-----------------------------------------------------
    /** @brief synchronize event after copy to device 0 */
    void wait_for_queries_copied();
    //---------------------------------------------------------------
    /** @brief record event after sketch creation on device 0 */
    void mark_query_finished(part_id gpuId);
    //---------------------------------------------------------------
    /** @brief make work stream wait for allhits copy (before starting new query) */
    void wait_for_allhits_copied(part_id gpuId);

    //---------------------------------------------------------------
    /**
     * @brief asynchronously sort in work stream,
     *        and copy allhits to host in copy stream if needed
     */
    void sort_and_copy_allhits_async(part_id hostId, part_id gpuId);
    /**
     * @brief asynchronously generate top candidates in work stream
     *        and copy top candidates in copy stream
     */
    void generate_and_copy_top_candidates_async(
        part_id hostId,
        part_id gpuId,
        const ranked_lineage * lineages,
        taxon_rank lowestRank);


    //---------------------------------------------------------------
private:
    part_id numGPUs_;

    index_type maxWindows_;
    size_type  maxSequenceLength_;
    size_type  maxSketchSize_;
    size_type  maxResultsPerWindow_;
    size_type  maxCandidatesPerQuery_;

    std::vector<part_id> gpus_;
    std::vector<query_host_data> hostData_;
    std::vector<query_gpu_data> gpuData_;
    std::vector<segmented_sort> sorters_;

    cudaStream_t h2dCopyStream_;

    cudaEvent_t queriesCopiedEvent_;
    cudaEvent_t maxWinCopiedEvent_;
};


} // namespace mc


#endif
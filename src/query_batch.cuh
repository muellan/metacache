#ifndef MC_QUERY_BATCH_H_
#define MC_QUERY_BATCH_H_

#include <functional>
#include <memory>

#include "cuda_runtime.h"

#include "config.h"
#include "candidate_generation.h"
#include "hash_dna.h"
#include "taxonomy.h"

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
    using index_type     = uint32_t;
    using size_type      = uint32_t;
    using location_type  = Location;
    using feature_type   = typename sketcher::feature_type;

    using taxon_rank     = taxonomy::rank;
    using ranked_lineage = taxonomy::ranked_lineage;

    //---------------------------------------------------------------
    class query_host_data
    {
    public:
        query_host_data(index_type maxQueries,
                        size_type maxSequenceLength,
                        size_type maxResultsPerQuery,
                        size_type maxCandidatesPerQuery);
        query_host_data(const query_host_data&) = delete;
        query_host_data(query_host_data&&);
        ~query_host_data();

        /*************************************************************************//**
        *
        * @brief add sequence pair to batch as windows of encoded characters
        *        if sequence pair does not fit into batch, don't add it
        *
        * @detail each window of a sequence is added as a separate query
        *
        * @return true if added, false otherwise
        *
        *****************************************************************************/
        template<class InputIterator>
        bool add_paired_read(
            InputIterator first1, InputIterator last1,
            InputIterator first2, InputIterator last2,
            const sketcher& querySketcher,
            size_t insertSizeMax,
            index_type maxQueries,
            uint32_t maxSequenceLength
        ) {
            using std::distance;

            const numk_t kmerSize = querySketcher.kmer_size();
            const size_t windowSize = querySketcher.window_size();
            const size_t windowStride = querySketcher.window_stride();

            const size_t seqLength1 = distance(first1, last1);
            const size_t seqLength2 = distance(first2, last2);

            // no kmers in sequence
            if(seqLength1 < kmerSize && seqLength2 < kmerSize) {
                // batch full, nothing processed
                if(numQueries_ + 1 > maxQueries) return false;

                // insert empty query
                queryIds_[numQueries_] = numSegments_;
                sequenceOffsets_[numQueries_+1] = sequenceOffsets_[numQueries_];

                ++numQueries_;
                ++numSegments_;
                return true;
            }

            const window_id numWindows1 = (seqLength1-kmerSize + windowStride) / windowStride;
            const window_id numWindows2 = (seqLength2-kmerSize + windowStride) / windowStride;

            // batch full, nothing processed
            if(numQueries_ + numWindows1 + numWindows2 > maxQueries) return false;

            const auto availableSize = maxSequenceLength - sequenceOffsets_[numQueries_];
            const auto windowSizePadded = (windowSize + 3) / 4 * 4;
            // batch full, nothing processed
            if((numWindows1 + numWindows2)*windowSizePadded > availableSize) return false;

            // insert first sequence into batch as separate windows
            for_each_window(first1, last1, windowSize, windowStride,
                [&] (InputIterator first, InputIterator last) {
                    auto length = distance(first, last);
                    if(length >= kmerSize) {
                        queryIds_[numQueries_] = numSegments_;
                        std::copy(first, last, sequences_ + sequenceOffsets_[numQueries_]);
                        auto lengthPadded = (length + 3) / 4 * 4;
                        std::fill(sequences_ + sequenceOffsets_[numQueries_] + length,
                                sequences_ + sequenceOffsets_[numQueries_] + lengthPadded,
                                'N');
                        sequenceOffsets_[numQueries_+1] = sequenceOffsets_[numQueries_] + lengthPadded;

                        ++numQueries_;
                    }
                }
            );

            // insert second sequence into batch as separate windows
            for_each_window(first2, last2, windowSize, windowStride,
                [&] (InputIterator first, InputIterator last) {
                    auto length = distance(first, last);
                    if(length >= kmerSize) {
                        queryIds_[numQueries_] = numSegments_;
                        std::copy(first, last, sequences_ + sequenceOffsets_[numQueries_]);
                        auto lengthPadded = (length + 3) / 4 * 4;
                        std::fill(sequences_ + sequenceOffsets_[numQueries_] + length,
                                sequences_ + sequenceOffsets_[numQueries_] + lengthPadded,
                                'N');
                        sequenceOffsets_[numQueries_+1] = sequenceOffsets_[numQueries_] + lengthPadded;

                        ++numQueries_;
                    }
                }
            );

            maxWindowsInRange_[numSegments_] = window_id( 2 +
                (std::max(seqLength1 + seqLength2, insertSizeMax) / windowStride ));

            ++numSegments_;

            return true;
        }

        //-----------------------------------------------------
        template<class Sequence>
        bool add_paired_read(
            Sequence seq1, Sequence seq2,
            const sketcher& querySketcher,
            size_t insertSizeMax,
            index_type maxQueries,
            uint32_t maxSequenceLength
        ) {
            using std::begin;
            using std::end;

            return add_paired_read(
                begin(seq1), end(seq1),
                begin(seq2), end(seq2),
                querySketcher,
                insertSizeMax,
                maxQueries,
                maxSequenceLength);
        }

        //---------------------------------------------------------------
        span<location_type> allhits(index_type id) const noexcept {
            if(id < num_segments())
                return span<location_type>{
                    query_results()+result_offsets()[id],
                    query_results()+result_offsets()[id+1]
                };
            else
                return span<location_type>{};
        }
        //---------------------------------------------------------------
        span<match_candidate> top_candidates(index_type id) const noexcept {
            if(id < num_segments())
                return span<match_candidate>{
                    top_candidates()+id*maxCandidatesPerQuery_,
                    top_candidates()+(id+1)*maxCandidatesPerQuery_
                };
            else
                return span<match_candidate>{};
        }

        //---------------------------------------------------------------
        index_type num_segments() const noexcept { return numSegments_; }
        index_type num_queries() const noexcept { return numQueries_; }

        index_type * query_ids() const noexcept { return queryIds_; };
        size_type  * sequence_offsets() const noexcept { return sequenceOffsets_; };
        char       * sequences() const noexcept { return sequences_; };
        window_id  * max_windows_in_range() const noexcept { return maxWindowsInRange_; };

        location_type   * query_results() const noexcept { return queryResults_; }
        int             * result_offsets() const noexcept { return resultOffsets_; }
        match_candidate * top_candidates() const noexcept { return topCandidates_; }

        cudaEvent_t& results_copied_event() noexcept { return resultsCopiedEvent_; }

        //---------------------------------------------------------------
        void wait_for_results();

        //---------------------------------------------------------------
        void clear() noexcept {
            numSegments_ = 0;
            numQueries_ = 0;
        }

    private:
        index_type numSegments_;
        index_type numQueries_;

        size_type  maxCandidatesPerQuery_;

        // input
        index_type * queryIds_;
        size_type  * sequenceOffsets_;
        char       * sequences_;
        window_id  * maxWindowsInRange_;
        // output
        location_type   * queryResults_;
        int             * resultOffsets_;
        match_candidate * topCandidates_;

        cudaEvent_t resultsCopiedEvent_;
    };

    //---------------------------------------------------------------
    struct query_gpu_data
    {
        query_gpu_data(index_type maxQueries,
                       size_type maxSequenceLength,
                       size_type maxSketchSize,
                       size_type maxResultsPerQuery,
                       size_type maxCandidatesPerQuery,
                       bool multiGPU = false,
                       gpu_id gpuId = 0);
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
        location_type   * queryResultsTmp_;
        int             * resultOffsets_;
        int             * resultCounts_;

        int             * segBinCounters_;

        match_candidate * topCandidates_;

        cudaStream_t workStream_;
        cudaStream_t copyStream_;

        cudaEvent_t queryFinishedEvent_;
        cudaEvent_t sketchesCopiedEvent_;
        cudaEvent_t offsetsReadyEvent_;
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
                size_type maxResultsPerQuery,
                size_type maxCandidatesPerQuery,
                gpu_id numHostThreads,
                gpu_id numGPUs);
    //-----------------------------------------------------
    query_batch(const query_batch&) = delete;
    //---------------------------------------------------------------
    /** @brief free memory allocation */
    ~query_batch();

    //---------------------------------------------------------------
    gpu_id num_gpus() const noexcept {
        return numGPUs_;
    }
    //---------------------------------------------------------------
    query_host_data& host_data(gpu_id hostId) noexcept {
        return hostData_[hostId];
    }
    //---------------------------------------------------------------
    const query_gpu_data& gpu_data(gpu_id gpuId) const noexcept {
        return gpuData_[gpuId];
    }

    //---------------------------------------------------------------
    void sync_work_stream(gpu_id gpu_id);
    void sync_copy_stream(gpu_id gpu_id);

    //---------------------------------------------------------------
    /** @brief add read to host batch */
    template<class... Args>
    bool add_paired_read(gpu_id hostId, Args&&... args)
    {
        return hostData_[hostId].add_paired_read(
            std::forward<Args>(args)...,
            maxQueries_,
            maxSequenceLength_);
    }

    //---------------------------------------------------------------
    /** @brief asynchronously copy queries to device */
    void copy_queries_to_device_async(gpu_id hostId);
    void copy_queries_to_next_device_async(gpu_id hostId, gpu_id gpuId);
    //-----------------------------------------------------
    /** @brief synchronize event after copy to device 0 */
    void wait_for_queries_copied();
    //---------------------------------------------------------------
    /** @brief record event after sketch creation on device 0 */
    void mark_query_finished();
    //---------------------------------------------------------------
    /** @brief make work stream wait for allhits copy (before starting new query) */
    void wait_for_allhits_copied(gpu_id gpuId);

private:
    //---------------------------------------------------------------
    /** @brief asynchronously compact results in work stream */
    void compact_results_async(gpu_id hostId, gpu_id gpuId);
public:
    //---------------------------------------------------------------
    /**
     * @brief asynchronously compact and sort in work stream,
     *        and copy allhits to host in copy stream if needed
     */
    void compact_sort_and_copy_allhits_async(
        gpu_id hostId,
        gpu_id gpuId,
        bool copyAllHits);

    /**
     * @brief asynchronously generate top candidates in work stream
     *        and copy top candidates in copy stream
     */
    void generate_and_copy_top_candidates_async(
        gpu_id hostId,
        gpu_id gpuId,
        const ranked_lineage * lineages,
        taxon_rank lowestRank);


    //---------------------------------------------------------------
private:
    index_type maxQueries_;
    size_type  maxSequenceLength_;
    size_type  maxSketchSize_;
    size_type  maxResultsPerQuery_;
    size_type  maxCandidatesPerQuery_;

    std::vector<query_host_data> hostData_;
    std::vector<query_gpu_data> gpuData_;

    std::vector<segmented_sort> sorters_;

    cudaStream_t h2dCopyStream_;

    cudaEvent_t queriesCopiedEvent_;
    cudaEvent_t maxWinCopiedEvent_;

    gpu_id numGPUs_;
};


} // namespace mc


#endif
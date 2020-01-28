#ifndef MC_QUERY_BATCH_H_
#define MC_QUERY_BATCH_H_

#include <functional>

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
 template<class result_type>
class query_batch
{
    using id_type = uint32_t;

    using taxon_rank     = taxonomy::rank;
    using ranked_lineage = taxonomy::ranked_lineage;

public:
    //---------------------------------------------------------------
    /** @brief allocate memory on host and device */
    query_batch(id_type maxQueries,
                size_t maxEncodeLength,
                size_t maxResultsPerQuery,
                uint32_t maxCandidatesPerQuery);
    //-----------------------------------------------------
    query_batch(const query_batch&) = delete;
    //---------------------------------------------------------------
    /** @brief free memory allocation */
    ~query_batch();

    //---------------------------------------------------------------
    id_type num_segments_host() const noexcept {
        return h_numSegments_;
    }
    //---------------------------------------------------------------
    id_type num_segments_device() const noexcept {
        return d_numSegments_;
    }
    //---------------------------------------------------------------
    id_type num_queries_host() const noexcept {
        return h_numQueries_;
    }
    //---------------------------------------------------------------
    id_type num_queries_device() const noexcept {
        return d_numQueries_;
    }
    //---------------------------------------------------------------
    encodinglen_t * encode_offsets_device() const noexcept {
        return d_encodeOffsets_;
    }
    //---------------------------------------------------------------
    encodedseq_t * encoded_seq_device() const noexcept {
        return d_encodedSeq_;
    }
    //---------------------------------------------------------------
    encodedambig_t * encoded_ambig_device() const noexcept {
        return d_encodedAmbig_;
    }
    //---------------------------------------------------------------
    result_type * query_results_device() const noexcept {
        return d_queryResults_;
    }
    //---------------------------------------------------------------
    int * result_offsets_host() const noexcept {
        return h_resultOffsets_;
    }
    //---------------------------------------------------------------
    int * result_counts_device() const noexcept {
        return d_resultCounts_;
    }
    //---------------------------------------------------------------
    span<result_type> results(id_type id) const noexcept {
        if(id < d_numSegments_)
            return span<result_type>{
                h_queryResults_+h_resultOffsets_[id],
                h_queryResults_+h_resultOffsets_[id+1]
            };
        else
            return span<result_type>{};
    }
    //---------------------------------------------------------------
    span<match_candidate> top_candidates(id_type id) const noexcept {
        if(id < d_numSegments_)
            return span<match_candidate>{
                h_topCandidates_+id*maxCandidatesPerQuery_,
                h_topCandidates_+(id+1)*maxCandidatesPerQuery_
            };
        else
            return span<match_candidate>{};
    }
    //---------------------------------------------------------------
    const cudaStream_t& stream() const noexcept {
        return stream_;
    }
    //---------------------------------------------------------------
    void sync_streams();
    void sync_result_stream();

    //---------------------------------------------------------------
    void clear_host() noexcept {
        h_numSegments_ = 0;
        h_numQueries_ = 0;
    }

    void clear_device() noexcept {
        d_numSegments_ = 0;
        d_numQueries_ = 0;
    }

    //---------------------------------------------------------------
    /** @brief asynchronously copy queries to device using stream_ */
    void copy_queries_to_device_async();
    void wait_for_queries_copied();
    //---------------------------------------------------------------
    /** @brief asynchronously compact, sort and copy results to host using stream_ */
    void compact_sort_and_copy_results_async();
    //---------------------------------------------------------------
    /** @brief asynchronously generate and copy top candidates using stream_ */
    void generate_and_copy_top_candidates_async(
        const ranked_lineage * lineages,
        taxon_rank lowestRank);

    /*************************************************************************//**
    *
    * @brief add sequence to batch as windows of encoded characters
    *        if sequence does not fit into batch, don't add it
    *
    * @detail each window of a sequence is added as a separate query
    *
    * @return true if added, false otherwise
    *
    *****************************************************************************/
    template<class InputIterator>
    bool add_read(
        InputIterator first, InputIterator last,
        const sketcher& querySketcher
    ) {
        const numk_t kmerSize = querySketcher.kmer_size();
        const size_t windowSize = querySketcher.window_size();
        const size_t windowStride = querySketcher.window_stride();

        const size_t seqLength = std::distance(first, last);

        // no kmers in sequence, nothing to do here
        if(seqLength < kmerSize) return true;

        const window_id numWindows = (seqLength-kmerSize + windowStride) / windowStride;

        // batch full, nothing processed
        if(h_numQueries_ + numWindows > maxQueries_) return false;

        const auto availableBlocks = maxEncodeLength_ - h_encodeOffsets_[h_numQueries_];
        constexpr auto lettersPerBlock = sizeof(encodedambig_t)*CHAR_BIT;
        const auto blocksPerWindow = (windowSize-1) / lettersPerBlock + 1;
        // batch full, nothing processed
        if(availableBlocks < numWindows*blocksPerWindow) return false;

        // insert sequence into batch as separate windows
        for_each_window(first, last, windowSize, windowStride,
            [&] (InputIterator first, InputIterator last) {
                h_queryIds_[h_numQueries_] = h_numSegments_;
                h_encodeOffsets_[h_numQueries_+1] = h_encodeOffsets_[h_numQueries_];

                for_each_consecutive_substring_2bit<encodedseq_t>(first, last,
                    [&, this] (encodedseq_t substring, encodedambig_t ambig) {
                        auto& index = h_encodeOffsets_[h_numQueries_+1];
                        h_encodedSeq_[index] = substring;
                        h_encodedAmbig_[index] = ambig;
                        ++index;
                    });

                ++h_numQueries_;
            });

        ++h_numSegments_;

        return true;
    }

    //-----------------------------------------------------
    template<class Sequence>
    bool add_read(Sequence seq, const sketcher& querySketcher) {
        using std::begin;
        using std::end;

        return add_read(begin(seq), end(seq), querySketcher);
    }

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
        size_t insertSizeMax
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
            if(h_numQueries_ + 1 > maxQueries_) return false;

            // insert empty query
            h_queryIds_[h_numQueries_] = h_numSegments_;
            h_encodeOffsets_[h_numQueries_+1] = h_encodeOffsets_[h_numQueries_];

            ++h_numQueries_;
            ++h_numSegments_;
            return true;
        }

        const window_id numWindows1 = (seqLength1-kmerSize + windowStride) / windowStride;
        const window_id numWindows2 = (seqLength2-kmerSize + windowStride) / windowStride;

        // batch full, nothing processed
        if(h_numQueries_ + numWindows1 + numWindows2 > maxQueries_) return false;

        const auto availableBlocks = maxEncodeLength_ - h_encodeOffsets_[h_numQueries_];
        constexpr auto lettersPerBlock = sizeof(encodedambig_t)*CHAR_BIT;
        const auto blocksPerWindow = (windowSize-1) / lettersPerBlock + 1;
        // batch full, nothing processed
        if(availableBlocks < (numWindows1 + numWindows2)*blocksPerWindow) return false;

        // insert first sequence into batch as separate windows
        for_each_window(first1, last1, windowSize, windowStride,
            [&] (InputIterator first, InputIterator last) {
                if(distance(first, last) >= kmerSize) {
                    h_queryIds_[h_numQueries_] = h_numSegments_;
                    h_encodeOffsets_[h_numQueries_+1] = h_encodeOffsets_[h_numQueries_];

                    for_each_consecutive_substring_2bit<encodedseq_t>(first, last,
                        [&, this] (encodedseq_t substring, encodedambig_t ambig) {
                            auto& index = h_encodeOffsets_[h_numQueries_+1];
                            h_encodedSeq_[index] = substring;
                            h_encodedAmbig_[index] = ambig;
                            ++index;
                        }
                    );

                    ++h_numQueries_;
                }
            }
        );

        // insert second sequence into batch as separate windows
        for_each_window(first2, last2, windowSize, windowStride,
            [&] (InputIterator first, InputIterator last) {
                if(distance(first, last) >= kmerSize) {
                    h_queryIds_[h_numQueries_] = h_numSegments_;
                    h_encodeOffsets_[h_numQueries_+1] = h_encodeOffsets_[h_numQueries_];

                    for_each_consecutive_substring_2bit<encodedseq_t>(first, last,
                        [&, this] (encodedseq_t substring, encodedambig_t ambig) {
                            auto& index = h_encodeOffsets_[h_numQueries_+1];
                            h_encodedSeq_[index] = substring;
                            h_encodedAmbig_[index] = ambig;
                            ++index;
                        }
                    );

                    ++h_numQueries_;
                }
            }
        );

        h_maxWindowsInRange_[h_numSegments_] = window_id( 2 + (
            std::max(seqLength1 + seqLength2, insertSizeMax) /
            windowStride ));

        ++h_numSegments_;

        return true;
    }

    //-----------------------------------------------------
    template<class Sequence>
    bool add_paired_read(
        Sequence seq1, Sequence seq2,
        const sketcher& querySketcher,
        size_t insertSizeMax
    ) {
        using std::begin;
        using std::end;

        return add_paired_read(begin(seq1), end(seq1),
                               begin(seq2), end(seq2),
                               querySketcher,
                               insertSizeMax);
    }


private:
    id_type h_numSegments_;
    id_type d_numSegments_;
    id_type h_numQueries_;
    id_type d_numQueries_;
    id_type maxQueries_;

    size_t maxEncodeLength_;
    size_t maxResultsPerQuery_;
    uint32_t maxCandidatesPerQuery_;

    id_type        * h_queryIds_;
    id_type        * d_queryIds_;

    encodinglen_t  * h_encodeOffsets_;
    encodinglen_t  * d_encodeOffsets_;
    encodedseq_t   * h_encodedSeq_;
    encodedseq_t   * d_encodedSeq_;
    encodedambig_t * h_encodedAmbig_;
    encodedambig_t * d_encodedAmbig_;

    window_id      * h_maxWindowsInRange_;
    window_id      * d_maxWindowsInRange_;

    result_type    * h_queryResults_;
    result_type    * d_queryResults_;
    result_type    * d_queryResultsTmp_;

    int            * h_resultOffsets_;
    int            * d_resultOffsets_;
    int            * d_resultCounts_;
    int            * d_binnedSegIds_;

    int            * h_segBinCounters_;
    int            * d_segBinCounters_;

    match_candidate * h_topCandidates_;
    match_candidate * d_topCandidates_;

    cudaStream_t stream_;
    cudaStream_t resultCopyStream_;
    cudaEvent_t queriesCopiedEvent_;
    cudaEvent_t offsetsCopiedEvent_;
    cudaEvent_t resultReadyEvent_;
};


} // namespace mc


#endif
#ifndef MC_QUERY_BATCH_H_
#define MC_QUERY_BATCH_H_

#include <functional>

#include "cuda_runtime.h"

#include "config.h"
#include "hash_dna.h"

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
public:
    //---------------------------------------------------------------
    /** @brief allocate memory on host and device */
    query_batch(id_type maxQueries = 0, size_t maxEncodeLength = 0, size_t maxResultsPerQuery = 0);
    //-----------------------------------------------------
    query_batch(const query_batch&) = delete;
    //-----------------------------------------------------
    query_batch(query_batch&& other) {
        numSegments_        = other.numSegments_;
        maxQueries_         = other.maxQueries_;
        maxEncodeLength_    = other.maxEncodeLength_;
        numQueries_         = other.numQueries_;
        maxResultsPerQuery_ = other.maxResultsPerQuery_;

        other.numSegments_        = 0;
        other.maxQueries_         = 0;
        other.maxEncodeLength_    = 0;
        other.numQueries_         = 0;
        other.maxResultsPerQuery_ = 0;

        h_queryIds_      = other.h_queryIds_;
        d_queryIds_      = other.d_queryIds_;

        h_encodeOffsets_ = other.h_encodeOffsets_;
        d_encodeOffsets_ = other.d_encodeOffsets_;

        h_encodedSeq_    = other.h_encodedSeq_;
        d_encodedSeq_    = other.d_encodedSeq_;

        h_encodedAmbig_  = other.h_encodedAmbig_;
        d_encodedAmbig_  = other.d_encodedAmbig_;

        h_queryResults_  = other.h_queryResults_;
        d_queryResults_  = other.d_queryResults_;

        d_queryResultsTmp_ = other.d_queryResultsTmp_;

        h_resultOffsets_ = other.h_resultOffsets_;
        d_resultOffsets_ = other.d_resultOffsets_;

        d_resultCounts_  = other.d_resultCounts_;

        stream_ = other.stream_;
        other.stream_ = 0;
    };

    //---------------------------------------------------------------
    /** @brief free memory allocation */
    ~query_batch();

    //---------------------------------------------------------------
    id_type num_segments() const noexcept {
        return numSegments_;
    }
    //---------------------------------------------------------------
    id_type num_queries() const noexcept {
        return numQueries_;
    }
    //---------------------------------------------------------------
    id_type max_queries() const noexcept {
        return maxQueries_;
    }
    //---------------------------------------------------------------
    size_t max_encode_length() const noexcept {
        return maxEncodeLength_;
    }
    //---------------------------------------------------------------
    size_t max_results_per_query() const noexcept {
        return maxResultsPerQuery_;
    }

    //---------------------------------------------------------------
    id_type * query_ids_host() const noexcept {
        return h_queryIds_;
    }
    //---------------------------------------------------------------
    id_type * query_ids_device() const noexcept {
        return d_queryIds_;
    }
    //---------------------------------------------------------------
    encodinglen_t * encode_offsets_host() const noexcept {
        return h_encodeOffsets_;
    }
    //---------------------------------------------------------------
    encodinglen_t * encode_offsets_device() const noexcept {
        return d_encodeOffsets_;
    }
    //---------------------------------------------------------------
    encodedseq_t * encoded_seq_host() const noexcept {
        return h_encodedSeq_;
    }
    //---------------------------------------------------------------
    encodedseq_t * encoded_seq_device() const noexcept {
        return d_encodedSeq_;
    }
    //---------------------------------------------------------------
    encodedambig_t * encoded_ambig_host() const noexcept {
        return h_encodedAmbig_;
    }
    //---------------------------------------------------------------
    encodedambig_t * encoded_ambig_device() const noexcept {
        return d_encodedAmbig_;
    }
    //---------------------------------------------------------------
    result_type * query_results_host() const noexcept {
        return h_queryResults_;
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
    int * result_offsets_device() const noexcept {
        return d_resultOffsets_;
    }
    //---------------------------------------------------------------
    int * result_counts_device() const noexcept {
        return d_resultCounts_;
    }
    //---------------------------------------------------------------
    const cudaStream_t& stream() const noexcept {
        return stream_;
    }

    //---------------------------------------------------------------
    void clear() noexcept {
        numSegments_ = 0;
        numQueries_ = 0;
    }

    //---------------------------------------------------------------
    /** @brief asynchronously copy queries to device using stream_ */
    void copy_queries_to_device_async();
    //---------------------------------------------------------------
    /** @brief asynchronously copy results to host using stream_ */
    void copy_results_to_host_async();
    void copy_results_to_host_tmp_async();
    //---------------------------------------------------------------
    void sync_stream();

    //---------------------------------------------------------------
    void sort_results();

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
        if(numQueries_ + numWindows > maxQueries_) return false;

        const auto availableBlocks = maxEncodeLength_ - h_encodeOffsets_[numQueries_];
        constexpr auto lettersPerBlock = sizeof(encodedambig_t)*CHAR_BIT;
        const auto blocksPerWindow = (windowSize-1) / lettersPerBlock + 1;
        // batch full, nothing processed
        if(availableBlocks < numWindows*blocksPerWindow) return false;

        // insert sequence into batch as separate windows
        for_each_window(first, last, windowSize, windowStride,
            [&] (InputIterator first, InputIterator last) {
                h_queryIds_[numQueries_] = numSegments_;
                h_encodeOffsets_[numQueries_+1] = h_encodeOffsets_[numQueries_];

                for_each_consecutive_substring_2bit<encodedseq_t>(first, last,
                    [&, this] (encodedseq_t substring, encodedambig_t ambig) {
                        auto& index = h_encodeOffsets_[numQueries_+1];
                        h_encodedSeq_[index] = substring;
                        h_encodedAmbig_[index] = ambig;
                        ++index;
                    });

                ++numQueries_;
            });

        ++numSegments_;

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
        const sketcher& querySketcher
    ) {
        using std::distance;

        const numk_t kmerSize = querySketcher.kmer_size();
        const size_t windowSize = querySketcher.window_size();
        const size_t windowStride = querySketcher.window_stride();

        const size_t seqLength1 = distance(first1, last1);
        const size_t seqLength2 = distance(first2, last2);

        // no kmers in sequence, nothing to do here
        if(seqLength1 < kmerSize && seqLength2 < kmerSize) return true;

        const window_id numWindows1 = (seqLength1-kmerSize + windowStride) / windowStride;
        const window_id numWindows2 = (seqLength2-kmerSize + windowStride) / windowStride;

        // batch full, nothing processed
        if(numQueries_ + numWindows1 + numWindows2 > maxQueries_) return false;

        const auto availableBlocks = maxEncodeLength_ - h_encodeOffsets_[numQueries_];
        constexpr auto lettersPerBlock = sizeof(encodedambig_t)*CHAR_BIT;
        const auto blocksPerWindow = (windowSize-1) / lettersPerBlock + 1;
        // batch full, nothing processed
        if(availableBlocks < (numWindows1 + numWindows2)*blocksPerWindow) return false;

        // insert first sequence into batch as separate windows
        for_each_window(first1, last1, windowSize, windowStride,
            [&] (InputIterator first, InputIterator last) {
                if(distance(first, last) >= kmerSize) {
                    h_queryIds_[numQueries_] = numSegments_;
                    h_encodeOffsets_[numQueries_+1] = h_encodeOffsets_[numQueries_];

                    for_each_consecutive_substring_2bit<encodedseq_t>(first, last,
                        [&, this] (encodedseq_t substring, encodedambig_t ambig) {
                            auto& index = h_encodeOffsets_[numQueries_+1];
                            h_encodedSeq_[index] = substring;
                            h_encodedAmbig_[index] = ambig;
                            ++index;
                        }
                    );

                    ++numQueries_;
                }
            }
        );

        // insert second sequence into batch as separate windows
        for_each_window(first2, last2, windowSize, windowStride,
            [&] (InputIterator first, InputIterator last) {
                if(distance(first, last) >= kmerSize) {
                    h_queryIds_[numQueries_] = numSegments_;
                    h_encodeOffsets_[numQueries_+1] = h_encodeOffsets_[numQueries_];

                    for_each_consecutive_substring_2bit<encodedseq_t>(first, last,
                        [&, this] (encodedseq_t substring, encodedambig_t ambig) {
                            auto& index = h_encodeOffsets_[numQueries_+1];
                            h_encodedSeq_[index] = substring;
                            h_encodedAmbig_[index] = ambig;
                            ++index;
                        }
                    );

                    ++numQueries_;
                }
            }
        );

        ++numSegments_;

        return true;
    }

    //-----------------------------------------------------
    template<class Sequence>
    bool add_paired_read(
        Sequence seq1, Sequence seq2,
        const sketcher& querySketcher
    ) {
        using std::begin;
        using std::end;

        return add_paired_read(begin(seq1), end(seq1),
                               begin(seq2), end(seq2),
                               querySketcher);
    }


private:
    id_type numSegments_;
    id_type numQueries_;
    id_type maxQueries_;
    size_t maxEncodeLength_;
    size_t maxResultsPerQuery_;

    id_type        * h_queryIds_;
    id_type        * d_queryIds_;

    encodinglen_t  * h_encodeOffsets_;
    encodinglen_t  * d_encodeOffsets_;
    encodedseq_t   * h_encodedSeq_;
    encodedseq_t   * d_encodedSeq_;
    encodedambig_t * h_encodedAmbig_;
    encodedambig_t * d_encodedAmbig_;

    result_type    * h_queryResults_;
    result_type    * d_queryResults_;
    result_type    * d_queryResultsTmp_;

    int            * h_resultOffsets_;
    int            * d_resultOffsets_;
    int            * d_resultCounts_;

    cudaStream_t stream_;
};


} // namespace mc


#endif
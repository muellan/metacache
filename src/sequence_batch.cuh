#ifndef MC_SEQUENCE_BATCH_H_
#define MC_SEQUENCE_BATCH_H_

#include <functional>

#include "config.h"

namespace mc {

enum class policy {Host, Device};

/*************************************************************************//**
 *
 * @brief batch contains sequence data of multiple targets
 *        allocated memory location depends on policy
 *
 *****************************************************************************/
template<policy P>
class sequence_batch
 {
public:
    //---------------------------------------------------------------
    /**
     * @brief allocate memory on host or device depending on policy
     */
    sequence_batch(size_t maxTargets = 0, size_t maxEncodeLength = 0);
    //-----------------------------------------------------
    sequence_batch(const sequence_batch&) = delete;
    //-----------------------------------------------------
    sequence_batch(sequence_batch&& other) {
        maxTargets_      = other.maxTargets_;
        maxEncodeLength_ = other.maxEncodeLength_;
        numTargets_      = other.numTargets_;

        other.maxTargets_      = 0;
        other.maxEncodeLength_ = 0;
        other.numTargets_      = 0;

        targetIds_     = other.targetIds_;
        windowOffsets_ = other.windowOffsets_;
        encodeOffsets_ = other.encodeOffsets_;
        encodedSeq_    = other.encodedSeq_;
        encodedAmbig_  = other.encodedAmbig_;
    };

    //---------------------------------------------------------------
    /**
     * @brief free memory allocation
     */
    ~sequence_batch();

    //---------------------------------------------------------------
    size_t max_targets() const noexcept {
        return maxTargets_;
    }
    //---------------------------------------------------------------
    size_t max_encode_length() const noexcept {
        return maxEncodeLength_;
    }
    //---------------------------------------------------------------
    size_t num_targets() const noexcept {
        return numTargets_;
    }
    //-----------------------------------------------------
    void num_targets(size_t n) noexcept {
        if(n > max_targets()) n = max_targets();
        numTargets_ = n;
    }
    //---------------------------------------------------------------
    target_id * target_ids() const noexcept {
        return targetIds_;
    }
    //---------------------------------------------------------------
    window_id * window_offsets() const noexcept {
        return windowOffsets_;
    }
    //---------------------------------------------------------------
    encodinglen_t * encode_offsets() const noexcept {
        return encodeOffsets_;
    }
    //---------------------------------------------------------------
    encodedseq_t * encoded_seq() const noexcept {
        return encodedSeq_;
    }
    //---------------------------------------------------------------
    encodedambig_t * encoded_ambig() const noexcept {
        return encodedAmbig_;
    }

    /*************************************************************************//**
    *
    * @brief encode target sequence and add it to batch
    *
    * @detail if sequence does not fit into batch, only some windows of it are added
    *
    * @return number of processed windows
    *
    *****************************************************************************/
    template<class InputIterator, policy U = P, std::enable_if_t<U==policy::Host, int> = 0>
    window_id
    add_target(
        InputIterator first, InputIterator last,
        target_id tgt, window_id win,
        const sketcher& targetSketcher
    ) {
        const numk_t kmerSize     = targetSketcher.kmer_size();
        const size_t windowSize   = targetSketcher.window_size();
        const size_t windowStride = targetSketcher.window_stride();

        window_id processedWindows = 0;

        const size_t seqLength = std::distance(first, last);

        // no kmers in sequence, nothing to do here
        //TODO different case than batch full
        if(seqLength < kmerSize) return processedWindows;

        // batch full, nothing processed
        if(numTargets_ == maxTargets_) return processedWindows;

        const auto availableBlocks = (maxEncodeLength_ - encodeOffsets_[numTargets_]);
        // batch full, nothing processed
        if(!availableBlocks) return processedWindows;

        constexpr auto lettersPerBlock = sizeof(encodedambig_t)*CHAR_BIT;
        const size_t availableLength = availableBlocks * lettersPerBlock;

        InputIterator end = last;
        if(seqLength > availableLength) {
            //sequence does not fit into batch as a whole
            //so need to process a whole window

            const window_id availableWindows = availableLength < windowSize ? 0 :
                (availableLength - windowSize) / windowStride + 1;
            // batch full, nothing processed
            if(!availableWindows) return processedWindows;

            const size_t providedLength = windowSize + (availableWindows-1) * windowStride;

            //split sequence into [first,end] and [next,last] with overlap
            end = first + providedLength;
            processedWindows += availableWindows;
        }
        else {
            const window_id numWindows =
                (seqLength - kmerSize + windowStride) / windowStride;
            processedWindows += numWindows;
        }

        // insert sequence into batch
        targetIds_[numTargets_] = tgt;
        windowOffsets_[numTargets_] = win;
        encodeOffsets_[numTargets_+1] = encodeOffsets_[numTargets_];

        for_each_consecutive_substring_2bit<encodedseq_t>(first, end,
            [&, this] (encodedseq_t substring, encodedambig_t ambig) {
                auto& index = encodeOffsets_[numTargets_+1];
                encodedSeq_[index] = substring;
                encodedAmbig_[index] = ambig;
                ++index;
            });

        ++numTargets_;

        return processedWindows;
    }

private:
    size_t maxTargets_;
    size_t maxEncodeLength_;
    size_t numTargets_;

    target_id      * targetIds_;
    window_id      * windowOffsets_;
    encodinglen_t  * encodeOffsets_;
    encodedseq_t   * encodedSeq_;
    encodedambig_t * encodedAmbig_;
};

template<>
sequence_batch<policy::Host>::sequence_batch(size_t maxTargets, size_t maxEncodeLength);

template<>
sequence_batch<policy::Host>::~sequence_batch();

template<>
sequence_batch<policy::Device>::sequence_batch(size_t maxTargets, size_t maxEncodeLength);

template<>
sequence_batch<policy::Device>::~sequence_batch();


} // namespace mc


#endif
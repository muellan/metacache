#ifndef MC_SEQUENCE_BATCH_H_
#define MC_SEQUENCE_BATCH_H_

#include <functional>

#include "cuda_runtime.h"

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
    using index_type = uint32_t;
    using size_type  = uint32_t;

    //---------------------------------------------------------------
    /**
     * @brief allocate memory on host or device depending on policy
     */
    sequence_batch(index_type maxTargets = 0, size_type maxEncodeLength = 0);
    //-----------------------------------------------------
    sequence_batch(const sequence_batch&) = delete;
    //-----------------------------------------------------
    sequence_batch(sequence_batch&& other) {
        maxTargets_       = other.maxTargets_;
        maxSequenceLength_ = other.maxSequenceLength_;
        numTargets_       = other.numTargets_;

        other.maxTargets_       = 0;
        other.maxSequenceLength_ = 0;
        other.numTargets_       = 0;

        targetIds_       = other.targetIds_;
        windowOffsets_   = other.windowOffsets_;
        sequenceOffsets_ = other.sequenceOffsets_;
        sequence_       = other.sequence_;
    };

    //---------------------------------------------------------------
    /**
     * @brief free memory allocation
     */
    ~sequence_batch();

    //---------------------------------------------------------------
    void clear() noexcept {
        num_targets(0);
    }
    //---------------------------------------------------------------
    index_type max_targets() const noexcept {
        return maxTargets_;
    }
    //---------------------------------------------------------------
    index_type num_targets() const noexcept {
        return numTargets_;
    }
    //-----------------------------------------------------
    void num_targets(index_type n) noexcept {
        if(n > max_targets()) n = max_targets();
        numTargets_ = n;
    }
    //---------------------------------------------------------------
    size_type max_sequence_length() const noexcept {
        return maxSequenceLength_;
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
    size_type * sequence_offsets() const noexcept {
        return sequenceOffsets_;
    }
    //-----------------------------------------------------
    size_type sequence_length() const noexcept {
        return sequenceOffsets_[numTargets_];
    }
    //---------------------------------------------------------------
    char * sequence() const noexcept {
        return sequence_;
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

        const auto availableLength = (maxSequenceLength_ - sequenceOffsets_[numTargets_]);
        // batch full, nothing processed
        if(!availableLength) return processedWindows;

        InputIterator end = last;
        if(seqLength > availableLength) {
            // sequence does not fit into batch as a whole
            // but we need to process whole windows

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

        std::copy(first, end, sequence_ + sequenceOffsets_[numTargets_]);
        sequenceOffsets_[numTargets_+1] = sequenceOffsets_[numTargets_] + distance(first, end);

        ++numTargets_;

        return processedWindows;
    }

    friend
    void copy_host_to_device_async(
        const sequence_batch<policy::Host>& hostBatch,
        sequence_batch<policy::Device>& deviceBatch,
        cudaStream_t stream);

private:
    index_type maxTargets_;
    index_type numTargets_;
    size_type  maxSequenceLength_;

    target_id * targetIds_;
    window_id * windowOffsets_;
    size_type * sequenceOffsets_;
    char      * sequence_;
};

template<>
sequence_batch<policy::Host>::sequence_batch(index_type maxTargets, size_type maxEncodeLength);

template<>
sequence_batch<policy::Host>::~sequence_batch();

template<>
sequence_batch<policy::Device>::sequence_batch(index_type maxTargets, size_type maxEncodeLength);

template<>
sequence_batch<policy::Device>::~sequence_batch();


} // namespace mc


#endif
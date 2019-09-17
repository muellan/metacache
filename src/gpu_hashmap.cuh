#ifndef MC_GPU_HASH_MAP_H_
#define MC_GPU_HASH_MAP_H_

#include <functional>
#include <limits>

#include "config.h"

namespace mc {

enum class policy {Host, Device};

template<policy P = policy::Host>
class sequence_batch
 {
public:
    //---------------------------------------------------------------
    sequence_batch(size_t maxTargets = 0, size_t maxEncodeLength = 0);
    //-----------------------------------------------------
    sequence_batch(const sequence_batch&) = delete;
    //-----------------------------------------------------
    sequence_batch(sequence_batch&& other) {
        maxTargets_      = other.max_targets();
        other.max_targets(0);
        maxEncodeLength_ = other.max_encode_length();
        other.max_encode_length(0);
        numTargets_      = other.num_targets();
        other.num_targets(0);

        targetIds_     = other.target_ids();
        windowOffsets_ = other.window_offsets();
        encodeOffsets_ = other.encode_offsets();
        encodedSeq_    = other.encoded_seq();
        encodedAmbig_  = other.encoded_ambig();
    };

    //---------------------------------------------------------------
    ~sequence_batch();

    //---------------------------------------------------------------
    size_t max_targets() const noexcept {
        return maxTargets_;
    }
    //-----------------------------------------------------
    void max_targets(size_t n) noexcept {
        maxTargets_ = n;
    }
    //---------------------------------------------------------------
    size_t max_encode_length() const noexcept {
        return maxEncodeLength_;
    }
    //-----------------------------------------------------
    void max_encode_length(size_t n) noexcept {
        maxEncodeLength_ = n;
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

    //---------------------------------------------------------------
    template<class InputIterator, policy U = P, std::enable_if_t<U==policy::Host, int> = 0>
    InputIterator
    add_target(
        InputIterator first, InputIterator last,
        target_id tgt, window_id win
    ) {
        using std::distance;

        //todo check if (sub) sequences longer than kmer size

        if(numTargets_ < maxTargets_) {
            targetIds_[numTargets_] = tgt;
            windowOffsets_[numTargets_] = win;
            encodeOffsets_[numTargets_+1] = encodeOffsets_[numTargets_];

            constexpr auto lettersPerBlock = sizeof(encodedseq_t)*CHAR_BIT;
            const auto availableBlocks = (maxEncodeLength_ - encodeOffsets_[numTargets_]);
            const size_t availableLength = availableBlocks * lettersPerBlock;
            const size_t seqLength = distance(first, last);

            InputIterator end = (seqLength <= availableLength) ?
                                last :
                                first + availableLength;

            for_each_consecutive_substring_2bit<encodedseq_t>(first, end,
                [&, this] (encodedseq_t substring, encodedambig_t ambig) {
                    encodedSeq_[encodeOffsets_[numTargets_+1]] = substring;
                    encodedAmbig_[encodeOffsets_[numTargets_+1]] = ambig;
                    ++encodeOffsets_[numTargets_+1];
                });

            ++numTargets_;
            return end;
        }
        else {
            return first;
        }
    }

    // template<class InputRange, policy U = P, std::enable_if_t<U==policy::Host, int> = 0>
    // typename InputRange::const_iterator
    // add_target(
    //     InputRange input,
    //     target_id tgt, window_id win
    // ) {
    //     using std::begin;
    //     using std::end;
    //     return add_target(begin(input), end(input), tgt, win);
    // }

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


/*************************************************************************//**
 *
 * @brief   (integer) key -> value hashed multimap
 *          optimized for many values per key (pay attention to max_bucket_size()!
 *          Each bucket contains only one key and all values mapped to that key.
 *          Buckets may also be occupied AND empty ('key only').
 *
 * @details Hash conflicts are resolved with open addressing.
 *
 * @tparam  Key:    key type
 * @tparam  ValueT: value type
 * @tparam  Hash:   hash function (object) type
 * @tparam  ProbingIterator: implements probing scheme
 *
 *****************************************************************************/
template<
    class Key,
    class ValueT,
    class Hash = std::hash<Key>,
    class KeyEqual = std::equal_to<Key>,
    class BucketSizeT = std::uint8_t
>
class gpu_hashmap
{
    static_assert(std::is_integral<BucketSizeT>::value &&
        std::is_unsigned<BucketSizeT>::value,
        "bucket size type must be an unsigned integer");

public:
    //---------------------------------------------------------------
    using key_type         = Key;
    using value_type       = ValueT;
    using mapped_type      = ValueT;
    using hasher           = Hash;
    using key_equal        = KeyEqual;
    using bucket_size_type = BucketSizeT;

    using size_type        = size_t;

    //---------------------------------------------------------------
    static constexpr std::size_t
    max_bucket_size() noexcept {
        return std::numeric_limits<bucket_size_type>::max();
    }

private:
    //-----------------------------------------------------
    struct feature_batch {
        Key    * features_;
        ValueT * values_;
        size_t * featureCounter_;
    };

public:
    //---------------------------------------------------------------
    gpu_hashmap()
    :
        numKeys_(0), numValues_(0), maxLoadFactor_(default_max_load_factor()),
        hash_{}, keyEqual_{},
        kmerLength_(16), sketchSize_(16), windowSize_(128), windowStride_(113),
        seqBatchesDevice_()
    {
        init();
    }

    //-----------------------------------------------------
    gpu_hashmap(const key_equal& keyComp)
    :
        numKeys_(0), numValues_(0), maxLoadFactor_(default_max_load_factor()),
        hash_{}, keyEqual_{keyComp},
        kmerLength_(16), sketchSize_(16), windowSize_(128), windowStride_(113),
        seqBatchesDevice_()
    {
        init();
    }

    //-----------------------------------------------------
    gpu_hashmap(
        const hasher& hash,
        const key_equal& keyComp)
    :
        numKeys_(0), numValues_(0), maxLoadFactor_(default_max_load_factor()),
        hash_{hash}, keyEqual_{keyComp},
        kmerLength_(16), sketchSize_(16), windowSize_(128), windowStride_(113),
        seqBatchesDevice_()
    {
        init();
    }

private:
    void init();

public:
    //---------------------------------------------------------------
    size_type key_count() const noexcept {
        return numKeys_;
    }
    //-----------------------------------------------------
    size_type value_count() const noexcept {
        return numValues_;
    }

    //-----------------------------------------------------
    bool empty() const noexcept {
        return (numKeys_ < 1);
    }

    //---------------------------------------------------------------
    static constexpr float default_max_load_factor() noexcept {
        return 0.95;
    }
    //-----------------------------------------------------
    float max_load_factor() const noexcept {
        return maxLoadFactor_;
    }

    //---------------------------------------------------------------
    std::vector<Key> insert(const sequence_batch<policy::Host>& seqBatchHost);

    //---------------------------------------------------------------
    size_type numKeys_;
    size_type numValues_;
    float maxLoadFactor_;
    hasher hash_;
    key_equal keyEqual_;

    numk_t kmerLength_;
    sketch_size_type sketchSize_;
    size_t windowSize_;
    size_t windowStride_;

    size_t maxBatchNum_;
    std::vector<sequence_batch<policy::Device>> seqBatchesDevice_;
    feature_batch featureBatch_;
};


} // namespace mc


#endif

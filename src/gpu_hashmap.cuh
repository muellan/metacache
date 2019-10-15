#ifndef MC_GPU_HASH_MAP_H_
#define MC_GPU_HASH_MAP_H_

#include <functional>
#include <limits>
#include <memory>

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
    sequence_batch(size_t maxTargets = 0, size_t maxEncodeLength = 0);
    //-----------------------------------------------------
    sequence_batch(const sequence_batch&) = delete;
    //-----------------------------------------------------
    sequence_batch(sequence_batch&& other) {
        maxTargets_      = other.max_targets();
        maxEncodeLength_ = other.max_encode_length();
        numTargets_      = other.num_targets();
        other.max_targets(0);
        other.max_encode_length(0);
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

    /*************************************************************************//**
    *
    * @brief batch contains features and locations of multiple targets
    *        allocated memory location depends on policy
    *
    * @return number of processed windows
    *
    *****************************************************************************/
    template<class InputIterator, policy U = P, std::enable_if_t<U==policy::Host, int> = 0>
    window_id
    add_target(
        InputIterator first, InputIterator last,
        target_id tgt, window_id win,
        numk_t kmerSize,
        size_t windowSize, size_t windowStride
    ) {
        using std::distance;

        window_id processedWindows = 0;

        const size_t seqLength = distance(first, last);

        // no kmers in sequence, nothing to do here
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
                encodedSeq_[encodeOffsets_[numTargets_+1]] = substring;
                encodedAmbig_[encodeOffsets_[numTargets_+1]] = ambig;
                ++(encodeOffsets_[numTargets_+1]);
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




/*************************************************************************//**
 * TODO
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
    class single_value_hash_table;

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

public:
    //---------------------------------------------------------------
    gpu_hashmap();
    ~gpu_hashmap();
    gpu_hashmap(gpu_hashmap&&);
    gpu_hashmap(const gpu_hashmap&) = delete;

private:
    /*************************************************************************//**
    *
    * @brief batch contains features and locations of multiple targets
    *        allocated memory location depends on policy
    *
    *****************************************************************************/
    class feature_batch
    {
    public:
        using counter_type = uint32_t;

        //---------------------------------------------------------------
        feature_batch(counter_type maxFeatures = 0);
        //-----------------------------------------------------
        feature_batch(const feature_batch&) = delete;
        //-----------------------------------------------------
        feature_batch(feature_batch&& other) {
            maxFeatures_ = other.max_features();
            other.max_features(0);

            features_       = other.features();
            values_         = other.values();
            featureCounter_ = other.feature_counter();
        };

        //---------------------------------------------------------------
        ~feature_batch();

        //---------------------------------------------------------------
        counter_type max_features() const noexcept {
            return maxFeatures_;
        }
        //-----------------------------------------------------
        void max_features(counter_type n) noexcept {
            maxFeatures_ = n;
        }
        //---------------------------------------------------------------
        key_type * features() const noexcept {
            return features_;
        }
        //---------------------------------------------------------------
        value_type * values() const noexcept {
            return values_;
        }
        //---------------------------------------------------------------
        counter_type * feature_counter() const noexcept {
            return featureCounter_;
        }

    private:
        counter_type maxFeatures_;

        key_type     * features_;
        value_type   * values_;
        counter_type * featureCounter_;
    };

public:
    //---------------------------------------------------------------
    static constexpr std::size_t
    max_bucket_size() noexcept {
        return std::numeric_limits<bucket_size_type>::max();
    }

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
    std::vector<key_type> insert(const sequence_batch<policy::Host>& seqBatchHost);


    //---------------------------------------------------------------
    void insert_all(key_type *keys_in, value_type *values_in, size_type size_in);

private:
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
    std::vector<sequence_batch<policy::Device>> seqBatches_;
    std::vector<feature_batch> featureBatches_;

    std::unique_ptr<single_value_hash_table> hashTable_;
};


} // namespace mc


#endif

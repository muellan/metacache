#ifndef MC_GPU_HASH_MAP_H_
#define MC_GPU_HASH_MAP_H_

#include <functional>
#include <limits>

#include "config.h"

namespace mc {

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

public:
    //---------------------------------------------------------------
    gpu_hashmap()
    :
        numKeys_(0), numValues_(0), maxLoadFactor_(default_max_load_factor()),
        hash_{}, keyEqual_{}
    {
        init();
    }

    //-----------------------------------------------------
    gpu_hashmap(const key_equal& keyComp)
    :
        numKeys_(0), numValues_(0), maxLoadFactor_(default_max_load_factor()),
        hash_{}, keyEqual_{keyComp}
    {
        init();
    }

    //-----------------------------------------------------
    gpu_hashmap(
        const hasher& hash,
        const key_equal& keyComp)
    :
        numKeys_(0), numValues_(0), maxLoadFactor_(default_max_load_factor()),
        hash_{hash}, keyEqual_{keyComp}
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
    void insert(target_id tgt,
                std::vector<encodedseq_t> encodedSeq,
                std::vector<encodedambig_t> encodedAmbig,
                numk_t k);


    //---------------------------------------------------------------
    size_type numKeys_;
    size_type numValues_;
    float maxLoadFactor_;
    hasher hash_;
    key_equal keyEqual_;
};


} // namespace mc


#endif

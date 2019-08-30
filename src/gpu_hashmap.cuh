#ifndef MC_GPU_HASH_MAP_H_
#define MC_GPU_HASH_MAP_H_

#include <cuda.h>


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
    explicit
    gpu_hashmap()
    :
        numKeys_(0), numValues_(0), maxLoadFactor_(default_max_load_factor()),
        hash_{}, keyEqual_{}
    {
        cudaSetDevice(0);
    }

    //-----------------------------------------------------
    explicit
    gpu_hashmap(const key_equal& keyComp)
    :
        numKeys_(0), numValues_(0), maxLoadFactor_(default_max_load_factor()),
        hash_{}, keyEqual_{keyComp}
    {
        cudaSetDevice(0);
    }

    //-----------------------------------------------------
    explicit
    gpu_hashmap(const hasher& hash,
                  const key_equal& keyComp = key_equal{})
    :
        numKeys_(0), numValues_(0), maxLoadFactor_(default_max_load_factor()),
        hash_{hash}, keyEqual_{keyComp}
    {
        cudaSetDevice(0);
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
    size_type numKeys_;
    size_type numValues_;
    float maxLoadFactor_;
    hasher hash_;
    key_equal keyEqual_;
};


} // namespace mc


#endif

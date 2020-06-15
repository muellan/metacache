/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2020 André Müller (muellan@uni-mainz.de)
 *                       & Robin Kobus  (kobus@uni-mainz.de)
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

#ifndef MC_HASH_MAP_H_
#define MC_HASH_MAP_H_

#include <vector>
#include <algorithm>
#include <utility>
#include <functional>
#include <atomic>
#include <mutex>
#include <iostream>
#include <type_traits>
#include <memory>

#include "chunk_allocator.h"
#include "io_serialize.h"
#include "cmdline_utility.h"


namespace mc {

/*************************************************************************//**
 *
 * @brief helper templates for allocator configuration
 *
 *****************************************************************************/
namespace detail {

template<class T>
constexpr auto
check_supports_reserve(int)
    -> decltype(std::declval<T>().reserve(std::size_t(1)), std::true_type{});

template<class>
constexpr std::false_type
check_supports_reserve(char);

template<class T>
struct supports_reserve : public decltype(check_supports_reserve<T>(0)) {};

} // namespace detail

//-------------------------------------------------------------------
template<class Alloc, bool = detail::supports_reserve<Alloc>::value>
struct allocator_config
{
    static bool reserve(Alloc&, std::size_t) { return false; }
};

template<class Alloc>
struct allocator_config<Alloc,true>
{
    static bool reserve(Alloc& alloc, std::size_t n) {
        return alloc.reserve(n);
    }
};




/*************************************************************************//**
 *
 * @brief  iterator adapter for linear probing within a given iterator range
 *
 * @tparam RAIterator :  random access iterator
 *
 *****************************************************************************/
struct linear_probing
{
    template<class RAIterator>
    class iterator
    {
    public:
        explicit
        iterator(RAIterator pos, RAIterator beg, RAIterator end):
           pos_(pos), fst_(pos), beg_(beg), end_(end)
        {}

        decltype(auto) operator -> () const {
            return pos_.operator->();
        }

        decltype(auto) operator * () const {
            return *pos_;
        }

        iterator& operator ++ () noexcept {
            ++pos_;
            if(pos_ >= end_) {
                if(beg_ == end_) {
                    pos_ = end_;
                    return *this;
                }
                pos_ = beg_;
                end_ = fst_;
                beg_ = fst_;
            }
            return *this;
        }

        explicit operator bool() const noexcept {
            return (pos_ < end_);
        }
        explicit operator RAIterator () const {
            return pos_;
        }

    private:
        RAIterator pos_;
        RAIterator fst_;
        RAIterator beg_;
        RAIterator end_;
    };
};



/*************************************************************************//**
 *
 * @brief  iterator adapter for quadratic probing within a given iterator range
 *
 * @tparam RAIterator :  random access iterator
 *
 *****************************************************************************/
struct single_pass_quadratic_probing
{
    template<class RAIterator>
    class iterator
    {
    public:
        explicit
        iterator(RAIterator pos, RAIterator beg, RAIterator end):
           hop_(1), pos_(pos), fst_(pos), beg_(beg), end_(end)
        {}

        decltype(auto) operator -> () const {
            return pos_.operator->();
        }

        decltype(auto) operator * () const {
            return *pos_;
        }

        iterator& operator ++ () noexcept {
            pos_ += hop_;
            if(pos_ >= end_) {
                if(beg_ == end_) {
                    pos_ = end_;
                    return *this;
                }
                pos_ -= (end_ - beg_);
                end_ = fst_;
                beg_ = fst_;
            }
            ++hop_;
            return *this;
        }

        explicit operator bool() const noexcept {
            return (pos_ < end_);
        }
        explicit operator RAIterator () const {
            return pos_;
        }

    private:
        std::size_t hop_;
        RAIterator pos_;
        RAIterator fst_;
        RAIterator beg_;
        RAIterator end_;
    };
};



/*************************************************************************//**
 *
 * @brief   (integer) key -> value hashed multimap
 *          optimized for many values per key (pay attention to max_bucket_size()!
 *          The main difference to std::[unordered_][multi]map is that
 *          iterators will not iterate over {key,value} pairs but over buckets.
 *          Each bucket contains only one key and all values mapped to that key.
 *          Buckets may also be occupied AND empty ('key only').
 *
 * @details Hash conflicts are resolved with open addressing.
 *
 * @tparam  Key:    key type
 * @tparam  ValueT: value type
 * @tparam  Hash:   hash function (object) type
 * @tparam  ValueAllocator:  controls allocation of values
 * @tparam  BucketAllocator: controls allocation of buckets/keys
 * @tparam  ProbingIterator: implements probing scheme
 *
 *****************************************************************************/
template<
    class Key,
    class ValueT,
    class Hash = std::hash<Key>,
    class KeyEqual = std::equal_to<Key>,
    class ValueAllocator = chunk_allocator<ValueT>,
    class BucketAllocator = std::allocator<Key>,
    class BucketSizeT = std::uint8_t,
    class ProbingScheme = single_pass_quadratic_probing
>
class hash_multimap
{
    static_assert(std::is_integral<BucketSizeT>::value &&
                  std::is_unsigned<BucketSizeT>::value,
                  "bucket size type must be an unsigned integer");

    static_assert(std::is_trivially_constructible<ValueT>::value &&
                  std::is_trivially_destructible<ValueT>::value,
                  "value type must be trivially con- & destructible");

    using value_alloc  = std::allocator_traits<ValueAllocator>;
    using alloc_config = allocator_config<ValueAllocator>;

public:
    //---------------------------------------------------------------
    using key_type         = Key;
    using value_type       = ValueT;
    using mapped_type      = ValueT;
    using hasher           = Hash;
    using key_equal        = KeyEqual;
    using value_allocator  = ValueAllocator;
    using bucket_size_type = BucketSizeT;


    //---------------------------------------------------------------
    static constexpr std::size_t
    max_bucket_size() noexcept {
        return std::numeric_limits<bucket_size_type>::max();
    }


    /****************************************************************
     * @brief bucket = key + dynamic array
     */
    class bucket_type
    {
        friend class hash_multimap;

    public:
        using size_type        = bucket_size_type;
        using key_type         = hash_multimap::key_type;
        using value_type       = hash_multimap::value_type;
        using reference        = value_type&;
        using const_reference  = const value_type&;
        using iterator         = value_type*;
        using const_iterator   = const value_type*;

        bucket_type():
            values_{nullptr}, key_{},
            size_(0), capacity_(0)
        {}
        bucket_type(key_type&& key):
            values_{nullptr}, key_{std::move(key)},
            size_(0), capacity_(0)
        {}

        bool unused() const noexcept { return !values_; }
        bool empty()  const noexcept { return (size_ < 1); }

        size_type size()     const noexcept { return size_; }
        size_type capacity() const noexcept { return capacity_; }


        const key_type& key() const noexcept { return key_; }

        reference
        operator [](size_type i) noexcept { return values_[i]; }

        const_reference
        operator [](size_type i) const noexcept { return values_[i]; }

        //-------------------------------------------
              iterator  begin()       noexcept { return values_; }
        const_iterator  begin() const noexcept { return values_; }
        const_iterator cbegin() const noexcept { return values_; }

              iterator  end()       noexcept { return values_ + size_; }
        const_iterator  end() const noexcept { return values_ + size_; }
        const_iterator cend() const noexcept { return values_ + size_; }

    private:
        //-----------------------------------------------------
        template<class V>
        bool insert(value_allocator& alloc, V&& v) {
            if(!reserve(alloc, size_+1)) return false;
            values_[size_] = std::forward<V>(v);
            ++size_;
            return true;
        }
        //-------------------------------------------
        template<class InputIterator, class EndSentinel>
        bool insert(value_allocator& alloc,
                    InputIterator first, EndSentinel last)
        {
            using std::distance;
            auto nsize = size_ + size_type(distance(first,last));
            if(!reserve(alloc, nsize)) return false;
            for(auto p = values_+size_; first != last; ++p, ++first) {
                *p = *first;
            }
            size_ = nsize;
            return true;
        }
        //-------------------------------------------
        bool insert(value_allocator&,
                    value_type* values, size_type size, size_type capacity)
        {
            values_ = values;
            size_ = size;
            capacity_ = capacity;
            return true;
        }

        //-------------------------------------------
        void free(value_allocator& alloc) {
            if(!values_) return;
            deallocate(alloc);
            values_ = nullptr;
            size_ = 0;
            capacity_ = 0;
        }
        //-------------------------------------------
        void clear() {
            size_ = 0;
        }

        //-----------------------------------------------------
        void key(const key_type& key) noexcept { key_ = key; }

        //-----------------------------------------------------
        bool resize(value_allocator& alloc, size_type n) {
            if(size_ < n) {
                if(!reserve(alloc, n)) return false;
            }
            size_ = n;
            return true;
        }

        //-----------------------------------------------------
        void deallocate(value_allocator& alloc) {
            value_alloc::deallocate(alloc, values_, capacity_);
        }

        //-----------------------------------------------------
        /// @brief does not change size!
        bool reserve(value_allocator& alloc, std::size_t n) {
            if(n > max_bucket_size()) return false;

            if(!unused()) {
                if(n > capacity_) {
                    auto ncap = std::size_t(n + 0.3*size_);
                    if(ncap > max_bucket_size()) ncap = max_bucket_size();
                    //make new array and copy old values
                    auto nvals = value_alloc::allocate(alloc, ncap);
                    if(!nvals) return false;
                    std::copy(values_, values_ + size_, nvals);
                    deallocate(alloc);
                    values_ = nvals;
                    capacity_ = size_type(ncap);
                }
            }
            else {
                //make new array
                auto nvals = value_alloc::allocate(alloc, n);
                if(!nvals) return false;
                values_ = nvals;
                capacity_ = size_type(n);
            }
            return true;
        }

        //-----------------------------------------------------
        value_type* values_;
        key_type key_;
        size_type size_;
        size_type capacity_;
    };


    //---------------------------------------------------------------
    using bucket_allocator = typename
        std::allocator_traits<BucketAllocator>::template rebind_alloc<bucket_type>;

private:
    //-----------------------------------------------------
    using bucket_store_t = std::vector<bucket_type,bucket_allocator>;


public:
    //-----------------------------------------------------
    using size_type       = typename bucket_store_t::size_type;
    using reference       = bucket_type&;
    using const_reference = const bucket_type&;
    //-----------------------------------------------------
    using iterator        = typename bucket_store_t::iterator;
    using const_iterator  = typename bucket_store_t::const_iterator;
    //-----------------------------------------------------
    using local_iterator       = typename bucket_type::iterator;
    using const_local_iterator = typename bucket_type::const_iterator;


private:
    //-----------------------------------------------------
    using probing_iterator = typename ProbingScheme::template iterator<iterator>;


public:
    //---------------------------------------------------------------
    explicit
    hash_multimap(const value_allocator& valloc = value_allocator{},
                  const bucket_allocator& kalloc = bucket_allocator{})
    :
        numKeys_(0), numValues_(0),
        batchSize_(default_batch_size()),
        maxLoadFactor_(default_max_load_factor()),
        hash_{}, keyEqual_{}, alloc_{valloc},
        buckets_{kalloc}
    {
        buckets_.resize(1500007);
    }

    //-----------------------------------------------------
    explicit
    hash_multimap(const key_equal& keyComp,
                  const value_allocator& valloc = value_allocator{},
                  const bucket_allocator& kalloc = bucket_allocator{})
    :
        numKeys_(0), numValues_(0),
        batchSize_(default_batch_size()),
        maxLoadFactor_(default_max_load_factor()),
        hash_{}, keyEqual_{keyComp}, alloc_{valloc},
        buckets_{kalloc}
    {
        buckets_.resize(1500007);
    }

    //-----------------------------------------------------
    explicit
    hash_multimap(const hasher& hash,
                  const key_equal& keyComp = key_equal{},
                  const value_allocator& valloc = value_allocator{},
                  const bucket_allocator& kalloc = bucket_allocator{})
    :
        numKeys_(0), numValues_(0),
        batchSize_(default_batch_size()),
        maxLoadFactor_(default_max_load_factor()),
        hash_{hash}, keyEqual_{keyComp}, alloc_{valloc},
        buckets_{kalloc}
    {
        buckets_.resize(1500007);
    }

private:
    //---------------------------------------------------------------
    explicit
    hash_multimap(size_type numKeys,
                  const value_allocator& valloc = value_allocator{},
                  const bucket_allocator& kalloc = bucket_allocator{})
    :
        numKeys_(0), numValues_(0),
        batchSize_(default_batch_size()),
        maxLoadFactor_(default_max_load_factor()),
        hash_{}, keyEqual_{}, alloc_{valloc},
        buckets_{kalloc}
    {
        buckets_.resize(numKeys);
    }

public:
    //-------------------------------------------------------------------
    hash_multimap(const hash_multimap& src):
        numKeys_(0), numValues_(0),
        maxLoadFactor_(src.maxLoadFactor_),
        hash_{src.hash_}, keyEqual_{src.keyEqual_},
        alloc_{value_alloc::select_on_container_copy_construction(src.alloc_)},
        buckets_{}
    {
        reserve_keys(src.numKeys_);
        reserve_values(src.numValues_);

        for(const auto& b : src.buckets_) {
            if(!b.unused()) insert(b.key(), b.begin(), b.end());
        }
    }


    //-----------------------------------------------------
    hash_multimap(hash_multimap&& src):
        numKeys_(src.numKeys_), numValues_(src.numValues_),
        maxLoadFactor_(src.maxLoadFactor_),
        hash_{std::move(src.hash_)},
        keyEqual_{std::move(src.keyEqual_)},
        alloc_{std::move(src.alloc_)},
        buckets_{std::move(src.buckets_)}
    { }


    //-------------------------------------------------------------------
    hash_multimap& operator = (const hash_multimap& src)
    {
        hash_multimap tmp(src);
        swap(src);
        return *this;
    }

    //-----------------------------------------------------
    hash_multimap& operator = (hash_multimap&& src)
    {
        numKeys_ = src.numKeys_;
        numValues_ = src.numValues_;
        maxLoadFactor_ = src.maxLoadFactor_;
        hash_ = std::move(src.hash_);
        keyEqual_ = std::move(src.keyEqual_);
        alloc_ = std::move(src.alloc_);
        buckets_ = std::move(src.buckets_);
        return *this;
    }


    //-------------------------------------------------------------------
    ~hash_multimap() {
        for(auto& b : buckets_) {
            b.deallocate(alloc_);
        }
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
    size_type bucket_count() const noexcept {
        return buckets_.size();
    }
    //-----------------------------------------------------
    size_type buckets_capacity() const noexcept {
        return buckets_.capacity();
    }
    //-----------------------------------------------------
    size_type max_bucket_count() const noexcept {
        return buckets_.max_size();
    }

    //-----------------------------------------------------
    /**
     * @return number of buckets (and keys) with at least one value
     */
    size_type non_empty_bucket_count() const
    {
        auto n = size_type(0);
        for(const auto& b : buckets_) {
            if(!b.empty()) ++n;
        }
        return n;
    }


    //---------------------------------------------------------------
    size_type bucket_size(size_type i) const noexcept {
        return buckets_[i].size();
    }


    /****************************************************************
     * @brief if the value_allocator supports pre-allocation
     *        storage space for 'n' values will be allocated
     *
     * @param  n: number of values for which memory should be reserved
     * @return true on success
     */
    bool reserve_values(size_type n)
    {
        return alloc_config::reserve(alloc_, n);
    }


    /****************************************************************
     * @brief  rehashes to accomodate a specific number of keys
     * @param  n: number of keys for which memory should be reserved
     * @return true on success
     */
    bool reserve_keys(size_type n)
    {
        return rehash(size_type(1 + (1/max_load_factor() * n)));
    }


    /****************************************************************
     * @brief forces a specific number of buckets
     *
     * @param n: number of buckets in the hash table
     *
     * @return false, if n is less than the key count or is too small
     *         to keep the load factor below the allowed maximum
     */
    bool rehash(size_type n)
    {
        if(!rehash_possible(n)) return false;

        //make temporary new map
        //buckets resize might throw
        hash_multimap newmap{n};
        newmap.maxLoadFactor_ = maxLoadFactor_;

        //move old bucket contents into new hash slots
        //this should use only non-throwing operations
        for(auto& b : buckets_) {
            if(!b.unused()) {
                newmap.insert_into_slot(std::move(b.key_),
                                        b.values_, b.size_, b.capacity_);
            }
        }

        //should all be noexcept
        buckets_ = std::move(newmap.buckets_);
        hash_ = std::move(newmap.hash_);
        return true;
    }


    //---------------------------------------------------------------
    iterator
    insert(const key_type& key, const value_type& value)
    {
        make_sure_enough_buckets_left(1);
        return insert_into_slot(key, value);
    }
    iterator
    insert(const key_type& key, value_type&& value)
    {
        make_sure_enough_buckets_left(1);
        return insert_into_slot(key, std::move(value));
    }
    iterator
    insert(key_type&& key, const value_type& value)
    {
        make_sure_enough_buckets_left(1);
        return insert_into_slot(std::move(key), value);
    }
    iterator
    insert(key_type&& key, value_type&& value)
    {
        make_sure_enough_buckets_left(1);
        return insert_into_slot(std::move(key), std::move(value));
    }

    //-----------------------------------------------------
    template<class InputIterator, class EndSentinel>
    iterator
    insert(const key_type& key, InputIterator first, EndSentinel last)
    {
        make_sure_enough_buckets_left(1);
        return insert_into_slot(key, first, last);
    }
    template<class InputIterator, class EndSentinel>
    iterator
    insert(key_type&& key, InputIterator first, EndSentinel last)
    {
        make_sure_enough_buckets_left(1);
        return insert_into_slot(std::move(key), first, last);
    }


    /****************************************************************
     * @brief discards values with indices n-1 ... size-1
     *        in the bucket associated with 'key'
     */
    void
    shrink(const key_type& key, bucket_size_type n)
    {
        auto it = find(key);
        if(it != end()) shrink(it,n);
    }
    //-----------------------------------------------------
    void
    shrink(const_iterator it, bucket_size_type n)
    {
        if(it->size() > n) {
            const auto old = it->size();
            const_cast<bucket_type*>(&(*it))->resize(alloc_, n);
            numValues_ -= (old - it->size());
        }
    }


    /****************************************************************
     * @brief discards all values in bucket, but keeps bucket with key
     */
    void
    clear(const_iterator it)
    {
        if(!it->empty()) {
            auto n = it->size();
            const_cast<bucket_type*>(&(*it))->clear();
            numValues_ -= n;
        }
    }


    //---------------------------------------------------------------
    void clear()
    {
        if(numKeys_ < 1) return;

        //free bucket memory
        for(auto& b : buckets_) {
            b.free(alloc_);
        }
        numKeys_ = 0;
        numValues_ = 0;
    }


    //---------------------------------------------------------------
    /**
     * @brief very dangerous!
     *        clears buckets without memory deallocation
     */
    void clear_without_deallocation()
    {
        for(auto& b : buckets_) {
            b.values_ = nullptr;
            b.size_ = 0;
            b.capacity_ = 0;
        }
        numKeys_ = 0;
        numValues_ = 0;
    }


    //---------------------------------------------------------------
    const_reference
    bucket(size_type i) const noexcept {
        return buckets_[i];
    }
    //-----------------------------------------------------
    reference
    bucket(size_type i) noexcept {
        return buckets_[i];
    }


    //---------------------------------------------------------------
    /**
     * @return  iterator to bucket with all values of a key
     *          or end iterator if key not found
     */
    iterator
    find(const key_type& key)
    {
        return find_occupied_slot(key);
    }
    //-----------------------------------------------------
    const_iterator
    find(const key_type& key) const
    {
        return const_iterator(
            const_cast<hash_multimap*>(this)->find_occupied_slot(key) );
    }

    //-----------------------------------------------------
    size_type
    count(const key_type& key) const {
        auto it = find_occupied_slot(key);
        return (it != end()) ? it->size() : 0;
    }


    //---------------------------------------------------------------
    float load_factor() const noexcept {
        return (numKeys_ / float(buckets_.size()));
    }

    //---------------------------------------------------------------
    static constexpr float default_max_load_factor() noexcept {
        return 0.95;
    }
    //-----------------------------------------------------
    float max_load_factor() const noexcept {
        return maxLoadFactor_;
    }
    //-----------------------------------------------------
    void max_load_factor(float lf) {
        if(lf > 1.0f) lf = 1.0f;
        if(lf < 0.1f) lf = 0.1f;

        using std::abs;
        if(abs(maxLoadFactor_ - lf) > 0.00001f) {
            maxLoadFactor_ = lf;
            if(load_factor() > maxLoadFactor_) {
                rehash(size_type(1.1 * buckets_.size() * maxLoadFactor_ + 2));
            }
        }
    }


    //---------------------------------------------------------------
    const value_allocator&
    get_value_allocator() const noexcept {
        return alloc_;
    }
    decltype(auto)
    get_bucket_allocator() const noexcept {
        return buckets_.get_allocator();
    }


    //---------------------------------------------------------------
    const hasher&
    hash_function() const noexcept {
        return hash_;
    }


    //---------------------------------------------------------------
    const key_equal&
    key_eq() const noexcept {
        return keyEqual_;
    }


    //---------------------------------------------------------------
    iterator
     begin() noexcept { return buckets_.begin(); }
    //-----------------------------------------------------
    const_iterator
     begin() const noexcept { return buckets_.begin(); }
    //-----------------------------------------------------
    const_iterator
    cbegin() const noexcept { return buckets_.begin(); }

    //-----------------------------------------------------
    iterator
     end() noexcept { return buckets_.end(); }
    //-----------------------------------------------------
    const_iterator
     end() const noexcept { return buckets_.end(); }
    //-----------------------------------------------------
    const_iterator
    cend() const noexcept { return buckets_.end(); }

    //---------------------------------------------------------------
    local_iterator
     begin(size_type i) noexcept { return buckets_[i].begin(); }
    //-----------------------------------------------------
    const_local_iterator
     begin(size_type i) const noexcept { return buckets_[i].begin(); }
    //-----------------------------------------------------
    const_local_iterator
    cbegin(size_type i) const noexcept { return buckets_[i].begin(); }

    //-----------------------------------------------------
    local_iterator
     end(size_type i) noexcept { return buckets_[i].end(); }
    //-----------------------------------------------------
    const_local_iterator
     end(size_type i) const noexcept { return buckets_[i].end(); }
    //-----------------------------------------------------
    const_local_iterator
    cend(size_type i) const noexcept { return buckets_[i].end(); }


    //---------------------------------------------------------------
    void swap(hash_multimap& other)
    {
        std::swap(numKeys_, other.numKeys_);
        std::swap(numValues_, other.numValues_);
        std::swap(maxLoadFactor_, other.maxLoadFactor_);
        std::swap(hash_, other.hash_);
        std::swap(keyEqual_, other.keyEqual_);
        std::swap(alloc_, other.alloc_);
        std::swap(buckets_, other.buckets_);
    }


    /****************************************************************
     * @brief deserialize hashmap from input stream
     */
    friend void read_binary(std::istream& is, hash_multimap& m) {
        m.deserialize(is);
    }

    /****************************************************************
     * @brief serialize hashmap to output stream
     */
    friend void write_binary(std::ostream& os, const hash_multimap& m) {
        m.serialize(os);
    }


    //---------------------------------------------------------------
    static constexpr size_type default_batch_size() noexcept {
        return size_type(1) << 20;
    }
    //-----------------------------------------------------
    size_type batch_size() const noexcept {
        return batchSize_;
    }


private:
    //---------------------------------------------------------------
    std::uint64_t deserialize_batch_of_buckets(
        std::istream& is,
        std::vector<key_type>& keyBuffer,
        std::vector<bucket_size_type>& sizeBuffer,
        std::uint64_t batchSize,
        value_type * valuesOffset)
    {
        auto batchValuesOffset = valuesOffset;

        //load batch
        read_binary(is, keyBuffer.data(), batchSize);
        read_binary(is, sizeBuffer.data(), batchSize);

        //insert batch
        for(std::uint64_t i = 0; i < batchSize; ++i) {
            const auto& bucketSize = sizeBuffer[i];

            if(bucketSize > 0) {
                const auto& key = keyBuffer[i];

                auto it = insert_into_slot(key, valuesOffset, bucketSize, bucketSize);
                if(it == buckets_.end())
                    std::cerr << "could not insert key " << key << '\n';

                valuesOffset += bucketSize;
            }
        }
        std::uint64_t batchValuesCount = valuesOffset - batchValuesOffset;
        read_binary(is, batchValuesOffset, batchValuesCount);

        return batchValuesCount;
    }


    //---------------------------------------------------------------
    void deserialize(std::istream& is)
    {
        using len_t = std::uint64_t;

        clear();
        std::cerr << '\n';
        show_progress_indicator(std::cerr, 0);

        len_t nkeys = 0;
        read_binary(is, nkeys);
        len_t nvalues = 0;
        read_binary(is, nvalues);

        len_t totalSize = nkeys*(sizeof(key_type)+sizeof(bucket_size_type))
                        + nvalues*sizeof(value_type);
        len_t indicator = 0;

        len_t batchSize = 0;
        read_binary(is, batchSize);

        if(nkeys > 0) {
            //if the allocator supports it: reserve one large memory chunk
            //for all values; individual buckets will then point into this
            //array; the default chunk_allocator does this
            reserve_keys(nkeys);
            reserve_values(nvalues);
            auto valuesOffset = alloc_.allocate(nvalues);

            {// read keys & bucket sizes & values in batches
                const len_t numFullBatches = nkeys / batchSize;
                const len_t lastBatchSize = nkeys % batchSize;

                std::vector<key_type> keyBuffer(batchSize);
                std::vector<bucket_size_type> sizeBuffer(batchSize);

                for(len_t b = 0; b < numFullBatches; ++b) {
                    auto batchValuesCount = deserialize_batch_of_buckets(
                        is, keyBuffer, sizeBuffer, batchSize, valuesOffset);

                    indicator += batchSize*(sizeof(key_type)+sizeof(bucket_size_type))
                               + batchValuesCount*sizeof(value_type);
                    show_progress_indicator(std::cerr, float(indicator)/totalSize);

                    valuesOffset += batchValuesCount;
                }

                deserialize_batch_of_buckets(
                    is, keyBuffer, sizeBuffer, lastBatchSize, valuesOffset);
            }

            numKeys_ = nkeys;
            numValues_ = nvalues;
            batchSize_ = batchSize;
        }

        clear_current_line(std::cerr);
    }


    //---------------------------------------------------------------
    /**
     * @brief binary serialization of all non-emtpy buckets
     */
    void serialize(std::ostream& os) const
    {
        using len_t = std::uint64_t;

        const len_t nonEmptyBucketCount = non_empty_bucket_count();

        write_binary(os, len_t(nonEmptyBucketCount));
        write_binary(os, len_t(value_count()));
        write_binary(os, len_t(batch_size()));

        const auto batchSize = batch_size();
        const len_t avgValueCount = value_count() / nonEmptyBucketCount;

        {// write keys & bucket sizes & values in batches
            std::vector<key_type> keyBuffer;
            keyBuffer.reserve(batchSize);
            std::vector<bucket_size_type> sizeBuffer;
            sizeBuffer.reserve(batchSize);
            std::vector<value_type> valBuffer;
            valBuffer.reserve(batchSize*avgValueCount);
            
            for(const auto& bucket : buckets_) {
                if(!bucket.empty()) {
                    keyBuffer.emplace_back(bucket.key());
                    sizeBuffer.emplace_back(bucket.size());
                    std::copy(bucket.begin(), bucket.end(), std::back_inserter(valBuffer));

                    if(keyBuffer.size() == batchSize) {
                        // store batch
                        write_binary(os, keyBuffer.data(), keyBuffer.size());
                        write_binary(os, sizeBuffer.data(), sizeBuffer.size());
                        write_binary(os, valBuffer.data(), valBuffer.size());
                        // reset batch
                        keyBuffer.clear();
                        sizeBuffer.clear();
                        valBuffer.clear();
                    }
                }
            }

            if(keyBuffer.size() > 0) {
                // store last batch
                write_binary(os, keyBuffer.data(), keyBuffer.size());
                write_binary(os, sizeBuffer.data(), sizeBuffer.size());
                write_binary(os, valBuffer.data(), valBuffer.size());
            }
        }
    }


    //---------------------------------------------------------------
    iterator
    find_occupied_slot(const key_type& key)
    {
        probing_iterator it {
            buckets_.begin() + (hash_(key) % buckets_.size()),
            buckets_.begin(), buckets_.end()};

        //find bucket
        do {
            if(it->unused()) return buckets_.end();
            if(keyEqual_(it->key(), key)) return iterator(it);
        } while(++it);

        return buckets_.end();
    }

    //-----------------------------------------------------
    template<class... Values>
    iterator
    insert_into_slot(key_type key, Values&&... newvalues)
    {
        probing_iterator it {
            buckets_.begin() + (hash_(key) % buckets_.size()),
            buckets_.begin(), buckets_.end()};

        do {
            //empty slot found
            if(it->unused()) {
                if(it->insert(alloc_, std::forward<Values>(newvalues)...)) {
                    it->key_ = std::move(key);
                    ++numKeys_;
                    numValues_ += it->size();
                    return iterator(it);
                }
                //could not insert
                return buckets_.end();
            }
            //key already inserted
            if(keyEqual_(it->key(), key)) {
                auto oldsize = it->size();
                if(it->insert(alloc_, std::forward<Values>(newvalues)...)) {
                    numValues_ += it->size() - oldsize;
                    return iterator(it);
                }
                return buckets_.end();
            }
        } while(++it);

        return buckets_.end();
    }

    //---------------------------------------------------------------
    void make_sure_enough_buckets_left(size_type more)
    {
        auto n = numKeys_ + more;

        if( (n / float(buckets_.size()) > maxLoadFactor_ ) ||
            (n >= buckets_.size()) )
        {
            rehash(std::max(
                size_type(1 + 1.8 * n),
                size_type(1 + 1.8 * (n / maxLoadFactor_)) ));
        }
    }


    //---------------------------------------------------------------
    bool rehash_possible(size_type n) const noexcept
    {
        //number of buckets must be greater or equal to the number of keys
        if(n == bucket_count() || n < key_count()) return false;

        //make sure we stay below the maximum load factor
        auto newload = (load_factor() * (float(n)/bucket_count()));
        if(n < bucket_count() && newload > max_load_factor()) return false;
        return true;
    }


    //---------------------------------------------------------------
    size_type numKeys_;
    size_type numValues_;
    size_type batchSize_;
    float maxLoadFactor_;
    hasher hash_;
    key_equal keyEqual_;
    value_allocator alloc_;
    bucket_store_t buckets_;
};



} // namespace mc


#endif

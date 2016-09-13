/*****************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * version 1.0
 *
 * Copyright (C) 2016 André Müller (muellan@uni-mainz.de)
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


namespace mc {

/*****************************************************************************
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
    static void reserve(Alloc&, std::size_t) {}
};

template<class Alloc>
struct allocator_config<Alloc,true>
{
    static void reserve(Alloc& alloc, std::size_t n) {
        alloc.reserve(n);
    }
};




/*****************************************************************************
 *
 * @brief   (integer) key -> value hashed multimap
 *          optimized for many values per key (pay attention to max_bucket_size()!
 *          The main difference to std::[unordered_][multi]map is that
 *          iterators will not iterate over {key,value} pairs but over buckets.
 *          Each bucket contains only one key and all values mapped to that key.
 *
 * @details Hash conflicts are resolved through open addressing with
 *          linear probing.
 *          Deletion uses 'back shifting' instead of 'tombstoning'.
 *
 * @tparam  Key:    key type
 * @tparam  ValueT: value type
 * @tparam  Hash:   hash function (object) type
 * @tparam  ValueAllocator:  controls allocation of values
 * @tparam  BucketAllocator: controls allocation of buckets/keys
 *
 *****************************************************************************/
template<
    class Key,
    class ValueT,
    class Hash = std::hash<Key>,
    class KeyEqual = std::equal_to<Key>,
    class ValueAllocator = chunk_allocator<ValueT>,
    class BucketAllocator = std::allocator<Key>,
    class BucketSizeT = std::uint8_t //256
>
class hash_multimap
{
    static_assert(std::is_integral<BucketSizeT>::value &&
                  std::is_unsigned<BucketSizeT>::value,
                  "bucket size type must be an unsigned integer");

    using value_alloc  = std::allocator_traits<ValueAllocator>;
    using alloc_config = allocator_config<ValueAllocator>;

    using probelen_t = std::uint8_t;


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

    //-----------------------------------------------------
    /**
     * @brief bucket = key + dynamic array
     *        uses 16-bit unsigned size and capcity counters
     *        to avoid padding losses when using a 32-bit key
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
            size_(0), capacity_(0)//, probelen_(0)
        {}
        bucket_type(key_type&& key):
            values_{nullptr}, key_{std::move(key)},
            size_(0), capacity_(0)//, probelen_(0)
        {}

        size_type size()     const noexcept { return size_; }
        size_type capacity() const noexcept { return capacity_; }

        bool empty() const noexcept { return !values_; }

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

//        probelen_t probe_length() const noexcept { return probelen_; }
        probelen_t probe_length() const noexcept { return 0; }

    private:
        //-----------------------------------------------------
        template<class V>
        bool insert(V&& v, value_allocator& alloc) {
            if(!reserve(size_+1, alloc)) return false;
            values_[size_] = std::forward<V>(v);
            ++size_;
            return true;
        }
        //-------------------------------------------
        template<class InputIterator, class EndSentinel>
        bool insert(InputIterator first, EndSentinel last,
                    value_allocator& alloc)
        {
            using std::distance;
            auto nsize = size_ + size_type(distance(first,last));
            if(!reserve(nsize, alloc)) return false;
            for(auto p = values_+size_; first != last; ++p, ++first) {
                *p = *first;
            }
            size_ = nsize;
            return true;
        }

        //-------------------------------------------
        void clear(value_allocator& alloc) {
            if(empty()) return;
            if(!std::is_pod<value_type>::value) {
                for(auto i = values_, e = i + size_; i < e; ++i) {
                    value_alloc::destroy(alloc, i);
                }
            }
            value_alloc::deallocate(alloc, values_, capacity_);
            values_ = nullptr;
            size_ = 0;
            capacity_ = 0;
//            probelen_ = 0;
        }

        //-----------------------------------------------------
        void key(const key_type& key) noexcept { key_ = key; }

        //-----------------------------------------------------
        bool resize(size_type n, value_allocator& alloc) {
            if(size_ < n) {
                return reserve(n, alloc);
            }
            else if(size_ > n) {
                size_ = n;
            }
            return true;
        }

        //-----------------------------------------------------
        /// @brief does not change size!
        bool reserve(std::size_t n, value_allocator& alloc) {
            if(n > max_bucket_size()) return false;

            if(values_) {
                if(n > capacity_) {
                    auto ncap = std::size_t(n + 0.3*size_);
                    if(ncap > max_bucket_size()) ncap = max_bucket_size();
                    //make new array
                    auto nvals = value_alloc::allocate(alloc, ncap);
                    if(!std::is_pod<value_type>::value) {
                        for(auto i = nvals, e = i + ncap; i < e; ++i) {
                            value_alloc::construct(alloc, i);
                        }
                    }
                    //move values over
                    auto nv = nvals;
                    auto ov = values_;
                    for(size_type i = 0; i < size_; ++i, ++ov, ++nv) {
                        *nv = std::move(*ov);
                    }
                    clear(alloc);
                    values_ = nvals;
                    capacity_ = size_type(ncap);
                }
                else if(!std::is_pod<value_type>::value) {
                    for(auto i = values_ + size_, e = values_ + n; i < e; ++i) {
                        value_alloc::construct(alloc, i);
                    }
                }
            }
            else {
                values_ = value_alloc::allocate(alloc, n);
                if(!std::is_pod<value_type>::value) {
                    for(auto i = values_, e = i + n; i < e; ++i) {
                        value_alloc::construct(alloc, i);
                    }
                }
                capacity_ = n;
            }
            return true;
        }

        //-----------------------------------------------------
        value_type* values_;
        key_type key_;
        size_type size_;
        size_type capacity_;
//        probelen_t probelen_;
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


    //---------------------------------------------------------------
    explicit
    hash_multimap(const value_allocator& valloc = value_allocator{},
                  const bucket_allocator& kalloc = bucket_allocator{})
    :
        numKeys_(0), numValues_(0), maxLoadFactor_(default_max_load_factor()),
        hash_{}, keyEqual_{}, alloc_{valloc},
        buckets_{kalloc}
    {
        buckets_.resize(5);
    }

    //-----------------------------------------------------
    explicit
    hash_multimap(const key_equal& keyComp,
                  const value_allocator& valloc = value_allocator{},
                  const bucket_allocator& kalloc = bucket_allocator{})
    :
        numKeys_(0), numValues_(0), maxLoadFactor_(default_max_load_factor()),
        hash_{}, keyEqual_{keyComp}, alloc_{valloc},
        buckets_{kalloc}
    {
        buckets_.resize(5);
    }

    //-----------------------------------------------------
    explicit
    hash_multimap(const hasher& hash,
                  const key_equal& keyComp = key_equal{},
                  const value_allocator& valloc = value_allocator{},
                  const bucket_allocator& kalloc = bucket_allocator{})
    :
        numKeys_(0), numValues_(0), maxLoadFactor_(default_max_load_factor()),
        hash_{hash}, keyEqual_{keyComp}, alloc_{valloc},
        buckets_{kalloc}
    {
        buckets_.resize(5);
    }


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
            if(!b.empty()) insert(b.key(), b.begin(), b.end());
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
        clear();
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
        return (numValues_ < 1);
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


    //---------------------------------------------------------------
    size_type bucket_size(size_type i) const noexcept {
        return buckets_[i].size();
    }


    /****************************************************************
     * @brief if the value_allocator supports pre-allocation
     *        storage space for 'n' values will be allocated
     *
     * @param n: number of values for which memory should be reserved
     */
    void reserve_values(size_type n)
    {
        alloc_config::reserve(alloc_, n);
    }


    //---------------------------------------------------------------
    void reserve_keys(size_type nkeys)
    {
        rehash(size_type(1 + (1/max_load_factor() * nkeys)));
    }


    /****************************************************************
     * @brief forces
     *
     * @param n: number of buckets in the hash table
     */
    bool rehash(size_type n)
    {
        //number of buckets must be greater or equal to the number of keys
        if(n == bucket_count() || n < key_count()) return false;

        //make sure we stay below the maximum load factor
        auto newload = (load_factor() * (float(n)/bucket_count()));
        if(n < bucket_count() && newload > max_load_factor()) return false;

//        std::cout << "\nREHASH to " << n << std::endl;

        //make temporary new map
        hash_multimap newmap{};
        newmap.maxLoadFactor_ = maxLoadFactor_;
        //this might throw
        newmap.buckets_.resize(n);

        //move old bucket contents into new hash slots
        //this should use only non-throwing operations
        for(auto& b : buckets_) {
            if(!b.empty()) {
                auto it = newmap.find_slot_for_insertion(std::move(b.key_));
                it->values_ = b.values_;
                it->size_ = b.size_;
                it->capacity_ = b.capacity_;
                //probe length was set by find_slot_for_insertion
            }
        }

        //should all be noexcept
        buckets_ = std::move(newmap.buckets_);
        hash_ = std::move(newmap.hash_);
        return true;
    }


    //---------------------------------------------------------------
    template<class K, class V>
    iterator
    insert(K&& key, V&& value)
    {
        make_sure_enough_buckets_left(1);

        auto it = find_slot_for_insertion(std::forward<K>(key));

        if(it->empty()) {
            if(!it->insert(std::forward<V>(value), alloc_)) {
//                it->probelen_ = 0;
                return end();
            }
            ++numKeys_;
            numValues_ += it->size();
            return it;
        }
        else { //key already inserted
            auto oldsize = it->size();
            if(!it->insert(std::forward<V>(value), alloc_)) return end();
            numValues_ += it->size() - oldsize;
            return it;
        }
    }


    //---------------------------------------------------------------
    template<class K, class InputIterator, class EndSentinel>
    iterator
    insert(K&& key, InputIterator first, EndSentinel last)
    {
        make_sure_enough_buckets_left(1);

        auto it = find_slot_for_insertion(std::forward<K>(key));

        if(it->empty()) {
            if(!it->insert(first, last, alloc_)) {
//                it->probelen_ = 0;
                return end();
            }
            ++numKeys_;
            numValues_ += it->size();
            return it;
        }
        else { //key already inserted
            auto oldsize = it->size();
            if(!it->insert(first, last, alloc_)) return end();
            numValues_ += it->size() - oldsize;
            return it;
        }
    }


    //---------------------------------------------------------------
    size_type
    erase(const key_type& key)
    {
        auto it = find(key);
        return (it != end()) ? erase_slot(it) : 0;
    }
    //-----------------------------------------------------
    size_type
    erase(iterator it)
    {
        using std::distance;
        return erase_slot(it);
    }
    //-----------------------------------------------------
    size_type
    erase(const_iterator it)
    {
        using std::distance;
        return erase_slot(cbegin() + distance(cbegin(),it));
    }


    /****************************************************************
     * @brief discards values with indices n-1 ... size-1
     *        in the bucket associated with 'key'
     */
    void
    shrink(const key_type& key, bucket_size_type n)
    {
        auto it = find(key);
        if(it != end()) shrink(it);
    }
    //-----------------------------------------------------
    void
    shrink(const_iterator it, bucket_size_type n)
    {
        if(it->size() > n) const_cast<bucket_type*>(&(*it))->resize(n, alloc_);
    }


    //---------------------------------------------------------------
    void clear()
    {
        if(numKeys_ < 1) return;

        //free bucket memory
        for(auto& b : buckets_) {
            b.clear(alloc_);
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
        return 0.80;
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
    auto
    get_bucket_allocator() const noexcept
        -> decltype(std::declval<bucket_store_t>().get_allocator())
    {
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


    //---------------------------------------------------------------
    friend void read_binary(std::istream& is, hash_multimap& m) {
        m.read_binary_(is);
    }
    friend void write_binary(std::ostream& os, const hash_multimap& m) {
        m.write_binary_(os);
    }


private:
    /****************************************************************
     * @brief deserialize hashmap from input stream
     * @param is:        input stream
     * @param memconfig: callback that can be used to configurate allocators
     */
    void read_binary_(std::istream& is)
    {
        using std::uint64_t;

        uint64_t nkeys = 0;
        read_binary(is, nkeys);
        if(nkeys < 1) return;
        uint64_t nvalues = 0;
        read_binary(is, nvalues);
        if(nvalues < 1) return;

        clear();
        reserve_values(nvalues);
        reserve_keys(nkeys);

        for(uint64_t i = 0; i < nkeys; ++i) {
            key_type key;
            bucket_size_type nvals;
            read_binary(is, key);
            read_binary(is, nvals);

            auto it = find_slot_for_insertion(key);
            it->key(key);
            it->resize(nvals, alloc_);

            for(auto v = it->values_, e = v+nvals; v < e; ++v) {
                read_binary(is, *v);
            }
        }
        numKeys_ = nkeys;
        numValues_ = nvalues;
    }


    //---------------------------------------------------------------
    /**
     * @brief serialize hashmap to output stream
     */
    void write_binary_(std::ostream& os) const
    {
        using std::uint64_t;

        write_binary(os, uint64_t(key_count()));
        write_binary(os, uint64_t(value_count()));

        for(const auto& bucket : buckets_) {
            //store non-empty buckets only
            if(!bucket.empty()) {
                write_binary(os, bucket.key());
                write_binary(os, bucket.size());

                for(const auto& ref : bucket) {
                    write_binary(os, ref);
                }
            }
        }
    }


    //---------------------------------------------------------------
    iterator
    find_occupied_slot(const key_type& key)
    {
//        probelen_t probelen = 0;

        auto fst = buckets_.begin() + (hash_(key) % buckets_.size());
        auto it = fst;
        auto end = buckets_.end();

        //find bucket by linear probing
        while(it < end) {
            if(it->empty()) return buckets_.end();
            if(keyEqual_(it->key(), key)) return it;
            //early out, possible due to Robin Hood invariant
//            if(probelen < std::numeric_limits<probelen_t>::max()) {
//                if(probelen > it->probelen_) return buckets_.end();
//                ++probelen;
//            }
            ++it;
            if(it == buckets_.end()) {
                it = buckets_.begin();
                end = fst;
            }
        }

        return buckets_.end();
    }

    //-----------------------------------------------------
    iterator
    find_slot_for_insertion(key_type key)
    {
        auto fst = buckets_.begin() + (hash_(key) % buckets_.size());
        auto end = buckets_.end();
//        auto ins = end;
        auto it = fst;

        //find bucket by linear probing
//        bucket_type dummy{std::move(key)};
        while(it < end) {
            //new slot found
            if(it->empty()) {
//                *it = std::move(dummy);
//                return ins != buckets_.end() ? ins : it;
                it->key(std::move(key));
                return it;
            }
            //key already inserted
//            if(keyEqual_(it->key(), dummy.key())) {
            if(keyEqual_(it->key(), key)) {
                return it;
            }
            //Robin Hood criterion
//            if(dummy.probelen_ < std::numeric_limits<probelen_t>::max()) {
//                if(it->probelen_ < dummy.probelen_) {
                    //check if no insert position has been found yet
                    //after the first swap (*ins) will have the key to be inserted
                    //and we'll continue with (*it)'s old key
//                    if(ins == buckets_.end()) ins = it;
//                    std::swap(*it,dummy);
//                }
//                ++dummy.probelen_;
//            }
            ++it;
            if(it == buckets_.end()) {
                it = buckets_.begin();
                end = fst;
            }
        }
        return it;
    }

    //---------------------------------------------------------------
    size_type
    erase_slot(iterator del)
    {
        //delete element
        auto s = del->size();
        del->clear(alloc_);
        --numKeys_;
        numValues_ -= s;

        //shift backwards
        auto nxt = del;
        auto end = buckets_.end();
        while(nxt < end) {
            ++nxt;
            if(nxt == buckets_.end()) {
                nxt = buckets_.begin();
                end = del;
            }
//            if(nxt->empty() || nxt->probelen_ == 0) {
            if(nxt->empty()) {
                break;
            }
            else {
                std::swap(*del, *nxt);
//                --del->probelen_;
            }
            ++del;
            if(del == buckets_.end()) del = buckets_.begin();
        }

        return s;
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
    size_type numKeys_;
    size_type numValues_;
    float maxLoadFactor_;
    hasher hash_;
    key_equal keyEqual_;
    value_allocator alloc_;
    bucket_store_t buckets_;
};



} // namespace mc


#endif

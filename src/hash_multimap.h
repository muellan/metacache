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

#include "memory.h"


namespace mc {


/*****************************************************************************
 *
 * @brief   (integer) key -> value hashed multimap
 *          optimized for many values per key
 *
 * @details Uses open addressing with linear probing to resolve hash conflicts.
 *          The main differences to std::[unordered_][multi]map
 *           - iterators will iterate over *buckets* not over all *values*
 *           - iterators can be invalidated by inserts/rehashes
 *
 * @tparam  Key:
 * @tparam  ValueT:
 * @tparam  Hash:
 * @tparam  ValueAllocator:  controls allocation of values, *not* of buckets
 * @tparam  BucketAllocator: controls allocation of buckets, *not* of values
 *
 *****************************************************************************/
template<
    class Key, class ValueT,
    class Hash = std::hash<Key>,
    class ValueAllocator = std::allocator<ValueT>,
    class BucketAllocator = std::allocator<Key>
>
class hash_multimap
{
//    static_assert(std::is_pod<Key>::value,    "Key must be a POD type");
//    static_assert(std::is_pod<ValueT>::value, "Value must be a POD type");

    using value_alloc = std::allocator_traits<ValueAllocator>;

public:
    //---------------------------------------------------------------
    using value_type      = ValueT;
    using key_type        = Key;
    using hasher          = Hash;
    using value_allocator = ValueAllocator;

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
        using size_type        = std::uint16_t;
        using key_type         = hash_multimap::key_type;
        using value_type       = hash_multimap::value_type;
        using reference        = value_type&;
        using const_reference  = const value_type&;
        using iterator         = value_type*;
        using const_iterator   = const value_type*;

        bucket_type():
            values_{nullptr}, key_{}, size_(0), capacity_(0)
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


    private:
        //-----------------------------------------------------
        void insert(const value_type& v, value_allocator& alloc) {
            grow_by(1, alloc);
            values_[size_-1] = v;
        }
        //-------------------------------------------
        void insert(value_type&& v, value_allocator& alloc) {
            grow_by(1, alloc);
            values_[size_-1] = std::move(v);
        }

        //-----------------------------------------------------
        void shrink_to(size_type n, value_allocator&) noexcept {
            if(size_ > n) size_ = n;
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
        }
        //-------------------------------------------
        void replace(bucket_type&& src, value_allocator& alloc) {
            key_ = std::move(src.key_);
            clear(alloc);
            size_ = src.size_;
            capacity_ = src.capacity_;
            values_ = src.values_;
            src.values_ = nullptr;
            src.size_ = 0;
            src.capacity_ = 0;
        }

        //-----------------------------------------------------
        void key(const key_type& key) noexcept { key_ = key; }

        //-----------------------------------------------------
        void resize(size_type n, value_allocator& alloc) {
            if(size_ < n) {
                grow_by(n - size_, alloc);
            }
            else if(size_ > n) {
                size_ = n;
            }
        }

        //-----------------------------------------------------
        void grow_by(size_type n, value_allocator& alloc) {
            if(values_) {
                auto nsize = std::uint64_t(size_) + n;
                if(nsize > 65535) throw std::runtime_error{
                    "hash_multimap: bucket size exceeded maximum"};

                if((size_ + n) > capacity_) {
                    auto ncap = std::uint64_t(size_ * 1.3 + n);
                    if(ncap > 65535) ncap = 65535;
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
                    for(auto i = values_ + size_, e = i + n; i < e; ++i) {
                        value_alloc::construct(alloc, i);
                    }
                }
                size_ = nsize;
            }
            else {
                values_ = value_alloc::allocate(alloc, n);
                if(!std::is_pod<value_type>::value) {
                    for(auto i = values_, e = i + n; i < e; ++i) {
                        value_alloc::construct(alloc, i);
                    }
                }
                capacity_ = n;
                size_ = n;
            }
        }

        //-----------------------------------------------------
        value_type* values_;  //64 bits
        key_type key_;        //usually only 32 bits
        size_type size_;      //16 bits
        size_type capacity_;  //16 bits
    };


    //---------------------------------------------------------------
    using bucket_allocator = typename
        std::allocator_traits<BucketAllocator>::template rebind_alloc<bucket_type>;

private:
    //-----------------------------------------------------
    using store_t = std::vector<bucket_type,bucket_allocator>;


public:
    //-----------------------------------------------------
    using size_type        = typename store_t::size_type;
    using bucket_size_type = typename bucket_type::size_type;
    using reference        = bucket_type&;
    using const_reference  = const bucket_type&;
    //-----------------------------------------------------
    using iterator        = typename store_t::iterator;
    using const_iterator  = typename store_t::const_iterator;
    //-----------------------------------------------------
    using local_iterator       = typename bucket_type::iterator;
    using const_local_iterator = typename bucket_type::const_iterator;


    //---------------------------------------------------------------
    hash_multimap(const value_allocator& alloc = value_allocator{}):
        numKeys_(0), numValues_(0), maxLoadFactor_(0.80),
        hash_{}, alloc_{alloc},
        mutables_{},
        buckets_{}
    {
        buckets_.resize(5);
    }
    //-----------------------------------------------------
    hash_multimap(const hasher& hash,
                  const value_allocator& alloc = value_allocator{})
    :
        numKeys_(0), numValues_(0), maxLoadFactor_(0.80),
        hash_{hash}, alloc_{alloc},
        mutables_{},
        buckets_{}
    {
        buckets_.resize(5);
    }


    //-------------------------------------------------------------------
    hash_multimap(const hash_multimap&) = delete;

    //-----------------------------------------------------
    hash_multimap(hash_multimap&& src):
        numKeys_(src.numKeys_), numValues_(src.numValues_),
        maxLoadFactor_(src.maxLoadFactor_),
        hash_{std::move(src.hash_)},
        alloc_{std::move(src.alloc_)},
        mutables_{},
        buckets_{std::move(src.buckets_)}
    { }


    //-------------------------------------------------------------------
    hash_multimap& operator = (const hash_multimap&) = delete;

    //-----------------------------------------------------
    hash_multimap& operator = (hash_multimap&& src)
    {
        numKeys_ = src.numKeys_;
        numValues_ = src.numValues_;
        maxLoadFactor_ = src.maxLoadFactor_;
        hash_ = std::move(src.hash_);
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


    //---------------------------------------------------------------
    const value_allocator&
    get_allocator() const noexcept {
        return alloc_;
    }


    //---------------------------------------------------------------
    const hasher&
    hash_function() const noexcept {
        return hash_;
    }


    //---------------------------------------------------------------
    void rehash(size_type n)
    {
        if(n == buckets_.size()) return;

        std::lock_guard<std::mutex> lock(mutables_);

//        std::cout << "\nREHASH to " << n << std::endl;

        //make temporary new map with our old allocator
        hash_multimap newmap{alloc_};
        newmap.maxLoadFactor_ = maxLoadFactor_;
        //this might throw
        newmap.buckets_.resize(n);

        //move old bucket contents into new hash slots
        //this should use only noexcept operations
        for(auto& b : buckets_) {
            if(!b.empty()) {
                newmap.replace(std::move(b));
            }
        }

        //should all be noexcept
        buckets_ = std::move(newmap.buckets_);
        numKeys_ = newmap.numKeys_;
        numValues_ = newmap.numValues_;
        hash_ = newmap.hash_;
    }


    //---------------------------------------------------------------
    template<class Value>
    iterator
    insert(const key_type& key, Value&& value)
    {
        make_sure_enough_buckets_left(1);  //might lock

        //TODO make concurrency-safe

        auto it = find_slot(key);

        if(it->empty()) {
            it->key(key);
            it->insert(std::forward<Value>(value), alloc_);
            ++numKeys_;
            numValues_ += it->size();
            return it;
        }
        else { //key already inserted
            auto oldsize = it->size();
            it->insert(std::forward<Value>(value), alloc_);
            numValues_ += it->size() - oldsize;
            return it;
        }
    }


    //---------------------------------------------------------------
    size_type
    erase(const key_type& key)
    {
        //TODO make concurrency-safe

        auto it = find(key);
        return (it != end()) ? erase(it) : 0;
    }
    //-----------------------------------------------------
    size_type
    erase(const_iterator it)
    {
        //TODO make concurrency-safe

        auto s = it->size();
        const_cast<bucket_type*>(&(*it))->clear(alloc_);
        --numKeys_;
        numValues_ -= s;

        return s;
    }


    //---------------------------------------------------------------
    void limit(const key_type& key, bucket_size_type n)
    {
        //TODO make concurrency-safe

        auto it = find(key);
        if(it != end()) limit(it, n);
    }
    //-----------------------------------------------------
    void limit(const_iterator bit, bucket_size_type n)
    {
        //TODO make concurrency-safe

        auto so = bit->size();
        const_cast<bucket_type*>(&(*bit))->shrink_to(n, alloc_);
        auto sn = bit->size();
        if(sn < so) numValues_ -= so - sn;
    }


    //---------------------------------------------------------------
    void clear()
    {
        if(numKeys_ < 1) return;

        std::lock_guard<std::mutex> lock(mutables_);

        //free bucket memory
        for(auto& b : buckets_) {
            b.clear(alloc_);
        }
        buckets_.clear();
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
        auto it = find_slot(key);
        return it->empty() ? end() : it;
    }
    //-----------------------------------------------------
    const_iterator
    find(const key_type& key) const
    {
        return const_iterator(const_cast<hash_multimap*>(this)->find(key));
    }

    //-----------------------------------------------------
    size_type
    count(const key_type& key) const {
        auto it = find_slot(key);
        return (!it->empty()) ? it->size() : 0;
    }


    //---------------------------------------------------------------
    float load_factor() const noexcept {
        return (numKeys_ / float(buckets_.size()));
    }

    //---------------------------------------------------------------
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


    /****************************************************************
     * @brief deserialize hashmap from input stream
     * @param is:        input stream
     * @param memconfig: callback that can be used to configurate allocators
     */
    template<class MemoryConfigurator>
    void read(std::istream& is,
              MemoryConfigurator&& memconfig = [](std::size_t,std::size_t){})
    {
        static_assert(std::is_pod<ValueT>::value,
                      "hash_multimap::read requires ValueT to be a POD type.");

        size_t nkeys = 0;
        is.read(reinterpret_cast<char*>(&nkeys), sizeof(nkeys));
        if(nkeys < 1) return;
        size_t nvalues = 0;
        is.read(reinterpret_cast<char*>(&nvalues), sizeof(nvalues));
        if(nvalues < 1) return;

        memconfig(nkeys, nvalues);

        //locks
        clear();
        rehash(typename store_t::size_type(10 + (1/max_load_factor() * nkeys)));

        std::lock_guard<std::mutex> lock(mutables_);

        for(size_t i = 0; i < nkeys; ++i) {
            key_type key;
            bucket_size_type nvals;
            is.read(reinterpret_cast<char*>(&key), sizeof(key));
            is.read(reinterpret_cast<char*>(&nvals), sizeof(nvals));

            auto it = find_slot(key);
            it->key(key);
            it->resize(nvals, alloc_);

            for(auto v = it->values_, e = v+nvals; v < e; ++v) {
                is.read(reinterpret_cast<char*>(&(*v)), sizeof(*v));
            }
        }
        numKeys_ = nkeys;
        numValues_ = nvalues;

//        std::cout << '\n'
//            << "keys:        " << key_count() << " <> " << nkeys << '\n'
//            << "values:      " << value_count() << " <> " << nvalues << '\n'
//            << "buckets:     " << bucket_count() << '\n'
//            << "bucket size: " << (8*sizeof(bucket(0))) << '\n'
//            << "key size:    " << (8*sizeof(bucket(0).key())) << '\n'
//            << "value size:  " << (8*sizeof(bucket(0)[0])) << '\n'
//            << "gid size:    " << (8*sizeof(bucket(0)[0].gid)) << '\n'
//            << "win size:    " << (8*sizeof(bucket(0)[0].win)) << '\n'
//            << "addr0      : " << &bucket(0)[0] << '\n'
//            << "addr1      : " << &bucket(0)[1] << '\n'
//            << std::flush;
//        char c;
//        std::cin >> c;

//        for(const auto& b : buckets_) {
//            std::cout << b.size() << '\n';
//        }
//        std::cin >> c;

    }


    //---------------------------------------------------------------
    /**
     * @brief serialize hashmap to output stream
     */
//    template<class = typename std::enable_if<std::is_pod<value_type>::value>::type>
    void write(std::ostream& os) const
    {
        std::lock_guard<std::mutex> lock(mutables_);

        static_assert(std::is_pod<ValueT>::value,
                      "hash_multimap::write requires ValueT to be a POD type.");

        //write sketch values (keys) and their origin references
        auto nkeys = key_count();
        os.write(reinterpret_cast<char*>(&nkeys), sizeof(nkeys));
        auto nvalues = value_count();
        os.write(reinterpret_cast<char*>(&nvalues), sizeof(nvalues));

        for(const auto& bucket : buckets_) {
            //store non-empty buckets only
            if(!bucket.empty()) {
                auto key = bucket.key();
                os.write(reinterpret_cast<const char*>(&key), sizeof(key));

                auto nvals = bucket.size();
                os.write(reinterpret_cast<const char*>(&nvals), sizeof(nvals));

                for(const auto& ref : bucket) {
                    os.write(reinterpret_cast<const char*>(&ref), sizeof(ref));
                }
            }
        }
    }


private:

    //---------------------------------------------------------------
    iterator
    replace(bucket_type&& b)
    {
        make_sure_enough_buckets_left(1);  //might lock

        //TODO make concurrency-safe

        auto it = find_slot(b.key());

        if(it->empty()) {
            it->replace(std::move(b), alloc_);
            ++numKeys_;
            numValues_ += it->size();
            return it;
        }
        else { //key already inserted
            auto oldsize = it->size();
            it->replace(std::move(b), alloc_);
            numValues_ += it->size() - oldsize;
            return it;
        }
    }

    //---------------------------------------------------------------
    iterator
    find_slot(const key_type& key)
    {
        //initial bucket index
        auto it = buckets_.begin() + (hash_(key) % buckets_.size());

//        char i = 0;
//        for(; i < 4 && it < buckets_.end(); ++it, ++i) {
//            if(it->empty() || it->key() == key) return it;
//        }
//        for() {
//
//        }

        //find correct bucket by linear probing
        //start at current bucket
        for(; it < buckets_.end(); ++it) {
            if(it->empty() || it->key() == key) return it;
        }
        //start at first bucket
        auto oldit = it;
        it = buckets_.begin();
        for(; it < oldit; ++it) {
            if(it->empty() || it->key() == key) return it;
        }

        return it;
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
    value_allocator alloc_;
    mutable std::mutex mutables_;
    store_t buckets_;
};



} // namespace mc


#endif

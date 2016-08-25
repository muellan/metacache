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
 *****************************************************************************/
template<
    class Key, class ValueT,
    class Hash = std::hash<Key>
>
class hash_multimap
{
public:
    //---------------------------------------------------------------
    using value_type = ValueT;
    using key_type   = Key;
    using hasher     = Hash;

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

        bool empty() const noexcept { return !values_.get(); }

        const key_type& key() const noexcept { return key_; }

        reference
        operator [](size_type i) noexcept { return values_[i]; }

        const_reference
        operator [](size_type i) const noexcept { return values_[i]; }

        //-------------------------------------------
              iterator  begin()       noexcept { return values_.get(); }
        const_iterator  begin() const noexcept { return values_.get(); }
        const_iterator cbegin() const noexcept { return values_.get(); }

              iterator  end()       noexcept { return values_.get() + size_; }
        const_iterator  end() const noexcept { return values_.get() + size_; }
        const_iterator cend() const noexcept { return values_.get() + size_; }


    private:
        //-----------------------------------------------------
        void insert(const value_type& v) {
            grow_by(1);
            values_[size_-1] = v;
        }
        //-------------------------------------------
        void insert(value_type&& v) {
            grow_by(1);
            values_[size_-1] = std::move(v);
        }

        //-----------------------------------------------------
        void shrink_to(size_type n) noexcept {
            if(size_ > n) size_ = n;
        }

        //-------------------------------------------
        void clear() {
            values_.reset(nullptr);
            size_ = 0;
            capacity_ = 0;
        }
        //-------------------------------------------
        void replace(size_type n, std::unique_ptr<value_type[]>&& values) {
            size_ = n;
            capacity_ = n;
            values_ = std::move(values);
        }

        //-----------------------------------------------------
        void key(const key_type& key) noexcept { key_ = key; }

        //-----------------------------------------------------
        void grow_by(size_type n) {
            if(values_) {
                auto nsize = std::uint64_t(size_) + n;
                if(nsize >  65535) throw std::runtime_error{
                        "hash_multimap: bucket size exceeds maximum"};

                if((size_ + n) > capacity_) {
                    auto ncap = std::uint64_t(size_ * 1.3 + n);
                    if(ncap > 65535) ncap = 65535;

                    auto nvals = std::unique_ptr<value_type[]>{
                                     new value_type[ncap]};

                    auto nv = nvals.get();
                    auto v = values_.get();
                    for(size_type i = 0; i < size_; ++i, ++v, ++nv) {
                        *nv = std::move(*v);
                    }
                    values_ = std::move(nvals);
                    capacity_ = size_type(ncap);
                }
                size_ = nsize;
            }
            else {
                values_.reset(new value_type[n]);
                capacity_ = n;
                size_ = n;
            }
        }
        //-------------------------------------------
        std::unique_ptr<value_type[]> values_; //64 bits
        key_type key_;        //usually only 32 bits
        size_type size_;      //16 bits
        size_type capacity_;  //16 bits
    };


private:
    //-----------------------------------------------------
    using store_t = std::vector<bucket_type>;


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
    using local_iterator        = typename bucket_type::iterator;
    using const_local_iterator  = typename bucket_type::const_iterator;


    //---------------------------------------------------------------
    hash_multimap():
        numKeys_(0), numValues_(0), maxLoadFactor_(0.80),
        hash_{},
//        mutables_{},
        buckets_{}
    {
        buckets_.resize(5);
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
    /**
     * @brief reserve memory for 'n' values in total
     */
    void reserve(size_type n)
    {
        if(n == capacity()) return;

//        std::lock_guard<std::mutex> lock(mutables_);

        //TODO
    }
    //-----------------------------------------------------
    size_type capacity() const noexcept {
        //TODO
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
    void rehash(size_type n)
    {
        if(n == buckets_.size()) return;

//        std::lock_guard<std::mutex> lock(mutables_);

//        std::cout << "\nREHASH to " << n << std::endl;

        //make temporary new map
        hash_multimap newmap;
        newmap.maxLoadFactor_ = maxLoadFactor_;
        //this might throw, but we haven't modified anything yet
        newmap.buckets_.resize(n);

        //move old bucket contents into new hash slots
        //this should use only noexcept operations
        for(auto& b : buckets_) {
            if(!b.empty()) {
                newmap.replace(std::move(b.key_), b.size_, std::move(b.values_));
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
        make_sure_enough_buckets_left(1);

        auto it = find_slot(key);

        if(it->empty()) {
            it->key(key);
            it->insert(std::forward<Value>(value));
            ++numKeys_;
            numValues_ += it->size();
            return it;
        }
        else { //key already inserted
            auto oldsize = it->size();
            it->insert(std::forward<Value>(value));
            numValues_ += it->size() - oldsize;
            return it;
        }
    }

    //-----------------------------------------------------
    template<class Values>
    iterator
    replace(const key_type& key, size_type size, Values&& values)
    {
        make_sure_enough_buckets_left(1);

        auto it = find_slot(key);

        if(it->empty()) {
            it->key(key);
            it->replace(size, std::forward<Values>(values));
            ++numKeys_;
            numValues_ += it->size();
            return it;
        }
        else { //key already inserted
            auto oldsize = it->size();
            it->replace(size, std::forward<Values>(values));
            numValues_ += it->size() - oldsize;
            return it;
        }
    }


    //---------------------------------------------------------------
    size_type
    erase(const key_type& key)
    {
        auto it = find(key);
        return (it != end()) ? erase(it) : 0;
    }
    //-----------------------------------------------------
    size_type
    erase(const_iterator it)
    {
        auto s = it->size();
        const_cast<bucket_type*>(&(*it))->clear();
        --numKeys_;
        numValues_ -= s;

        return s;
    }


    //---------------------------------------------------------------
    void limit(const key_type& key, bucket_size_type n) {
        auto it = find(key);
        if(it != end()) limit(it, n);
    }
    //-----------------------------------------------------
    void limit(const_iterator bit, bucket_size_type n) {
        auto so = bit->size();
        const_cast<bucket_type*>(&(*bit))->shrink_to(n);
        auto sn = bit->size();
        if(sn < so) numValues_ -= so - sn;
    }


    //---------------------------------------------------------------
    void clear() {
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


    //---------------------------------------------------------------
    /**
     * @brief deserialize hashmap from input stream
     */
//    template<class = typename std::enable_if<std::is_pod<value_type>::value>::type>
    void read(std::istream& is)
    {
        static_assert(std::is_pod<ValueT>::value,
                      "hash_multimap::read requires ValueT to be a POD type.");

        size_t nkeys = 0;
        is.read(reinterpret_cast<char*>(&nkeys), sizeof(nkeys));
        if(nkeys < 1) return;
        size_t nvalues = 0;
        is.read(reinterpret_cast<char*>(&nvalues), sizeof(nvalues));
        if(nvalues < 1) return;

        clear();
        rehash( typename store_t::size_type(10 + (1/max_load_factor() * nkeys)));

        for(size_t i = 0; i < nkeys; ++i) {
            key_type key;
            bucket_size_type nvals;
            is.read(reinterpret_cast<char*>(&key), sizeof(key));
            is.read(reinterpret_cast<char*>(&nvals), sizeof(nvals));

            auto vals = std::unique_ptr<value_type[]>{new value_type[nvals]};

            for(auto v = vals.get(), e = v+nvals; v < e; ++v) {
                is.read(reinterpret_cast<char*>(&(*v)), sizeof(*v));
            }

            replace(key, nvals, std::move(vals));
        }

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
    }


    //---------------------------------------------------------------
    /**
     * @brief serialize hashmap to output stream
     */
//    template<class = typename std::enable_if<std::is_pod<value_type>::value>::type>
    void write(std::ostream& os) const
    {
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
//    mutable std::mutex mutables_;
    store_t buckets_;
};



} // namespace mc


#endif

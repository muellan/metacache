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

#ifndef MC_SKETCHERS_H_
#define MC_SKETCHERS_H_


#include <array>
#include <iterator>
#include <algorithm>
#include <limits>
//#include <random>

#include "dna_encoding.h"
#include "hash_family.h"


namespace mc {


/*****************************************************************************
 *
 *
 *****************************************************************************/
class single_function_min_hasher
{
public:
    //---------------------------------------------------------------
    using kmer_type = std::uint32_t;
    using hash_type = std::uint32_t;
    using result_type = std::vector<hash_type>;
    //-----------------------------------------------------
    using kmer_size_type   = numk_t;
    using sketch_size_type = typename result_type::size_type;


    //---------------------------------------------------------------
    static constexpr std::uint8_t max_kmer_size() noexcept {
        return max_word_size<kmer_type,2>::value;
    }
    static constexpr sketch_size_type max_sketch_size() noexcept {
        return std::numeric_limits<sketch_size_type>::max();
    }


    //---------------------------------------------------------------
    single_function_min_hasher():
        kmerSize_(16), maxSketchSize_(32)
    {}


    //---------------------------------------------------------------
    kmer_size_type
    kmer_size() const noexcept {
        return kmerSize_;
    }
    void
    kmer_size(kmer_size_type k) noexcept {
        if(k < 1) k = 1;
        if(k > max_kmer_size()) k = max_kmer_size();
        kmerSize_ = k;
    }

    //---------------------------------------------------------------
    sketch_size_type
    sketch_size() const noexcept {
        return maxSketchSize_;
    }
    void
    sketch_size(sketch_size_type s) noexcept {
        if(s < 1) s = 1;
        maxSketchSize_ = s;
    }


    //---------------------------------------------------------------
    template<class Sequence>
    result_type
    operator () (const Sequence& s) const {
        using std::begin;
        using std::end;
        return operator()(begin(s), end(s));
    }

    //-----------------------------------------------------
    template<class InputIterator>
    result_type
    operator () (InputIterator first, InputIterator last) const
    {
        using std::distance;
        using std::begin;
        using std::end;

        const auto numkmers = sketch_size_type(distance(first, last) - kmerSize_ + 1);
        const auto sketchsize = std::min(maxSketchSize_, numkmers);

        auto sketch = result_type(sketchsize, hash_type(~0));

        for_each_kmer_2bit<kmer_type>(kmerSize_, first, last,
            [&] (kmer_type kmer, half_size_t<kmer_type> ambiguous) {
                if(!ambiguous) {
                    //make canonical (strand neutral) k-mer and hash it
                    auto h = thomas_mueller_hash(make_canonical(kmer, kmerSize_));

                    if(h < sketch.back()) {
                        auto pos = std::upper_bound(sketch.begin(), sketch.end(), h);
                        if(pos != sketch.end()) {
                            sketch.pop_back();
                            sketch.insert(pos, h);
                        }
                    }
                }
            });

//        for(auto x : sketch) std::cout << x << ' '; std::cout << '\n';

        return sketch;
    }

private:
    //---------------------------------------------------------------
    kmer_size_type kmerSize_;
    sketch_size_type maxSketchSize_;

    //---------------------------------------------------------------
    static std::uint32_t
    thomas_mueller_hash(std::uint32_t x) noexcept {
        x = ((x >> 16) ^ x) * 0x45d9f3b;
        x = ((x >> 16) ^ x) * 0x45d9f3b;
        x = ((x >> 16) ^ x);
        return x;
    }
    static void thomas_mueller_hash(std::uint64_t) = delete;
    static void thomas_mueller_hash(std::uint16_t) = delete;
    static void thomas_mueller_hash(std::uint8_t) = delete;

    //---------------------------------------------------------------
    static std::uint32_t
    murmur3_int_init(std::uint32_t x) noexcept {
        x ^= x >> 16;
        x *= 0x85ebca6b;
        x ^= x >> 13;
        x *= 0xc2b2ae35;
        x ^= x >> 16;
        return x;
    }
    static void murmur3_int_init(std::uint64_t) = delete;
    static void murmur3_int_init(std::uint16_t) = delete;
    static void murmur3_int_init(std::uint8_t) = delete;

    //-----------------------------------------------------
    static std::uint32_t
    nvidia_hash(std::uint32_t x) noexcept {
        x = (x + 0x7ed55d16) + (x << 12);
        x = (x ^ 0xc761c23c) ^ (x >> 19);
        x = (x + 0x165667b1) + (x <<  5);
        x = (x + 0xd3a2646c) ^ (x <<  9);
        x = (x + 0xfd7046c5) + (x <<  3);
        x = (x ^ 0xb55a4f09) ^ (x >> 16);
        return x;
    }
    static void nvidia_hash(std::uint64_t) = delete;
    static void nvidia_hash(std::uint16_t) = delete;
    static void nvidia_hash(std::uint8_t) = delete;

    //-----------------------------------------------------
    static std::uint32_t
    hash_64_to_32(std::uint64_t x) noexcept {
        x = (~x) + (x << 18); // x = (x << 18) - x - 1;
        x = x ^ (x >> 31);
        x = x * 21; // x = (x + (x << 2)) + (x << 4);
        x = x ^ (x >> 11);
        x = x + (x << 6);
        x = x ^ (x >> 22);
        return std::uint32_t(x);
    }
//    std::hash<hash_t> hash_;
};




/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
class multi_function_min_hasher
{

public:
    //---------------------------------------------------------------
    using kmer_type = std::uint32_t;
    using hash_type = std::uint64_t;
    using result_type = std::vector<hash_type>;
    //-----------------------------------------------------
    using kmer_size_type   = numk_t;
    using sketch_size_type = typename result_type::size_type;


    //---------------------------------------------------------------
    static constexpr std::uint8_t max_kmer_size() noexcept {
        return max_word_size<kmer_type,2>::value;
    }
    static constexpr sketch_size_type max_sketch_size() noexcept {
        return std::numeric_limits<sketch_size_type>::max();
    }


    //---------------------------------------------------------------
    multi_function_min_hasher():
        kmerSize_(16), maxSketchSize_(32), hash_{}
    {}


    //---------------------------------------------------------------
    kmer_size_type
    kmer_size() const noexcept {
        return kmerSize_;
    }
    void
    kmer_size(kmer_size_type k) noexcept {
        if(k < 1) k = 1;
        if(k > max_kmer_size()) k = max_kmer_size();
        kmerSize_ = k;
    }

    //---------------------------------------------------------------
    sketch_size_type
    sketch_size() const noexcept {
        return maxSketchSize_;
    }
    void
    sketch_size(sketch_size_type s) noexcept {
        if(s < 1) s = 1;
        if(s > 128) s = 128;
        maxSketchSize_ = s;
    }


    //---------------------------------------------------------------
    template<class Sequence>
    result_type
    operator () (const Sequence& s) const {
        using std::begin;
        using std::end;
        return operator()(begin(s), end(s));
    }

    //-----------------------------------------------------
    template<class InputIterator>
    result_type
    operator () (InputIterator first, InputIterator last) const
    {
        using std::distance;
        using std::begin;
        using std::end;

        const auto numkmers = sketch_size_type(distance(first, last) - kmerSize_ + 1);
        const auto sketchsize = std::min(maxSketchSize_, numkmers);

        auto sketch = result_type(sketchsize, hash_type(~0));

        //least significant 32 bits of sketch values = kmer hash
        //most significant bits of sketch values = hash function index
        //=> sketch values of different hash functions won't collide
        for_each_kmer_2bit<kmer_type>(kmerSize_, first, last,
            [&] (kmer_type kmer, half_size_t<kmer_type> ambiguous) {
                if(!ambiguous) {
                    //make canonical (strand neutral) k-mer and hash it
                    auto cankmer = make_canonical(kmer, kmerSize_);
                    //for each function keep track of smallest hash value
                    for(std::size_t i = 0; i < sketch.size(); ++i) {
                        auto h = hash_(i, cankmer);
                        if(h < sketch[i]) sketch[i] = h;
                    }
                }
            });

        //set function indices
        for(std::size_t i = 0; i < sketch.size(); ++i) {
            sketch[i] += std::uint64_t(i) << 32;
        }

//        for(auto x : sketch) std::cout << x << ' '; std::cout << '\n';

        return sketch;
    }

private:
    //---------------------------------------------------------------
    kmer_size_type kmerSize_;
    sketch_size_type maxSketchSize_;
    hash32_family128 hash_;
};



} // namespace mc


#endif

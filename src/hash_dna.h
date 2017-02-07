/*****************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
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


#include <algorithm>
#include <array>
#include <cstdint>
#include <iterator>
#include <limits>
#include <utility>
#include <type_traits>
//#include <random>

#include "dna_encoding.h"
#include "hash_int.h"
#include "io_serialize.h"


namespace mc {


/*****************************************************************************
 *
 * @brief default min-hasher that uses the 'sketch_size' lexicographically
 *
 *        smallest hash values of *one* hash function
 *
 *****************************************************************************/
template<class KmerT, class Hash = same_size_hash<KmerT>>
class single_function_min_hasher
{
public:
    //---------------------------------------------------------------
    using kmer_type    = KmerT;
    using hasher       = Hash;
    using feature_type = typename std::result_of<hasher(kmer_type)>::type;
    using sketch_type  = std::vector<feature_type>;
    //-----------------------------------------------------
    using kmer_size_type   = numk_t;
    using sketch_size_type = typename sketch_type::size_type;


    //---------------------------------------------------------------
    static constexpr std::uint8_t max_kmer_size() noexcept {
        return max_word_size<kmer_type,2>::value;
    }
    static constexpr sketch_size_type max_sketch_size() noexcept {
        return std::numeric_limits<sketch_size_type>::max();
    }


    //---------------------------------------------------------------
    explicit
    single_function_min_hasher(hasher hash = hasher{}):
        hash_(std::move(hash)), k_(16), sketchSize_(16)
    {}


    //---------------------------------------------------------------
    kmer_size_type
    kmer_size() const noexcept {
        return k_;
    }
    void
    kmer_size(kmer_size_type k) noexcept {
        if(k < 1) k = 1;
        if(k > max_kmer_size()) k = max_kmer_size();
        k_ = k;
    }

    //---------------------------------------------------------------
    sketch_size_type
    sketch_size() const noexcept {
        return sketchSize_;
    }
    void
    sketch_size(sketch_size_type s) noexcept {
        if(s < 1) s = 1;
        sketchSize_ = s;
    }


    //---------------------------------------------------------------
    template<class Sequence>
    sketch_type
    operator () (const Sequence& s) const {
        using std::begin;
        using std::end;
        return operator()(begin(s), end(s));
    }

    //-----------------------------------------------------
    template<class InputIterator>
    sketch_type
    operator () (InputIterator first, InputIterator last) const
    {
        using std::distance;
        using std::begin;
        using std::end;

        const auto n = distance(first,last);
        if(n < k_) return sketch_type{};
        const auto s = std::min(sketchSize_, sketch_size_type(n - k_ + 1));
        if(s < 1) return sketch_type{};

        auto sketch = sketch_type(s, feature_type(~0));

        for_each_unambiguous_canonical_kmer_2bit<kmer_type>(k_, first, last,
            [&] (kmer_type kmer) {
                auto h = hash_(kmer);
                if(h < sketch.back()) {
                    auto pos = std::upper_bound(sketch.begin(), sketch.end(), h);
                    if(pos != sketch.end()) {
                        sketch.pop_back();
                        sketch.insert(pos, h);
                    }
                }
            });

        //check if some features are invalid (in case of many ambiguous kmers)
        if(!sketch.empty() && sketch.back() == feature_type(~0)) {
            for(auto i = sketch.begin(), e = sketch.end(); i != e; ++i) {
                if(*i == feature_type(~0)) {
                    sketch.erase(i,sketch.end());
                    break;
                }
            }
        }

        return sketch;
    }


    //---------------------------------------------------------------
    friend void
    write_binary(std::ostream& os, const single_function_min_hasher& h)
    {
        write_binary(os, std::uint64_t(h.k_));
        write_binary(os, std::uint64_t(h.sketchSize_));
    }

    //---------------------------------------------------------------
    friend void
    read_binary(std::istream& is, single_function_min_hasher& h)
    {
        std::uint64_t n = 0;
        read_binary(is, n);
        h.k_ = (n <= max_kmer_size()) ? n : max_kmer_size();

        n = 0;
        read_binary(is, n);
        h.sketchSize_ = (n <= max_sketch_size()) ? n : max_sketch_size();
    }


private:
    //---------------------------------------------------------------
    hasher hash_;
    kmer_size_type k_;
    sketch_size_type sketchSize_;
};


} //namespace mc

#endif

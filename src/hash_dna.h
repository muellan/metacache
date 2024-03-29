/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2024 André Müller (muellan@uni-mainz.de)
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

#ifndef MC_SKETCHERS_H_
#define MC_SKETCHERS_H_


#include "dna_encoding.h"
#include "hash_int.h"
#include "io_serialize.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <iterator>
#include <limits>
#include <type_traits>
#include <utility>


namespace mc {

/*************************************************************************//**
* @brief loops through all subranges of [first,last) of lenght 'len' with a
*        stride of 'stride'
*
* @param first    iterator to start of range
* @param last     iterator to one after end of range
* @param len      length of windows (subranges)
* @param stride   distance between first elements of two subsequent windows
* @param consume  function object / lambda consuming one window
*****************************************************************************/
template<class InputIterator, class Consumer>
inline void
for_each_window(InputIterator first, InputIterator last,
                const size_t len, const size_t stride,
                Consumer&& consume)
{
    using std::distance;
    // sequence not longer than window?
    if (size_t(distance(first,last)) <= len) {
        consume(first,last);
    }
    else {
        for (auto wend = first + len;
            wend <= last;
            first += stride, wend += stride)
        {
            consume(first, wend);
        }
        if (first < last) consume(first, last);
    }
}
//-------------------------------------------------------------------
template<class Sequence, class Consumer>
inline void
for_each_window(const Sequence& s,
                const size_t len, const size_t stride,
                Consumer&& consume)
{
    using std::begin;
    using std::end;
    for_each_window(begin(s), end(s), len, stride,
                    std::forward<Consumer>(consume));
}



/*************************************************************************//**
 * @brief sequence sketching parameters
 *****************************************************************************/
template<class KmerT>
struct sketching_options
{
    using kmer_type        = KmerT;
    //---------------------------------------------------------------
    using kmer_size_type   = numk_t;
    using sketch_size_type = std::uint32_t;
    using window_size_type = std::uint32_t;
    //---------------------------------------------------------------
    static constexpr std::uint8_t max_kmer_size() noexcept {
        return max_word_size<kmer_type,2>::value;
    }
    static constexpr sketch_size_type max_sketch_size() noexcept {
        return std::numeric_limits<sketch_size_type>::max();
    }
    static constexpr window_size_type max_window_size() noexcept {
        return std::numeric_limits<window_size_type>::max();
    }
    static constexpr window_size_type max_window_stride() noexcept {
        return std::numeric_limits<window_size_type>::max();
    }


    //---------------------------------------------------------------
    friend void
    write_binary(std::ostream& os, const sketching_options<kmer_type>& h)
    {
        write_binary(os, std::uint64_t(h.kmerlen));
        write_binary(os, std::uint64_t(h.sketchlen));
        write_binary(os, std::uint64_t(h.winlen));
        write_binary(os, std::uint64_t(h.winstride));
    }

    //---------------------------------------------------------------
    friend void
    read_binary(std::istream& is, sketching_options<kmer_type>& h)
    {
        std::uint64_t n = 0;
        read_binary(is, n);
        h.kmerlen = (n <= max_kmer_size()) ? n : max_kmer_size();

        n = 0;
        read_binary(is, n);
        h.sketchlen = (n <= max_sketch_size()) ? n : max_sketch_size();

        n = 0;
        read_binary(is, n);
        h.winlen = (n <= max_window_size()) ? n : max_window_size();

        n = 0;
        read_binary(is, n);
        h.winstride = (n <= max_window_size()) ? n : max_window_size();
    }


    //---------------------------------------------------------------
    // number of characters in substring
    kmer_size_type kmerlen;
    // number of features (kmer hashes) in a sketch of ONE window
    sketch_size_type sketchlen;
    // number of characters in one window
    window_size_type winlen;
    // difference between two successive window start positions
    window_size_type winstride;
};



/*************************************************************************//**
 *
 * @brief default min-hasher that uses the 'sketch_size' lexicographically
 *
 *        smallest *unique* hash values of *one* hash function
 *
 *****************************************************************************/
template<class KmerT, class Hash = same_size_hash<KmerT>>
class single_function_unique_min_hasher
{
public:
    //---------------------------------------------------------------
    using kmer_type    = KmerT;
    using hasher       = Hash;
    using feature_type = typename std::result_of<hasher(kmer_type)>::type;
    //---------------------------------------------------------------
    using sketch_type      = std::vector<feature_type>;


    //---------------------------------------------------------------
    explicit
    single_function_unique_min_hasher(hasher hash = hasher{}):
        sketch{},
        hash_(std::move(hash))
    {}


    //---------------------------------------------------------------
    template<class Sequence, class Consumer>
    void
    for_each_sketch(const Sequence& s,
                    const sketching_options<kmer_type>& opt,
                    Consumer&& consume)
    {
        using std::begin;
        using std::end;
        for_each_sketch(begin(s), end(s), opt,
                        std::forward<Consumer>(consume));
    }

    //-----------------------------------------------------
    template<class InputIterator, class Consumer>
    void
    for_each_sketch(InputIterator first, InputIterator last,
                    const sketching_options<kmer_type>& opt,
                    Consumer&& consume)
    {
        using sketch_size_type = typename sketching_options<kmer_type>::sketch_size_type;

        for_each_window(first, last, opt.winlen, opt.winstride,
            [&] (InputIterator first, InputIterator last) {
                using std::distance;
                using std::begin;
                using std::end;

                const auto n = distance(first,last);
                if (n < opt.kmerlen) return;

                const auto s = std::min(opt.sketchlen, sketch_size_type(n - opt.kmerlen + 1));
                if (s < 1) return;

                sketch.clear();
                sketch.resize(s, feature_type(~0));

                for_each_unambiguous_canonical_kmer_2bit<kmer_type>(opt.kmerlen, first, last,
                    [&] (kmer_type kmer) {
                        auto h = hash_(kmer);
                        if (h < sketch.back()) {
                            auto pos = std::lower_bound(sketch.begin(), sketch.end(), h);
                            // make sure we don't insert the same feature more than once
                            if (pos != sketch.end() && *pos != h) {
                                sketch.pop_back();
                                sketch.insert(pos, h);
                            }
                        }
                    });

                // check if some features are invalid (in case of many ambiguous kmers)
                if (!sketch.empty() && sketch.back() == feature_type(~0)) {
                    for (auto i = sketch.begin(), e = sketch.end(); i != e; ++i) {
                        if (*i == feature_type(~0)) {
                            sketch.erase(i,sketch.end());
                            break;
                        }
                    }
                }

                consume(sketch);
            });
    }


private:
    //---------------------------------------------------------------
    sketch_type sketch;
    hasher hash_;
};


} // namespace mc

#endif

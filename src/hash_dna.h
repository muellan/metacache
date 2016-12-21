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
#include "hash_family.h"
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
    using result_type  = std::vector<feature_type>;
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

        const auto n = distance(first,last);
        if(n < k_) return result_type{};
        const auto s = std::min(sketchSize_, sketch_size_type(n - k_ + 1));

        auto sketch = result_type(s, feature_type(~0));

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




/*****************************************************************************
 *
 * @brief default min-hasher that uses the 'sketch_size' lexicographically
 *
 *        smallest hash values of *one* hash function
 *
 *****************************************************************************/
class single_function_multi_kmer_min_hasher
{
public:
    //---------------------------------------------------------------
    using kmer_type    = std::uint64_t;
    using hasher       = same_size_hash<kmer_type>;
    using feature_type = typename std::result_of<hasher(kmer_type)>::type;
    using result_type  = std::vector<feature_type>;
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
    explicit
    single_function_multi_kmer_min_hasher():
        hash_(), sketchSize_(16)
    {}


    //---------------------------------------------------------------
    kmer_size_type
    kmer_size() const noexcept { return 32; }
    void kmer_size(kmer_size_type) noexcept { }

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

        const auto n = distance(first,last);
        if(n < 32) return result_type{};
        const auto s = std::min(sketchSize_, sketch_size_type(n - 31));

        auto sketch = result_type(s, feature_type(~0));
        for(kmer_size_type k = 16; k <= 32; k += 4) {
            for_each_unambiguous_canonical_kmer_2bit<kmer_type>(k, first, last,
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
        }

        return sketch;
    }


    //---------------------------------------------------------------
    friend void
    write_binary(std::ostream& os, const single_function_multi_kmer_min_hasher& h)
    {
        write_binary(os, std::uint64_t(h.sketchSize_));
    }

    //---------------------------------------------------------------
    friend void
    read_binary(std::istream& is, single_function_multi_kmer_min_hasher& h)
    {
        std::uint64_t n = 0;
        read_binary(is, n);
        h.sketchSize_ = (n <= max_sketch_size()) ? n : max_sketch_size();
    }


private:
    //---------------------------------------------------------------
    template<class InputIterator>
    result_type make_sketch(sketch_size_type s,
                            kmer_size_type k,
                            InputIterator first, InputIterator last) const
    {
        auto sketch = result_type(s, feature_type(~0));

        for_each_unambiguous_canonical_kmer_2bit<kmer_type>(k, first, last,
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

        return sketch;
    }

    //---------------------------------------------------------------
    hasher hash_;
    sketch_size_type sketchSize_;
};




/*****************************************************************************
 *
 * @brief min-hasher that uses a different hash function for each feature
 *
 *****************************************************************************/
class multi_function_min_hasher
{

public:
    //---------------------------------------------------------------
    using kmer_type    = std::uint32_t;
    using feature_type = std::uint64_t;
    using result_type  = std::vector<feature_type>;
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
        k_(16), sketchSize_(32), hash_{}
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
        if(s > 128) s = 128;
        sketchSize_ = s;
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

        const auto n = distance(first,last);
        if(n < k_) return result_type{};
        const auto s = std::min(sketchSize_, sketch_size_type(n - k_ + 1));

        auto sketch = result_type(s, feature_type(~0));

        //least significant 32 bits of features = kmer hash
        //most significant bits of features = hash function index
        //=> features of different hash functions won't collide
        for_each_unambiguous_canonical_kmer_2bit<kmer_type>(k_, first, last,
            [&] (kmer_type kmer) {
                //for each function keep track of smallest hash value
                for(std::size_t i = 0; i < sketch.size(); ++i) {
                    auto h = hash_(i, kmer);
                    if(h < sketch[i]) sketch[i] = h;
                }
            });

        //set function indices
        for(std::size_t i = 0; i < sketch.size(); ++i) {
            sketch[i] += std::uint64_t(i) << 32;
        }

//        for(auto x : sketch) std::cout << x << ' '; std::cout << '\n';

        return sketch;
    }


    //---------------------------------------------------------------
    friend void
    write_binary(std::ostream& os, const multi_function_min_hasher& h)
    {
        write_binary(os, std::uint64_t(h.k_));
        write_binary(os, std::uint64_t(h.sketchSize_));
    }

    //---------------------------------------------------------------
    friend void
    read_binary(std::istream& is, multi_function_min_hasher& h)
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
    kmer_size_type k_;
    sketch_size_type sketchSize_;
    hash32_family128 hash_;
};




/*****************************************************************************
 *
 * @brief
 *
 *****************************************************************************/
template<class KmerT = std::uint64_t, class Hash = same_size_hash<KmerT>>
class minimizer_hasher
{
public:
    //---------------------------------------------------------------
    using hasher       = Hash;
    using kmer_type    = KmerT;
    using feature_type = typename std::result_of<hasher(kmer_type)>::type;
    using result_type  = std::array<feature_type,1>;
    //-----------------------------------------------------
    using kmer_size_type   = numk_t;
    using sketch_size_type = typename result_type::size_type;


    //---------------------------------------------------------------
    static constexpr std::uint8_t max_kmer_size() noexcept {
        return max_word_size<kmer_type,2>::value;
    }
    static constexpr sketch_size_type max_sketch_size() noexcept {
        return 1;
    }


    //---------------------------------------------------------------
    explicit
    minimizer_hasher(hasher hash = hasher{}):
        hash_{std::move(hash)},
        k_(max_kmer_size() <= 16 ? 16 : 30)
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
    static sketch_size_type sketch_size() noexcept { return 1; }
    static void sketch_size(sketch_size_type) noexcept { }


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

        result_type sketch {feature_type(~0)};

        for_each_unambiguous_canonical_kmer_2bit<kmer_type>(k_, first, last,
            [&] (kmer_type kmer) {
                if(kmer < sketch[0]) sketch[0] = kmer;
            });

        sketch[0] = hash_(sketch[0]);
        return sketch;
    }


    //---------------------------------------------------------------
    friend void
    write_binary(std::ostream& os, const minimizer_hasher& h)
    {
        write_binary(os, std::uint64_t(h.k_));
    }

    //---------------------------------------------------------------
    friend void
    read_binary(std::istream& is, minimizer_hasher& h)
    {
        std::uint64_t n = 0;
        read_binary(is, n);
        h.k_ = (n <= max_kmer_size()) ? n : max_kmer_size();
    }


private:
    //---------------------------------------------------------------
    hasher hash_;
    kmer_size_type k_;
};






/*****************************************************************************
 *
 * @brief
 *
 *****************************************************************************/
template<class KmerT = std::uint32_t, class Hash = same_size_hash<KmerT>>
class entropy_hasher
{
public:
    //---------------------------------------------------------------
    using kmer_type    = KmerT;
    using hasher       = Hash;
    using feature_type = kmer_type;
    using result_type  = std::vector<feature_type>;
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
    entropy_hasher(hasher hash = hasher{}):
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

        const auto n = distance(first,last);
        if(n < k_) return result_type{};
        const auto s = std::min(sketchSize_, sketch_size_type(n - k_ + 1));

        using pair_t = std::pair<float,kmer_type>;
        std::vector<pair_t> kmers;
        kmers.reserve(s);

        for_each_unambiguous_canonical_kmer_2bit<kmer_type>(k_, first, last,
            [&] (kmer_type kmer) {
                kmers.emplace_back(char_variance(kmer),kmer);
            });

        //sort kmers by entropy
        std::sort(kmers.begin(), kmers.end(),
            [](const pair_t& a, const pair_t& b) { return a.first > b.first; });

        result_type sketch;
        sketch.reserve(s);

        //sketch = s first kmers
        for(sketch_size_type i = 0; i < s; ++i) {
            sketch.push_back(hash_(kmers[i].second));
        }

        return sketch;
    }


    //---------------------------------------------------------------
    friend void
    write_binary(std::ostream& os, const entropy_hasher& h)
    {
        write_binary(os, std::uint64_t(h.k_));
        write_binary(os, std::uint64_t(h.sketchSize_));
    }

    //---------------------------------------------------------------
    friend void
    read_binary(std::istream& is, entropy_hasher& h)
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
    kmer_type
    char_variance(kmer_type kmer) const noexcept
    {
        std::uint8_t c[] {0,0,0,0};

        auto msk = kmer_type(3);
        for(kmer_size_type i = 0; i < k_; ++i) {
            ++c[(kmer & msk)>>(2*i)];
            msk <<= 2;
        }
        auto m = double(c[0] + c[1] + c[2] + c[3]) / 4.0;

        return (c[0]-m)*(c[0]-m) + (c[1]-m)*(c[1]-m) +
               (c[2]-m)*(c[2]-m) + (c[3]-m)*(c[3]-m);
    }


    //---------------------------------------------------------------
    hasher hash_;
    kmer_size_type k_;
    sketch_size_type sketchSize_;
};




/*****************************************************************************
 *
 *
 *****************************************************************************/
template<class T, std::size_t n>
struct kmer_histogram
{
    kmer_histogram(): f_{} { for(T* x = f_; x < f_+n; ++x) *x = 0; }

    T  operator [] (std::size_t i) const noexcept { return f_[i]; }
    T& operator [] (std::size_t i)       noexcept { return f_[i]; }

          T* begin()       noexcept { return f_; }
    const T* begin() const noexcept { return f_; }

          T* end()         noexcept { return f_ + n; }
    const T* end()   const noexcept { return f_ + n; }

private:
    T f_[n];
};

//-------------------------------------------------------------------
template<class T, std::size_t n>
inline void
write_binary(std::ostream& os, const kmer_histogram<T,n>& a)
{
    std::uint64_t l = n;
    os.write(reinterpret_cast<const char*>(&l), sizeof(l));
    for(const auto& x : a) {
        write_binary(os, x);
    }
}

//-------------------------------------------------------------------
template<class T, std::size_t n>
inline void
read_binary(std::istream& is, kmer_histogram<T,n>& a)
{
    std::uint64_t l = 0;
    is.read(reinterpret_cast<char*>(&l), sizeof(l));
    for(auto& x : a) {
        read_binary(is, x);
    }
}

//-------------------------------------------------------------------
template<class T, std::size_t n>
inline std::ostream&
operator << (std::ostream& os, const kmer_histogram<T,n>& a)
{
    os << '[';
    for(const auto& x : a) {
        os << x << ',';
    }
    os << ']';
    return os;
}

}  //namespace mc




/*****************************************************************************
 *
 *
 *****************************************************************************/
namespace std {

template<class T, std::size_t n>
struct hash<mc::kmer_histogram<T,n>> {
    std::size_t
    operator () (const mc::kmer_histogram<T,n>& h) const noexcept {
        std::size_t x = h[0];
        for(std::size_t i = 1; i < n; ++i) {
            x ^= h[i];
        }
        return x;
    }
};

template<class T, std::size_t n>
struct equal_to<mc::kmer_histogram<T,n>> {
    bool
    operator () (const mc::kmer_histogram<T,n>& h1,
                 const mc::kmer_histogram<T,n>& h2) const noexcept
    {
        for(std::size_t i = 0; i < n; ++i) {
            if(h1[i] != h2[i]) return false;
        }
        return true;
    }
};

} // namespace std




namespace mc {

/*****************************************************************************
 *
 * @brief default min-hasher that uses the 'sketch_size' lexicographically
 *
 *        smallest hash values of *one* hash function
 *
 *****************************************************************************/
class single_function_statistics_min_hasher
{
public:
    //---------------------------------------------------------------
    using kmer_type    = std::uint64_t;
    using feature_type = std::uint32_t;
    using hasher       = same_size_hash<feature_type>;
    using result_type  = std::vector<feature_type>;
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
    explicit
    single_function_statistics_min_hasher(hasher hash = hasher{}):
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

        const auto n = distance(first,last);
        if(n < k_) return result_type{};
        const auto s = std::min(sketchSize_, sketch_size_type(n - k_ + 1));

        auto sketch = result_type(s, feature_type(~0));

        for_each_2mer_stat_of_kmer_2bit(k_, first, last,
            [&] (feature_type stat) {
                auto h = hash_(stat);
                if(h < sketch.back()) {
                    auto pos = std::upper_bound(sketch.begin(), sketch.end(), h);
                    if(pos != sketch.end()) {
                        sketch.pop_back();
                        sketch.insert(pos, h);
                    }
                }
            });

        return sketch;
    }


    //---------------------------------------------------------------
    friend void
    write_binary(std::ostream& os, const single_function_statistics_min_hasher& h)
    {
        write_binary(os, std::uint64_t(h.k_));
        write_binary(os, std::uint64_t(h.sketchSize_));
    }

    //---------------------------------------------------------------
    friend void
    read_binary(std::istream& is, single_function_statistics_min_hasher& h)
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

    //---------------------------------------------------------------
    template<class InputIterator>
    static constexpr char
    char2mer2num(InputIterator c) {
        return c2n(c[0]) + 4 * c2n(c[1]);
    }
    template<class InputIterator>
    static constexpr char
    char2mer2num_rev(InputIterator c) {
        return c2n_rev(c[1]) + 4 * c2n_rev(c[0]);
    }

    //---------------------------------------------------------------
    template<class InputIterator>
    static constexpr char
    char3mer2num(InputIterator c) {
        return c2n(c[0]) + 4 * c2n(c[1]) + 16 * c2n(c[2]);
    }
    template<class InputIterator>
    static constexpr char
    char3mer2num_rev(InputIterator c) {
        return c2n_rev(c[2]) + 4 * c2n_rev(c[1]) + 16 * c2n_rev(c[0]);
    }

    //---------------------------------------------------------------
    static constexpr char c2n(char c) noexcept {
        return   (c == 'A') ? 0
               : (c == 'C') ? 1
               : (c == 'G') ? 2
               : (c == 'T') ? 3
               : (c == 'a') ? 0
               : (c == 'c') ? 1
               : (c == 'g') ? 2
               : (c == 't') ? 3
               : 0;
    }
    //---------------------------------------------------------------
    static constexpr char c2n_rev(char c) noexcept {
        return   (c == 'A') ? 3
               : (c == 'C') ? 2
               : (c == 'G') ? 1
               : (c == 'T') ? 0
               : (c == 'a') ? 3
               : (c == 'c') ? 2
               : (c == 'g') ? 1
               : (c == 't') ? 0
               : 3;
    }
};



} // namespace mc


#endif

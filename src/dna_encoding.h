/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2020 André Müller (muellan@uni-mainz.de)
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

#ifndef MC_DNA_ENCODING_H_
#define MC_DNA_ENCODING_H_


#include <vector>
#include <algorithm>

#include "bitmanip.h"


namespace mc {


/** @brief kmer size type */
using numk_t = std::uint8_t;



/*************************************************************************//**
 *
 *
 *****************************************************************************/
template<class CharT>
inline void
reverse_complement(std::basic_string<CharT>& str)
{
    std::reverse(begin(str), end(str));

    for(auto& c : str) {
        switch(c) {
            case 'A': c = 'T'; break;
            case 'a': c = 't'; break;
            case 'C': c = 'G'; break;
            case 'c': c = 'g'; break;
            case 'G': c = 'C'; break;
            case 'g': c = 'c'; break;
            case 'T': c = 'A'; break;
            case 't': c = 'a'; break;
            default: break;
        }
    }
}

//-------------------------------------------------------------------
template<class CharT>
inline std::basic_string<CharT>
make_reverse_complement(std::basic_string<CharT> str)
{
    reverse_complement(str);
    return str;
}

//-------------------------------------------------------------------
// disable all other overloads
template<class T1, class T2>
T1 make_reverse_complement(T1,T2) = delete;



/*************************************************************************//**
 * @return number of kmers in a squence
 *****************************************************************************/
inline constexpr size_t
num_kmers(size_t k, size_t sequenceLen)
{
    return (sequenceLen - k) + 1;
}



/*************************************************************************//**
 *
 *
 * facilities for creating/manipulating 2bit encoding of DNA
 *
 *
 *****************************************************************************/
namespace detail {

template<numk_t k>
struct dna_2bit_encoding {
    static_assert(k > 0,  "k must be at least 1");
    static_assert(k < 33, "k can be at most 32");
};

template<> struct dna_2bit_encoding< 1> { using type = std::uint8_t; };
template<> struct dna_2bit_encoding< 2> { using type = std::uint8_t; };
template<> struct dna_2bit_encoding< 3> { using type = std::uint8_t; };
template<> struct dna_2bit_encoding< 4> { using type = std::uint8_t; };

template<> struct dna_2bit_encoding< 5> { using type = std::uint16_t; };
template<> struct dna_2bit_encoding< 6> { using type = std::uint16_t; };
template<> struct dna_2bit_encoding< 7> { using type = std::uint16_t; };
template<> struct dna_2bit_encoding< 8> { using type = std::uint16_t; };

template<> struct dna_2bit_encoding< 9> { using type = std::uint32_t; };
template<> struct dna_2bit_encoding<10> { using type = std::uint32_t; };
template<> struct dna_2bit_encoding<11> { using type = std::uint32_t; };
template<> struct dna_2bit_encoding<12> { using type = std::uint32_t; };
template<> struct dna_2bit_encoding<13> { using type = std::uint32_t; };
template<> struct dna_2bit_encoding<14> { using type = std::uint32_t; };
template<> struct dna_2bit_encoding<15> { using type = std::uint32_t; };
template<> struct dna_2bit_encoding<16> { using type = std::uint32_t; };

template<> struct dna_2bit_encoding<17> { using type = std::uint64_t; };
template<> struct dna_2bit_encoding<18> { using type = std::uint64_t; };
template<> struct dna_2bit_encoding<19> { using type = std::uint64_t; };
template<> struct dna_2bit_encoding<20> { using type = std::uint64_t; };
template<> struct dna_2bit_encoding<21> { using type = std::uint64_t; };
template<> struct dna_2bit_encoding<22> { using type = std::uint64_t; };
template<> struct dna_2bit_encoding<23> { using type = std::uint64_t; };
template<> struct dna_2bit_encoding<24> { using type = std::uint64_t; };
template<> struct dna_2bit_encoding<25> { using type = std::uint64_t; };
template<> struct dna_2bit_encoding<26> { using type = std::uint64_t; };
template<> struct dna_2bit_encoding<27> { using type = std::uint64_t; };
template<> struct dna_2bit_encoding<28> { using type = std::uint64_t; };
template<> struct dna_2bit_encoding<29> { using type = std::uint64_t; };
template<> struct dna_2bit_encoding<30> { using type = std::uint64_t; };
template<> struct dna_2bit_encoding<31> { using type = std::uint64_t; };
template<> struct dna_2bit_encoding<32> { using type = std::uint64_t; };

} // namespace detail

template<numk_t k>
using dna_2bit_encoding_t = typename detail::dna_2bit_encoding<k>::type;



namespace { //internal linkage

/*************************************************************************//**
 * @param  k : length to consider (in 2-bit letters, so #bits = 2*k)
 * @return reverse complement of kmer
 * original code from Kraken source
 *****************************************************************************/
inline constexpr std::uint64_t
make_reverse_complement_2bit(std::uint64_t s, numk_t k) noexcept
{
    s = ((s >> 2)  & 0x3333333333333333ul) | ((s & 0x3333333333333333ul) << 2);
    s = ((s >> 4)  & 0x0F0F0F0F0F0F0F0Ful) | ((s & 0x0F0F0F0F0F0F0F0Ful) << 4);
    s = ((s >> 8)  & 0x00FF00FF00FF00FFul) | ((s & 0x00FF00FF00FF00FFul) << 8);
    s = ((s >> 16) & 0x0000FFFF0000FFFFul) | ((s & 0x0000FFFF0000FFFFul) << 16);
    s = ((s >> 32) & 0x00000000FFFFFFFFul) | ((s & 0x00000000FFFFFFFFul) << 32);
    return (std::uint64_t(-1) - s) >> (8 * sizeof(s) - (k << 1));
}

//-------------------------------------------------------------------
inline constexpr std::uint32_t
make_reverse_complement_2bit(std::uint32_t s, numk_t k) noexcept
{
    s = ((s >> 2)  & 0x33333333u) | ((s & 0x33333333u) << 2);
    s = ((s >> 4)  & 0x0F0F0F0Fu) | ((s & 0x0F0F0F0Fu) << 4);
    s = ((s >> 8)  & 0x00FF00FFu) | ((s & 0x00FF00FFu) << 8);
    s = ((s >> 16) & 0x0000FFFFu) | ((s & 0x0000FFFFu) << 16);
    return (std::uint32_t(-1) - s) >> (8 * sizeof(s) - (k << 1));
}

//-------------------------------------------------------------------
inline constexpr std::uint16_t
make_reverse_complement_2bit(std::uint16_t s, numk_t k) noexcept
{
    s = ((s >> 2)  & 0x3333u) | ((s & 0x3333u) << 2);
    s = ((s >> 4)  & 0x0F0Fu) | ((s & 0x0F0Fu) << 4);
    s = ((s >> 8)  & 0x00FFu) | ((s & 0x00FFu) << 8);
    return (std::uint16_t(-1) - s) >> (8 * sizeof(s) - (k << 1));
}

//-------------------------------------------------------------------
inline constexpr std::uint8_t
make_reverse_complement_2bit(std::uint8_t s, numk_t k) noexcept
{
    s = ((s >> 2)  & 0x3333u) | ((s & 0x3333u) << 2);
    s = ((s >> 4)  & 0x0F0Fu) | ((s & 0x0F0Fu) << 4);
    return (std::uint8_t(-1) - s) >> (8 * sizeof(s) - (k << 1));
}

} // internal linkage


// disable all other overloads
template<class T1, class T2>
T1 make_reverse_complement_2bit(T1,T2) = delete;



/*************************************************************************//**
 * @param  k : length to consider (in 2-bit letters!, so #bits = 2*k)
 * @return lexicographically smallest of k-mer and reverse complement of k-mer
 *****************************************************************************/
template<class UInt>
inline constexpr UInt
make_canonical_2bit(UInt s, numk_t k) noexcept
{
    static_assert(std::is_integral<UInt>::value &&
                  std::is_unsigned<UInt>::value,
                  "only unsigned integer types are supported");

    auto revcom = make_reverse_complement_2bit(s, k);
    return s < revcom ? s : revcom;
}

//-------------------------------------------------------------------
template<class UInt>
inline constexpr UInt
make_canonical_2bit(UInt s) noexcept
{
    static_assert(std::is_integral<UInt>::value &&
                  std::is_unsigned<UInt>::value,
                  "only unsigned integer types are supported");

    auto revcom = make_reverse_complement_2bit(s, numk_t(sizeof(UInt) * 4));
    return s < revcom ? s : revcom;
}



/*************************************************************************//**
 * @brief canonical kmer comparator
 *****************************************************************************/
struct canonical_2bit_less
{
    using result_type = bool;

    explicit constexpr
    canonical_2bit_less(numk_t k) : k_(k) {}

    constexpr bool
    operator () (std::uint64_t a, std::uint64_t b) const noexcept {
        return (make_canonical_2bit(a,k_) < make_canonical_2bit(b,k_));
    }

    constexpr bool
    operator () (std::uint32_t a, std::uint32_t b) const noexcept {
        return (make_canonical_2bit(a,k_) < make_canonical_2bit(b,k_));
    }

private:
    numk_t k_;
};



/*************************************************************************//**
 * @brief loops through all 2-bit encoded k-mers in a sequence of characters
 *
 * @tparam UInt    result type, must be an unsigned integer type
 *
 * @param k        number of characters in a k-mer
 * @param first    iterator to the first character of the input sequence
 * @param last     iterator to one after the last character of the input sequence
 * @param consume  function object lambda consuming the k-mers
 *****************************************************************************/
template<class UInt, class InputIterator, class Consumer>
inline void
for_each_kmer_2bit(numk_t k,
                   InputIterator first, InputIterator last,
                   Consumer&& consume)
{
    static_assert(std::is_integral<UInt>::value &&
                  std::is_unsigned<UInt>::value,
                  "only unsigned integer types are supported");

    using std::next;

    using ambig_t = half_size_t<UInt>;

    auto kmer    = UInt(0);
    auto kmerMsk = UInt(~0);
    kmerMsk >>= (sizeof(kmerMsk) * CHAR_BIT) - (k * 2);

    auto ambig    = ambig_t(0);  //bitfield marking ambiguous nucleotides
    auto ambigMsk = ambig_t(~0);
    ambigMsk >>= (sizeof(ambigMsk) * CHAR_BIT) - k;

    ++last;
    for(auto ssend = next(first); ssend != last; ++first, ++ssend) {
        //encode next letter
        kmer <<= 2;
        ambig <<= 1;
        switch(*first) {
            case 'A': case 'a': break;
            case 'C': case 'c': kmer |= 1; break;
            case 'G': case 'g': kmer |= 2; break;
            case 'T': case 't': kmer |= 3; break;
            default: ambig |= 1; break;
        }
        --k;
        //make sure we load k letters at the beginning
        if(k == 0) {
            kmer  &= kmerMsk;   //stamp out 2*k lower bits
            ambig &= ambigMsk;  //stamp out k lower bits

            //do something with the kmer (and the ambiguous letters flag)
            consume(kmer, ambig);
            ++k; //we want only one letter next time
        }
    }
}



/*************************************************************************//**
 * @brief loops through all 2-bit encoded k-mers in a sequence of characters
 *
 * @tparam UInt    result type, must be an unsigned integer type
 *
 * @param k        number of characters in a k-mer
 * @param input    input sequence
 * @param consume  function object/lambda consuming (k-mer, ambiguity bitmask)
 *****************************************************************************/
template<class UInt, class InputRange, class Consumer>
inline void
for_each_kmer_2bit(const numk_t k, InputRange input, Consumer&& consume)
{
    using std::begin;
    using std::end;
    for_each_kmer_2bit<UInt>(k, begin(input), end(input),
                                std::forward<Consumer>(consume));
}



/*************************************************************************//**
 * @brief loops through all non-ambiguous 2-bit encoded k-mers in a
 *        sequence of characters
 *
 * @tparam UInt    result type, must be an unsigned integer type
 *
 * @param k        number of characters in a k-mer
 * @param first    iterator to the first character of the input sequence
 * @param last     iterator to one after the last character of the input sequence
 * @param consume  function object (lambda, etc.) consuming the k-mers
 *****************************************************************************/
template<class UInt, class InputIterator, class Consumer>
inline void
for_each_unambiguous_kmer_2bit(
    const numk_t k, InputIterator first, InputIterator last, Consumer&& consume)
{
    for_each_kmer_2bit<UInt>(k, first, last,
        [&] (UInt kmer, half_size_t<UInt> ambig) {
            if(!ambig) consume(kmer);
        });
}

//-------------------------------------------------------------------
template<class UInt, class InputRange, class Consumer>
inline void
for_each_unambiguous_kmer_2bit(
    const numk_t k, InputRange input, Consumer&& consume)
{
    using std::begin;
    using std::end;
    for_each_unambiguous_kmer_2bit<UInt>(k, begin(input), end(input),
                                           std::forward<Consumer>(consume));
}



/*************************************************************************//**
 * @brief loops through all 2-bit encoded k-mers in a sequence of characters
 *        kmers are canonical = min(kmer, reverse_complement(kmer))
 *
 * @tparam UInt    result type, must be an unsigned integer type
 *
 * @param k        number of characters in a k-mer
 * @param first    iterator to the first character of the input sequence
 * @param last     iterator to one after the last character of the input sequence
 * @param consume  function object/lambda consuming (k-mer, ambiguity bitmask)
 *****************************************************************************/
template<class UInt, class InputIterator, class Consumer>
inline void
for_each_canonical_kmer_2bit(const numk_t k,
                             InputIterator first, InputIterator last,
                             Consumer&& consume)
{
    for_each_kmer_2bit<UInt>(k, first, last,
        [&] (UInt kmer, half_size_t<UInt> ambig) {
            consume(make_canonical_2bit(kmer, k), ambig);
        });
}

//-------------------------------------------------------------------
template<class UInt, class InputRange, class Consumer>
inline void
for_each_canonical_kmer_2bit(const numk_t k, InputRange input,
                             Consumer&& consume)
{
    using std::begin;
    using std::end;
    for_each_canonical_kmer_2bit<UInt>(k, begin(input), end(input),
                                       std::forward<Consumer>(consume));
}



/*************************************************************************//**
 * @brief loops through all non-ambiguous 2-bit encoded k-mers in a
 *        sequence of characters
 *        kmers are canonical = min(kmer, reverse_complement(kmer))
 *
 * @tparam UInt    result type, must be an unsigned integer type
 *
 * @param k        number of characters in a k-mer
 * @param input    input sequence
 * @param consume  function object (lambda, etc.) consuming the k-mers
 *****************************************************************************/
template<class UInt, class InputIterator, class Consumer>
inline void
for_each_unambiguous_canonical_kmer_2bit(
    const numk_t k, InputIterator first, InputIterator last, Consumer&& consume)
{
    for_each_kmer_2bit<UInt>(k, first, last,
        [&] (UInt kmer, half_size_t<UInt> ambig) {
            if(!ambig) consume(make_canonical_2bit(kmer, k));
        });
}

//-------------------------------------------------------------------
template<class UInt, class InputRange, class Consumer>
inline void
for_each_unambiguous_canonical_kmer_2bit(
    const numk_t k, InputRange input, Consumer&& consume)
{
    using std::begin;
    using std::end;
    for_each_unambiguous_canonical_kmer_2bit<UInt>(
        k, begin(input), end(input), std::forward<Consumer>(consume));
}

} // namespace mc


#endif

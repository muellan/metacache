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

#ifndef MC_BITMANIP_H_
#define MC_BITMANIP_H_

#include <cstdint>
#include <climits>
#include <string>
#include <type_traits>


namespace mc {


/*************************************************************************//**
 *
 *
 *
 *****************************************************************************/
template<class T, std::uint8_t bitsPerSymbol>
struct max_word_size
{
    static constexpr std::uint8_t value = ((sizeof(T) * CHAR_BIT) / bitsPerSymbol);
};





/*************************************************************************//**
 *
 * @param  x : any integer
 * @return binary representation of 'x' as a string of characters '0' and '1'
 *
 *****************************************************************************/
template<class Int>
inline std::string
bits_as_string(Int x, std::uint8_t groupsize = CHAR_BIT)
{
    static_assert(std::is_integral<Int>::value,
                  "only integer types are supported");

    constexpr auto bits = sizeof(Int) * CHAR_BIT;

    auto s = std::string{};
    s.reserve(bits + (bits / groupsize));

    auto m = Int(1); m <<= (bits-1);

    for(std::uint8_t i = 0; i < bits; ++i, m >>= 1) {
          if(i > 0 && !(i % groupsize)) s.push_back('.');
        s.push_back((x & m) ? '1' : '0');
    }

    return s;
}



/*************************************************************************//**
 *
 * @return integer with 1s on the lowest 'bits' bits and 0s on all other bits
 *
 *****************************************************************************/
template<class Int>
inline Int
lowbitmask(int bits)
{
    static_assert(std::is_integral<Int>::value,
                  "only integer types are supported");

    auto m = Int(1) << (bits - 1);
    return m ^ (m - 1);
}



/*************************************************************************//**
 *
 *
 *
 *****************************************************************************/
namespace detail {

template<class T>
struct half_size;

template<> struct half_size<std::uint8_t>  { using type = std::uint8_t; };
template<> struct half_size<std::uint16_t> { using type = std::uint8_t; };
template<> struct half_size<std::uint32_t> { using type = std::uint16_t; };
template<> struct half_size<std::uint64_t> { using type = std::uint32_t; };

} // namespace detail


template<class UInt>
using half_size_t = typename detail::half_size<UInt>::type;


} // namespace mc


#endif

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

#ifndef MC_IO_HELPERS_H_
#define MC_IO_HELPERS_H_

#include <cstdint>
#include <string>
#include <vector>
#include <array>
#include <iostream>


namespace mc {


/*************************************************************************//**
 *
 *
 *****************************************************************************/
template<class T, class = std::enable_if_t<std::is_fundamental<T>::value,T>>
inline void
write_binary(std::ostream& os, T x)
{
    os.write(reinterpret_cast<const char*>(&x), sizeof(x));
}


//-------------------------------------------------------------------
template<class CharT>
inline void
write_binary(std::ostream& os, const std::basic_string<CharT>& str)
{
    std::uint64_t n = str.size();
    os.write(reinterpret_cast<const char*>(&n), sizeof(n));
    os.write(str.c_str(), str.size() * sizeof(CharT));
}


//-------------------------------------------------------------------
template<class T>
inline void
write_binary(std::ostream& os, const std::vector<T>& v)
{
    std::uint64_t n = v.size();
    os.write(reinterpret_cast<const char*>(&n), sizeof(n));
    if(n > 0)
        os.write(reinterpret_cast<const char*>(v.data()), n * sizeof(T));
}


//-------------------------------------------------------------------
template<class T, std::size_t n>
inline void
write_binary(std::ostream& os, const std::array<T,n>& a)
{
    //dummy
    std::uint64_t l = n;
    os.write(reinterpret_cast<const char*>(&l), sizeof(l));
    if(l > 0)
        os.write(reinterpret_cast<const char*>(a.data()), l * sizeof(T));
}


//-------------------------------------------------------------------
template<class T>
inline void
write_binary(std::ostream& os, const T* a, std::uint64_t n)
{
    if(n > 0)
        os.write(reinterpret_cast<const char*>(a), n * sizeof(T));
}



/*************************************************************************//**
 *
 *
 *****************************************************************************/
template<class T, class = std::enable_if_t<std::is_fundamental<T>::value,T>>
inline void
read_binary(std::istream& is, T& x)
{
    is.read(reinterpret_cast<char*>(&x), sizeof(x));
}


//-------------------------------------------------------------------
template<class CharT>
inline void
read_binary(std::istream& is, std::basic_string<CharT>& str)
{
    std::uint64_t n = 0;
    is.read(reinterpret_cast<char*>(&n), sizeof(n));
    str.clear();
    if(n > 0) {
        str.resize(n);
        is.read(reinterpret_cast<char*>(&(*str.begin())), n * sizeof(CharT));
    }
}


//-------------------------------------------------------------------
template<class T>
inline void
read_binary(std::istream& is, std::vector<T>& v)
{

    std::uint64_t n = 0;
    is.read(reinterpret_cast<char*>(&n), sizeof(n));
    v.resize(n);
    if(n > 0)
        is.read(reinterpret_cast<char*>(v.data()), n * sizeof(T));
}


//-------------------------------------------------------------------
template<class T, std::size_t n>
inline void
read_binary(std::istream& is, std::array<T,n>& a)
{
    //dummy
    std::uint64_t l = 0;
    is.read(reinterpret_cast<char*>(&l), sizeof(l));
    //what if l > n?
    if(l > 0 && l <= n)
        is.read(reinterpret_cast<char*>(a.data()), l * sizeof(T));
}


//-------------------------------------------------------------------
template<class T>
inline void
read_binary(std::istream& is, T* a, std::uint64_t n)
{
    if(n > 0)
        is.read(reinterpret_cast<char*>(a), n * sizeof(T));
}


} // namespace mc

#endif

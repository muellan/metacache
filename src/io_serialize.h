#ifndef MC_IO_HELPERS_H_
#define MC_IO_HELPERS_H_

#include <cstdint>
#include <string>
#include <vector>
#include <array>
#include <iostream>


namespace mc {


/*****************************************************************************
 *
 *
 *****************************************************************************/
template<class T, class = typename
    std::enable_if<std::is_fundamental<T>::value,T>::type>
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
    for(const auto& x : v) {
        write_binary(os, x);
    }
}


//-------------------------------------------------------------------
template<class T, std::size_t n>
inline void
write_binary(std::ostream& os, const std::array<T,n>& a)
{
    //dummy
    std::uint64_t l = n;
    os.write(reinterpret_cast<const char*>(&l), sizeof(l));
    for(const auto& x : a) {
        write_binary(os, x);
    }
}



/*****************************************************************************
 *
 *
 *****************************************************************************/
template<class T, class = typename
    std::enable_if<std::is_fundamental<T>::value,T>::type>
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
    str.resize(n);
    is.read(reinterpret_cast<char*>(&(*str.begin())), n * sizeof(CharT));
}


//-------------------------------------------------------------------
template<class T>
inline void
read_binary(std::istream& is, std::vector<T>& v)
{

    std::uint64_t l = 0;
    is.read(reinterpret_cast<char*>(&l), sizeof(l));
    v.clear();
    v.resize(l);
    for(auto& x : v) {
        read_binary(is, x);
    }
}


//-------------------------------------------------------------------
template<class T, std::size_t n>
inline void
read_binary(std::istream& is, std::array<T,n>& a)
{
    //dummy
    std::uint64_t l = 0;
    is.read(reinterpret_cast<char*>(&l), sizeof(l));
    for(auto& x : a) {
        read_binary(is, x);
    }
}


} // namespace mc

#endif

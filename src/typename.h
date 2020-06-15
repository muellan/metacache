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

#ifndef MC_TYPENAME_H_
#define MC_TYPENAME_H_

#include <type_traits>
#include <typeinfo>
#include <string>


#ifdef __GNUC__
    #include <cxxabi.h>
#endif


namespace mc {


/*************************************************************************//**
 *
 * @brief returnes demangled type name
 *
 *****************************************************************************/
template<class T>
std::string type_name()
{
    if(std::is_same<char,T>::value) return "char";
    if(std::is_same<short int,T>::value) return "short int";
    if(std::is_same<int,T>::value) return "int";
    if(std::is_same<long int,T>::value) return "long int";
    if(std::is_same<long long int,T>::value) return "long long int";

    if(std::is_same<unsigned char,T>::value) return "unsigned char";
    if(std::is_same<unsigned short int,T>::value) return "unsigned short int";
    if(std::is_same<unsigned int,T>::value) return "unsigned int";
    if(std::is_same<unsigned long int,T>::value) return "unsigned long int";
    if(std::is_same<unsigned long long int,T>::value) return "unsigned long long int";

    if(std::is_same<float,T>::value) return "float";
    if(std::is_same<double,T>::value) return "double";
    if(std::is_same<long double,T>::value) return "long double";

#ifdef __GNUC__
    int status;
    char* realname = abi::__cxa_demangle(typeid(T).name(), 0, 0, &status);
    auto str = std::string(realname);
    free(realname);
    return str;
#else
    return std::string( typeid(T).name() );
#endif
}



//-------------------------------------------------------------------
template<class T>
std::string type_name(const T&) {
    return type_name<T>();
}


}  // namespace mc


#endif

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

#ifndef MC_GENERIC_H_
#define MC_GENERIC_H_


#include <array>
#include <cstddef>
#include <utility>
#include <tuple>
#include <memory>


namespace mc {

/*****************************************************************************
 *
 * @brief make_unique from N3656 by Stephan T. Lavavej <stl@microsoft.com>
 *        make_unique will be in C++14
 *
 * @details make_unique<T>(args...)
 *          make_unique<T[]>(n)
 *          make_unique<T[N]>(args...) = delete
 *
 *****************************************************************************/
namespace detail {

template<class T>
struct _Unique_if
{
    using _Single_object = std::unique_ptr<T>;
};

template<class T>
struct _Unique_if<T[]>
{
    using _Unknown_bound = std::unique_ptr<T[]>;
};

template<class T, size_t N>
struct _Unique_if<T[N]>
{
    using _Known_bound = void;
};

}  // namespace detail


//-------------------------------------------------------------------
template<class T, class... Args>
typename detail::_Unique_if<T>::_Single_object
make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

//---------------------------------------------------------
template<class T>
typename detail::_Unique_if<T>::_Unknown_bound
make_unique(size_t n)
{
    using U_t_ = typename std::remove_extent<T>::type;
    return std::unique_ptr<T>(new U_t_[n]());
}

//---------------------------------------------------------
template<class T, class... Args>
typename detail::_Unique_if<T>::_Known_bound
make_unique(Args&&...) = delete;




/*****************************************************************************
 *
 * @brief integer sequence
 *
 *****************************************************************************/
template<class Int, Int... ints>
struct integer_sequence
{
    using value_type = Int;

    static constexpr
    std::size_t size() noexcept {
        return sizeof...(ints);
    }

};


//-------------------------------------------------------------------
template<std::size_t... indices>
using index_sequence = integer_sequence<std::size_t,indices...>;


//-------------------------------------------------------------------
template<class Ostream, class Int>
inline
Ostream&
operator << (Ostream& os, integer_sequence<Int>)
{
    return os;
}

//---------------------------------------------------------
template<class Ostream, class Int, Int i>
inline
Ostream&
operator << (Ostream& os, integer_sequence<Int,i>)
{
    return (os << i);
}

//---------------------------------------------------------
template<class Ostream, class Int, Int i, Int... is>
inline
Ostream&
operator << (Ostream& os, integer_sequence<Int,i,is...>)
{
    return (os << i <<' '<< integer_sequence<Int,is...>());
}






namespace detail {

//-------------------------------------------------------------------
template<class IntegerSequence, typename IntegerSequence::value_type right>
struct integer_appender;

//---------------------------------------------------------
template<class Int, Int... init, Int last>
struct integer_appender<integer_sequence<Int,init...>,last>
{
    using type = integer_sequence<Int,init...,last>;
};

} //namespace detail

//---------------------------------------------------------
template<class IntegerSequence, typename IntegerSequence::value_type right>
using append_integer =
    typename detail::integer_appender<IntegerSequence,right>::type;



namespace detail {


/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
template<class Int, Int min, Int max, bool empty = (min > max)>
struct ascending_integer_sequence_maker
{
    using type = append_integer<
        typename ascending_integer_sequence_maker<Int,min,max-1>::type, max>;
};
//-----------------------------------------------------
template<class Int, Int min, Int max>
struct ascending_integer_sequence_maker<Int,min,max,true>
{
    using type = integer_sequence<Int>;
};
//-----------------------------------------------------
template<class Int, Int min>
struct ascending_integer_sequence_maker<Int,min,min,false>
{
    using type = integer_sequence<Int,min>;
};



/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
template<class Int, Int max, Int min, bool empty = (min > max)>
struct descending_integer_sequence_maker
{
    using type = append_integer<
        typename descending_integer_sequence_maker<Int,max,min+1>::type, min>;
};
//-----------------------------------------------------
template<class Int, Int min, Int max>
struct descending_integer_sequence_maker<Int,min,max,true>
{
    using type = integer_sequence<Int>;
};
//-----------------------------------------------------
template<class Int, Int max>
struct descending_integer_sequence_maker<Int,max,max,false>
{
    using type = integer_sequence<Int,max>;
};



/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
template<class Int, std::size_t size, Int value>
struct uniform_integer_sequence_maker
{
    using type = append_integer<
            typename uniform_integer_sequence_maker<Int,size-1,value>::type,value>;
};
//-----------------------------------------------------
template<class Int, Int value>
struct uniform_integer_sequence_maker<Int,0,value>
{
    using type = integer_sequence<Int>;
};



/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
template<
    class Int,
    std::size_t pos,
    std::size_t size = pos,
    int one = 1,
    int zero = 0
>
struct integer_sequence_mask_maker
{
    using type = append_integer<
            typename integer_sequence_mask_maker<Int,pos,size-1,one,zero>::type,
            (((size-1) == pos) ? one : zero)
        >;
};
//-----------------------------------------------------
template<class Int, std::size_t pos, int one, int zero>
struct integer_sequence_mask_maker<Int,pos,1,one,zero>
{
    using type = integer_sequence<Int,((pos == 0) ? one : zero)>;
};
//-----------------------------------------------------
template<class Int, std::size_t pos, int one, int zero>
struct integer_sequence_mask_maker<Int,pos,0,one,zero>
{
    using type = integer_sequence<Int>;
};

}  // namespace detail




/*****************************************************************************
 *
 *
 *
 *
 *****************************************************************************/
template<class Int, Int maxExcl>
using make_integer_sequence =
    typename detail::ascending_integer_sequence_maker<Int,0,maxExcl-1>::type;

//---------------------------------------------------------
template<std::size_t maxExcl>
using make_index_sequence =
    typename detail::ascending_integer_sequence_maker<
        std::size_t,0,maxExcl-1>::type;


//-------------------------------------------------------------------
template<class Int, Int min, Int max>
using make_ascending_integer_sequence =
    typename detail::ascending_integer_sequence_maker<Int,min,max>::type;

//-----------------------------------------------------
template<class Int, Int max, Int min>
using make_descending_integer_sequence =
    typename detail::descending_integer_sequence_maker<Int,max,min>::type;

//-----------------------------------------------------
template<class Int, std::size_t size, Int value>
using make_uniform_integer_sequence =
    typename detail::uniform_integer_sequence_maker<Int,size,value>::type;

//-----------------------------------------------------
template<
    class Int, std::size_t pos, std::size_t size = pos,
    int one = 1, int zero = 0
>
using make_integer_sequence_mask =
    typename detail::integer_sequence_mask_maker<Int,pos,size,one,zero>::type;




/*****************************************************************************
 *
 *
 *
 *
 *****************************************************************************/
namespace detail {

//-------------------------------------------------------------------
template<class T>
inline constexpr const T&
constant(const T& t, int) {
    return t;
}

//---------------------------------------------------------
/// @brief general case
template<class Target, class Arg, std::size_t... ints>
inline constexpr Target
make_uniform_t_(Target*, const Arg& a, index_sequence<ints...>)
{
    return Target{constant(a,ints)...};
}

//---------------------------------------------------------
/// @brief specialization for std::array
template<class T, std::size_t n, class Arg, std::size_t... ints>
inline constexpr std::array<T,n>
make_uniform_t_(std::array<T,n>*, const Arg& a, index_sequence<ints...>)
{
    return std::array<T,n>{{constant(a,ints)...}};
}


}  // namespace detail

//-----------------------------------------------------
/// @brief calls multi-param ctor of 'Target' with n times 'a'
template<class Target, std::size_t n = std::tuple_size<Target>::value, class Arg>
inline constexpr Target
make_uniform(Arg&& a)
{
    return detail::make_uniform_t_(static_cast<Target*>(nullptr),
        std::forward<Arg>(a), make_index_sequence<n>{});
}




/*****************************************************************************
 *
 * @brief helper function exploding the tuple arguments into a function call
 *
 *****************************************************************************/
namespace detail {
template<class F, class Tuple, std::size_t...ns>
inline auto
apply_helper(F&& f, Tuple&& t, index_sequence<ns...>)
    -> decltype(f(
        std::forward<
            typename std::tuple_element<ns,typename std::decay<Tuple>::type>::type
        >(std::get<ns>(t)) ... ) )
{
    return f(
            std::forward<
            typename std::tuple_element<ns,typename std::decay<Tuple>::type>::type
        >(std::get<ns>(t)) ... );
}
}  // namespace detail



/*****************************************************************************
 *
 * @brief invokes a callable object with a tuple of arguments
 *
 *****************************************************************************/
template<class F>
inline auto
apply(F&& f, std::tuple<>) -> decltype(f())
{
    return f();
}

//---------------------------------------------------------
template<class F, class...Ts, class = typename
         std::enable_if<(sizeof...(Ts) > 0),F>::type>
inline auto
apply(F&& f, std::tuple<Ts...>& t)
    -> decltype(detail::apply_helper(std::forward<F>(f),t,make_index_sequence<sizeof...(Ts)>{}))
{
    return detail::apply_helper(std::forward<F>(f),t,make_index_sequence<sizeof...(Ts)>{});
}

//-----------------------------------------------------
template<class F, class...Ts, class = typename
         std::enable_if<(sizeof...(Ts) > 0),F>::type>
inline auto
apply(F&& f, const std::tuple<Ts...>& t)
    -> decltype(detail::apply_helper(std::forward<F>(f),t,make_index_sequence<sizeof...(Ts)>{}))
{
    return detail::apply_helper(std::forward<F>(f),t,make_index_sequence<sizeof...(Ts)>{});
}


}  // namespace mc


#endif

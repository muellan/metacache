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

#ifndef MC_STATISTICS_MOMENTS_H_
#define MC_STATISTICS_MOMENTS_H_


#include <vector>
#include <initializer_list>
#include <array>
#include <algorithm>
#include <numeric>
#include <string>
#include <utility>
#include <type_traits>
#include <functional>
#include <limits>
#include <iterator>
#include <cstdint>
#include <cmath>


namespace mc {


/*************************************************************************//**
 *
 *
 *
 *****************************************************************************/
using real_t = double;


//-------------------------------------------------------------------
template<class T>
inline real_t
make_real(const T& x) {
    return real_t(x);
}


//-------------------------------------------------------------------
template<class Accumulator>
inline auto
result(const Accumulator& a)
{
    return a.result();
}




/*************************************************************************//**
 *
 * @brief 1st moments
 *
 *****************************************************************************/
template<class InputIterator>
inline auto
mean(InputIterator begin, InputIterator end)
{
    using std::distance;

    using val_t = typename std::iterator_traits<InputIterator>::value_type;

    return make_real(std::accumulate(begin,end,val_t(0)) /
           make_real(distance(begin,end)));
}

//---------------------------------------------------------
template<class InputIterator>
inline auto
raw_moment_1(InputIterator begin, InputIterator end)
{
    using std::distance;

    using val_t = typename std::iterator_traits<InputIterator>::value_type;

    return make_real(std::accumulate(begin,end,val_t(0)) /
           make_real(distance(begin,end)));
}



/*************************************************************************//**
 *
 * @brief 2nd moments
 *
 *****************************************************************************/
template<class InputIterator>
inline auto
raw_moment_2(InputIterator begin, InputIterator end)
{
    using std::distance;

    using fp_t = decltype(make_real(*begin));

    const auto n = fp_t(distance(begin,end));

    if(n < fp_t(1)) return fp_t(0);

    auto s2 = fp_t(0);

    for(; begin != end; ++begin) {
        s2 += (*begin) * (*begin);
    }

    return fp_t(s2 / n);
}

//---------------------------------------------------------
template<class InputIterator>
inline auto
variance(InputIterator begin, InputIterator end)
{
    using std::distance;

    using fp_t = decltype(make_real(*begin));

    const auto n = fp_t(distance(begin,end));

    if(n < fp_t(1)) return fp_t(0);

    auto s = fp_t(0);
    auto s2 = fp_t(0);

    for(; begin != end; ++begin) {
        s += *begin;
        s2 += (*begin) * (*begin);
    }

    return fp_t((s2 - (s*s) /n) /(n - fp_t(1)));
}

//---------------------------------------------------------
template<class InputIterator>
inline auto
central_moment_2(InputIterator begin, InputIterator end)
{
    return variance(begin,end);
}

//---------------------------------------------------------
template<class InputIterator>
inline auto
stddev(InputIterator begin, InputIterator end)
{
    using std::sqrt;

    return sqrt(variance(begin,end));
}



/*************************************************************************//**
 *
 * @brief 3rd moments
 *
 *****************************************************************************/
template<class InputIterator>
inline auto
raw_moment_3(InputIterator begin, InputIterator end)
{
    using std::distance;

    using fp_t = decltype(make_real(*begin));

    const auto n = fp_t(distance(begin,end));

    if(n < fp_t(1)) return fp_t(0);

    auto s3 = fp_t(0);

    for(; begin != end; ++begin) {
        s3 += (*begin) * (*begin) * (*begin);
    }

    return fp_t(s3/n);
}

//---------------------------------------------------------
template<class InputIterator>
inline auto
central_moment_3(InputIterator begin, InputIterator end)
{
    using std::distance;

    using fp_t = decltype(make_real(*begin));

    const auto n = fp_t(distance(begin,end));

    if(n < fp_t(1)) return fp_t(0);

    auto s = fp_t(0);
    auto s2 = fp_t(0);
    auto s3 = fp_t(0);

    for(; begin != end; ++begin) {
        const auto x = *begin;
        s += x;
        x *= x;
        s2 += x;
        s3 += x * (*begin);
    }
    const auto n2 = n*n;
    return fp_t((n2*s3 - fp_t(3)*n*(s2*s2) + fp_t(2)*(s*s))/(n*n2));
}

//---------------------------------------------------------
template<class InputIterator>
inline auto
skewness(InputIterator begin, InputIterator end)
{
    using std::distance;
    using std::pow;

    using fp_t = decltype(make_real(*begin));

    const auto n = fp_t(distance(begin,end));

    if(n < fp_t(1)) return fp_t(0);

    auto s = fp_t(0);
    auto s2 = fp_t(0);
    auto s3 = fp_t(0);

    for(; begin != end; ++begin) {
        const auto x = *begin;
        s += x;
        x *= x;
        s2 += x;
        s3 += x * (*begin);
    }
    const auto n2 = n*n;
    s *= s;

    //2nd central moment
    auto cm2 = (s2 - s /n) /(n - fp_t(1));

    //3rd central moment
    auto cm3 = (n2*s3 - fp_t(3)*n*(s2*s2) + fp_t(2)*s)/(n*n2);

    return fp_t(cm3 / pow(cm2, fp_t(3)/fp_t(2)) );
}



/*************************************************************************//**
 *
 * @brief 4th moments
 *
 *****************************************************************************/
template<class InputIterator>
inline auto
raw_moment_4(InputIterator begin, InputIterator end)
{
    using std::distance;

    using fp_t = decltype(make_real(*begin));

    const auto n = fp_t(distance(begin,end));

    if(n < fp_t(1)) return fp_t(0);

    auto s4 = fp_t(0);

    for(; begin != end; ++begin) {
        s4 += (*begin) * (*begin) * (*begin) * (*begin);
    }

    return fp_t(s4 / n);
}
//---------------------------------------------------------
template<class InputIterator>
inline auto
central_moment_4(InputIterator begin, InputIterator end)
{
    using std::distance;

    using fp_t = decltype(make_real(*begin));

    const auto n = fp_t(distance(begin,end));

    if(n < fp_t(1)) return fp_t(0);

    auto s = fp_t(0);
    auto s2 = fp_t(0);
    auto s3 = fp_t(0);
    auto s4 = fp_t(0);

    for(; begin != end; ++begin) {
        const auto x = *begin;
        s += x;
        x *= x;
        s2 += x;
        s3 += x * (*begin);
        s4 += x * x;
    }
    const auto n2 = n*n;
    const auto ss = s*s ;

    return fp_t( (n2*n*s4 - fp_t(4)*n2*(s*s3) +
                 fp_t(6)*n*(ss*s2) - fp_t(3)*(ss*s)) / (n2*n2) );
}

//---------------------------------------------------------
template<class InputIterator>
inline auto
kurtosis(InputIterator begin, InputIterator end)
{
    using std::distance;

    using fp_t = decltype(make_real(*begin));

    const auto n = fp_t(distance(begin,end));

    if(n < fp_t(1)) return fp_t(0);

    auto s = fp_t(0);
    auto s2 = fp_t(0);
    auto s3 = fp_t(0);
    auto s4 = fp_t(0);

    for(; begin != end; ++begin) {
        const auto x = *begin;
        s += x;
        x *= x;
        s2 += x;
        s3 += x * (*begin);
        s4 += x * x;
    }
    const auto n2 = n*n;
    const auto ss = s*s ;

    //2nd central moment
    const auto cm2 = ((s2 - ss /n) /(n - fp_t(1)));

    //4th central moment
    const auto cm4 = ((n2*n*s4 -
                 fp_t(4)*n2*(s*s3) + fp_t(6)*n*(ss*s2) - fp_t(3)*(ss*s)) /
                (n2*n2));

    return fp_t(cm4 / (cm2*cm2));
}

//---------------------------------------------------------
template<class InputIterator>
inline auto
kurtosis_excess(InputIterator begin, InputIterator end)
{
    return (kurtosis(begin,end) - 3);
}





/****************************************************************************
 *
 *
 * @brief moments accumulator; provides running moments
 *
 *
 ****************************************************************************/
template<
    class Arg = real_t,
    int maxMoment = 2,
    class Size = std::uint_least64_t
>
class moments_accumulator
{};




/*************************************************************************//**
 *
 * @brief specialization for 0 moments (serves as common base class)
 *
 *****************************************************************************/
template<class Arg, class Size>
class moments_accumulator<Arg,0,Size>
{
    using this_t_ = moments_accumulator<Arg,0,Size>;

public:
    //---------------------------------------------------------------
    using argument_type = Arg;
    using result_type = typename std::common_type<real_t,argument_type>::type;
    //-----------------------------------------------------
    using size_type = Size;


    //---------------------------------------------------------------
    constexpr
    moments_accumulator():
        n_(0)
    {}


    //---------------------------------------------------------------
    void
    clear() {
        n_ = 0;
    }


    //---------------------------------------------------------------
    size_type
    size() const noexcept {
        return n_;
    }
    //-----------------------------------------------------
    bool
    empty() const {
        return (n_ < 1);
    }

    //-----------------------------------------------------
    static result_type
    central_moment_0() {
        return 1;
    }


    //---------------------------------------------------------------
    void
    push() {
        ++n_;
    }
    //-----------------------------------------------------
    void
    pop() {
        --n_;
    }


protected:
    //---------------------------------------------------------------
    result_type n() const {
        return n_;
    }


private:
    size_type n_;
};




/*************************************************************************//**
 *
 * @brief first moments only
 *
 *****************************************************************************/
template<class Arg, class Size>
class moments_accumulator<Arg,1,Size> :
    public moments_accumulator<Arg,0,Size>
{
    //---------------------------------------------------------------
    using base_t_ = moments_accumulator<Arg,0,Size>;
    using this_t_ = moments_accumulator<Arg,1,Size>;

public:
    //---------------------------------------------------------------
    using base_t_::n;
    using base_t_::size;


    //---------------------------------------------------------------
    using argument_type = typename base_t_::argument_type;
    using result_type = typename base_t_::result_type;


    //---------------------------------------------------------------
    constexpr
    moments_accumulator():
        base_t_(), sum_(0)
    {}
    //-----------------------------------------------------
    explicit constexpr
    moments_accumulator(const argument_type& t):
        base_t_(), sum_(t)
    {}


    //---------------------------------------------------------------
    // INITIALIZE
    //---------------------------------------------------------------
    void
    clear() {
        base_t_::clear();
        sum_ = 0;
    }
    //-----------------------------------------------------
    this_t_&
    operator = (const argument_type& x) {
        base_t_::clear();
        sum_ = x;
        return *this;
    }


    //---------------------------------------------------------------
    // COLLECT
    //---------------------------------------------------------------
    this_t_&
    operator += (const argument_type& x) {
        push(x);
        return *this;
    }
    //-----------------------------------------------------
    this_t_&
    operator -= (const argument_type& x) {
        pop(x);
        return *this;
    }

    //-----------------------------------------------------
    void
    push(const argument_type& x) {
        base_t_::push();
        sum_ += x;
    }
    //-----------------------------------------------------
    void
    pop(const argument_type& x) {
        base_t_::pop();
        sum_ -= x;
    }


    //---------------------------------------------------------------
    // RESULTS
    //---------------------------------------------------------------
    const result_type&
    sum() const {
        return sum_;
    }
    //-----------------------------------------------------
    static result_type
    central_moment_1() {
        return result_type(0);
    }
    //-----------------------------------------------------
    result_type
    raw_moment_1() const {
        return result_type((size() < 1) ? sum_ : (sum_ / n()));
    }
    //-----------------------------------------------------
    result_type
    mean() const {
        return raw_moment_1();
    }


private:
    result_type sum_;
};




/*************************************************************************//**
 *
 * @brief 1st & 2nd moments
 *
 *****************************************************************************/
template<class Arg, class Size>
class moments_accumulator<Arg,2,Size> :
    public moments_accumulator<Arg,1,Size>
{
    //---------------------------------------------------------------
    using base_t_ = moments_accumulator<Arg,1,Size>;
    using this_t_ = moments_accumulator<Arg,2,Size>;

public:
    //---------------------------------------------------------------
    using base_t_::n;
    using base_t_::size;
    using base_t_::sum;


    //---------------------------------------------------------------
    using argument_type = typename base_t_::argument_type;
    using result_type = typename base_t_::result_type;


    //---------------------------------------------------------------
    constexpr
    moments_accumulator():
        base_t_(), sum2_(0)
    {}
    //-----------------------------------------------------
    explicit constexpr
    moments_accumulator(const argument_type& t):
        base_t_{t}, sum2_(0)
    {}


    //---------------------------------------------------------------
    // INITIALIZE
    //---------------------------------------------------------------
    void
    clear() {
        base_t_::clear();
        sum2_ = 0;
    }
    //-----------------------------------------------------
    this_t_&
    operator = (const argument_type& x) {
        base_t_::operator = (x);
        sum2_ = 0;
        return *this;
    }


    //---------------------------------------------------------------
    // COLLECT
    //---------------------------------------------------------------
    this_t_&
    operator += (const argument_type& x) {
        push(x);
        return *this;
    }
    //-----------------------------------------------------
    this_t_&
    operator -= (const argument_type& x) {
        pop(x);
        return *this;
    }

    //-----------------------------------------------------
    void
    push(const argument_type& x) {
        base_t_::push(x);
        sum2_ += x*x;
    }
    //-----------------------------------------------------
    void
    pop(const argument_type& x) {
        base_t_::pop(x);
        sum2_ -= x*x;
    }


    //---------------------------------------------------------------
    // RESULTS
    //---------------------------------------------------------------
    const result_type&
    sum_2() const {
        return sum2_;
    }
    //-----------------------------------------------------
    result_type
    raw_moment_2() const {
        return result_type((size() < 1) ? sum2_ : (sum2_ / n()) );
    }
    //-----------------------------------------------------
    result_type
    central_moment_2() const {
        auto s2 = sum();
        s2 *= sum();
        return result_type((size() < 1) ? 0 : ((sum2_ - s2 /n()) / (n()-1)) );
    }
    //-----------------------------------------------------
    result_type
    variance() const {
        return central_moment_2();
    }

    //-----------------------------------------------------
    result_type
    stddev() const {
        using std::sqrt;
        return sqrt(central_moment_2());
    }


private:
    result_type sum2_;
};




/****************************************************************************
 *
 * @brief 1st, 2nd & 3rd moments
 *
 *****************************************************************************/
template<class Arg, class Size>
class moments_accumulator<Arg,3,Size> :
    public moments_accumulator<Arg,2,Size>
{
    //---------------------------------------------------------------
    using base_t_ = moments_accumulator<Arg,2,Size>;
    using this_t_ = moments_accumulator<Arg,3,Size>;

public:
    //---------------------------------------------------------------
    using base_t_::n;
    using base_t_::size;
    using base_t_::sum;
    using base_t_::sum_2;
    using base_t_::central_moment_2;


    //---------------------------------------------------------------
    using argument_type = typename base_t_::argument_type;
    using result_type = typename base_t_::result_type;


    //---------------------------------------------------------------
    constexpr
    moments_accumulator():
        base_t_(), sum3_(0)
    {}
    //-----------------------------------------------------
    explicit constexpr
    moments_accumulator(const argument_type& t):
        base_t_{t}, sum3_(0)
    {}


    //---------------------------------------------------------------
    // INITIALIZE
    //---------------------------------------------------------------
    void
    clear() {
        base_t_::clear();
        sum3_ = 0;
    }
    //-----------------------------------------------------
    this_t_&
    operator = (const argument_type& x) {
        base_t_::operator = (x);
        sum3_ = 0;
        return *this;
    }


    //---------------------------------------------------------------
    // COLLECT
    //---------------------------------------------------------------
    this_t_&
    operator += (const argument_type& x) {
        push(x);
        return *this;
    }
    //-----------------------------------------------------
    this_t_&
    operator -= (const argument_type& x) {
        pop(x);
        return *this;
    }

    //-----------------------------------------------------
    void
    pop(const argument_type& x) {
        base_t_::pop(x);
        sum3_ -= x*x*x;
    }
    //-----------------------------------------------------
    void
    push(const argument_type& x) {
        base_t_::push(x);
        sum3_ += x*x*x;
    }


    //---------------------------------------------------------------
    // RESULTS
    //---------------------------------------------------------------
    const result_type&
    sum_3() const {
        return sum3_;
    }
    //-----------------------------------------------------
    result_type
    raw_moment_3() const {
        return result_type((size() < 1) ? sum3_ : (sum3_ / n()) );
    }
    //-----------------------------------------------------
    result_type
    central_moment_3() const {
        if(size() < 2) return result_type(0);

        auto n2 = n() * n();

        return ((    n2 * sum3_
                    - 3 * n() * (sum() * sum_2())
                    + 2 * (sum() * sum() * sum())
                ) / (n() * n2)
            );
    }
    //-----------------------------------------------------
    result_type
    skewness() const {
        using std::pow;

        if(size() < 2) return result_type(0);

        auto n2 = n() * n();
        auto s2 = sum() * sum();

        auto cm2 = (sum_2() - s2 /n()) / (n()-1);

        auto cm3 = (
                (    n2 * sum3_
                    - 3 * n() * (sum() * sum_2())
                    + 2 * (s2 * sum())
                ) / (n() * n2)
            );

        return result_type(cm3 / pow(cm2, 3/2.) );
    }


private:
    result_type sum3_;
};




/****************************************************************************
 *
 * @brief 1st - 4th moments
 *
 *****************************************************************************/
template<class Arg, class Size>
class moments_accumulator<Arg,4,Size> :
    public moments_accumulator<Arg,3,Size>
{
    //---------------------------------------------------------------
    using base_t_ = moments_accumulator<Arg,3,Size>;
    using this_t_ = moments_accumulator<Arg,4,Size>;

public:
    //---------------------------------------------------------------
    using base_t_::n;
    using base_t_::size;
    using base_t_::empty;
    using base_t_::sum;
    using base_t_::sum_2;
    using base_t_::sum_3;
    using base_t_::central_moment_2;


    //---------------------------------------------------------------
    using argument_type = typename base_t_::argument_type;
    using result_type = typename base_t_::result_type;


    //---------------------------------------------------------------
    constexpr
    moments_accumulator():
        base_t_(), sum4_(static_cast<argument_type>(0))
    {}
    //-----------------------------------------------------
    explicit constexpr
    moments_accumulator(const argument_type& t):
        base_t_{t}, sum4_(static_cast<argument_type>(0))
    {}


    //---------------------------------------------------------------
    // INITIALIZE
    //---------------------------------------------------------------
    void
    clear() {
        base_t_::clear();
        sum4_ = 0;
    }
    //-----------------------------------------------------
    this_t_&
    operator = (const argument_type& x) {
        base_t_::operator = (x);
        sum4_ = 0;
        return *this;
    }


    //---------------------------------------------------------------
    // COLLECT
    //---------------------------------------------------------------
    this_t_&
    operator += (const argument_type& x) {
        push(x);
        return *this;
    }
    //-----------------------------------------------------
    this_t_&
    operator -= (const argument_type& x) {
        pop(x);
        return *this;
    }

    //-----------------------------------------------------
    void
    pop(argument_type x) {
        base_t_::pop(x);

        x *= x;
        x *= x;
        sum4_ -= x;
    }
    //-----------------------------------------------------
    void
    push(argument_type x) {
        base_t_::push(x);

        x *= x;
        x *= x;
        sum4_ += x;
    }


    //---------------------------------------------------------------
    // RESULTS
    //---------------------------------------------------------------
    const result_type&
    sum_4() const {
        return sum4_;
    }

    //-----------------------------------------------------
    result_type
    raw_moment_4() const {
        return result_type(n() < 2) ? sum4_ : (sum4_ / n());
    }
    //-----------------------------------------------------
    result_type
    central_moment_4() const {
        if(size() < 2) return result_type(0);

        auto n2 = n()*n();
        auto s2 = sum() * sum();

        return ((    n2 * n() * sum4_
                    - 4 * n2 * sum() * sum_3()
                    + 6 * n() * s2 * sum_2()
                    - 3 * s2 * s2
                ) / (n2*n2)
            );
    }
    //-----------------------------------------------------
    result_type
    kurtosis() const {
        auto n2 = n()*n();
        auto s2 = sum() * sum();

        auto cm2 = (sum_2() - s2 /n()) / (n()-1);
        cm2 *= cm2;

        auto cm4 = ((    n2 * n() * sum4_
                    - 4 * n2 * sum() * sum_3()
                    + 6 * n() * s2 * sum_2()
                    - 3 * s2 * s2
                ) / (n2*n2)
            );

        return result_type(cm4 / cm2);
    }
    //-----------------------------------------------------
    result_type
    kurtosis_excess() const {
        return result_type(kurtosis() - 3);
    }


private:
    result_type sum4_;
};




/*************************************************************************//**
 *
 * convenience definitions
 *
 *****************************************************************************/
template<class Arg = real_t, class Size = std::uint_least64_t>
using mean_accumulator = moments_accumulator<Arg,1,Size>;


template<class Arg = real_t, class Size = std::uint_least64_t>
using variance_accumulator = moments_accumulator<Arg,2,Size>;


template<class Arg = real_t, class Size = std::uint_least64_t>
using skewness_accumulator = moments_accumulator<Arg,3,Size>;


template<class Arg = real_t, class Size = std::uint_least64_t>
using kurtosis_accumulator = moments_accumulator<Arg,4,Size>;


} //namespace mc


#endif

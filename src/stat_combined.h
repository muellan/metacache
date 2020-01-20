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

#ifndef MC_STATISTICS_COMBINED_H_
#define MC_STATISTICS_COMBINED_H_


#include "stat_moments.h"


namespace mc {


/*************************************************************************//**
 *
 * @brief gathers multiple statistical measures from a sample of numbers
 *
 *****************************************************************************/
class statistics_accumulator
{
public:
    //---------------------------------------------------------------
    using argument_type = double;
    using result_type   = argument_type;
    using size_type     = std::uint_least64_t;


    //---------------------------------------------------------------
    statistics_accumulator():
        max_{std::numeric_limits<argument_type>::lowest()}, moments_{}
    {}
    //-----------------------------------------------------
    explicit
    statistics_accumulator(argument_type init):
        max_{init}, moments_{init}
    {}


    //---------------------------------------------------------------
    statistics_accumulator&
    operator = (const statistics_accumulator&) = default;
    //-----------------------------------------------------
    statistics_accumulator&
    operator = (argument_type x) {
        max_ = x;
        moments_ = x;
        return *this;
    }


    //---------------------------------------------------------------
    statistics_accumulator&
    operator += (argument_type x) {
        push(x);
        return *this;
    }
    //-----------------------------------------------------
    void
    push(argument_type x) {
        if(x > max_) max_ = x;
        moments_.push(x);
    }


    //---------------------------------------------------------------
    size_type size() const noexcept {
        return moments_.size();
    }
    //-----------------------------------------------------
    bool empty() const noexcept {
        return moments_.empty();
    }


    //---------------------------------------------------------------
    result_type max()      const { return max_; }

    result_type sum()      const { return moments_.sum(); }

    result_type mean()     const { return moments_.mean(); }
    result_type stddev()   const { return moments_.stddev(); }
    result_type variance() const { return moments_.variance(); }

    result_type skewness() const { return moments_.skewness(); }

private:
    argument_type max_;
    moments_accumulator<argument_type,3,size_type> moments_;

};




} //namespace mc


#endif

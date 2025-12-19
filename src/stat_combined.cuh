/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2022 Robin Kobus  (github.com/funatiq)
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

#ifndef MC_STATISTICS_GPU_HPP_
#define MC_STATISTICS_GPU_HPP_


#include "config.hpp"

#include "../dep/hpc_helpers/include/cuda_helpers.cuh"

#include <cmath>
#include <numeric>


namespace mc {


//-----------------------------------------------------------------------------
/**
 * @brief gathers multiple statistical measures from a sample of numbers
 */
template <policy P>
class statistics_accumulator_gpu
{
    template <policy Q>
    friend class statistics_accumulator_gpu;

public:
    //---------------------------------------------------------------
    using argument_type = uint32_t;
    using result_type   = double;
    using size_type     = uint64_t;


    //---------------------------------------------------------------
    /** @brief allocate memory on host or device depending on policy  */
    statistics_accumulator_gpu ();
    
    statistics_accumulator_gpu (const statistics_accumulator_gpu& other) noexcept :
        data_(other.data_), isCopy_(true)
    {}
    
    statistics_accumulator_gpu (statistics_accumulator_gpu&& other) noexcept :
        data_(other.data_), isCopy_(other.isCopy_)
    {
        other.isCopy_ = true;
    }
    
    /** @brief free memory allocation */
    ~statistics_accumulator_gpu();


    //---------------------------------------------------------------
    statistics_accumulator_gpu&
    operator = (const statistics_accumulator_gpu<policy::Host>& other);

    statistics_accumulator_gpu&
    operator = (const statistics_accumulator_gpu<policy::Device>& other);
    
    statistics_accumulator_gpu&
    operator += (const statistics_accumulator_gpu& other);


    //---------------------------------------------------------------
    HOSTQUALIFIER
    void update (
        size_type n, argument_type max,
        result_type sum, result_type sum2, result_type sum3)
    {
        *size_() += n;
        if (max > *max_()) *max_() = max;
        *sum_()  += sum;
        *sum2_() += sum2;
        *sum3_() += sum3;
    }

    DEVICEQUALIFIER
    void atomic_update(
        size_type n, argument_type max,
        result_type sum, result_type sum2, result_type sum3)
    {
        atomicAdd(size_(), n);
        atomicMax(max_(), max);
        atomicAdd(sum_(), sum);
        atomicAdd(sum2_(), sum2);
        atomicAdd(sum3_(), sum3);
    }


    //---------------------------------------------------------------
    template <class Value>
    void accumulate(Value * values, size_type numValues);


    //---------------------------------------------------------------
    size_type size () const noexcept {
        return *size_();
    }
    
    bool empty () const noexcept {
        return *size_() < 1;
    }
    
    
    //---------------------------------------------------------------
    result_type max () const noexcept {
         return *max_();
    }
   
    result_type sum () const noexcept {
        return *sum_();
    }
   
    result_type mean () const noexcept {
        return result_type((size() < 1) ? sum() : (sum() / n()));
    }
  
    result_type variance () const noexcept {
        auto s2 = sum();
        s2 *= sum();
        return result_type((size() < 1) ? 0 : ((sum_2() - s2 /n()) / (n()-1)) );
    }
   
    result_type stddev () const noexcept {
        using std::sqrt;
        return sqrt(variance());
    }
    
    result_type skewness () const noexcept {
        using std::pow;

        if (size() < 2) return result_type(0);

        auto n2 = n() * n();
        auto s2 = sum() * sum();

        auto cm2 = (sum_2() - s2 /n()) / (n()-1);

        auto cm3 = (
                (    n2 * sum_3()
                    - 3 * n() * (sum() * sum_2())
                    + 2 * (s2 * sum())
                ) / (n() * n2)
            );

        return result_type(cm3 / pow(cm2, 3/2.) );
    }


private:
    //---------------------------------------------------------------
    result_type n () const noexcept {
        return *size_();
    }
    
    result_type sum_2 () const noexcept {
        return *sum2_();
    }
  
    result_type sum_3 () const noexcept {
        return *sum3_();
    }

  
    //---------------------------------------------------------------
    HOSTDEVICEQUALIFIER
    argument_type * max_ () const noexcept {
        return reinterpret_cast<argument_type *>(data_);
    }
 
    HOSTDEVICEQUALIFIER
    size_type * size_ () const noexcept {
        return reinterpret_cast<size_type *>(data_+1);
    }
 
    HOSTDEVICEQUALIFIER
    result_type * sum_ () const noexcept {
        return reinterpret_cast<result_type *>(data_+2);
    }
   
    HOSTDEVICEQUALIFIER
    result_type * sum2_ () const noexcept {
        return reinterpret_cast<result_type *>(data_+3);
    }
   
    HOSTDEVICEQUALIFIER
    result_type * sum3_ () const noexcept {
        return reinterpret_cast<result_type *>(data_+4);
    }

private:
    uint64_t * data_;
    bool isCopy_;
};


} // namespace mc


#endif

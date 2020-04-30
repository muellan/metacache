/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2019 André Müller (muellan@uni-mainz.de)
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


#include "stat_combined.cuh"
#include "../dep/cudahelpers/cuda_helpers.cuh"


namespace mc {


template<>
statistics_accumulator_gpu<policy::Host>::statistics_accumulator_gpu() : isCopy_(false) {
    cudaMallocHost(&data_, 5*sizeof(uint64_t));

    *max_()  = 0;
    *size_() = 0;
    *sum_()  = 0;
    *sum2_() = 0;
    *sum3_() = 0;
}
template<>
statistics_accumulator_gpu<policy::Device>::statistics_accumulator_gpu() : isCopy_(false) {
    cudaMalloc(&data_, 5*sizeof(uint64_t));

    cudaMemset(data_, 0, 5*sizeof(uint64_t));
}


template<>
statistics_accumulator_gpu<policy::Host>::~statistics_accumulator_gpu() {
    if(!isCopy_)
        cudaFreeHost(data_);
}
template<>
statistics_accumulator_gpu<policy::Device>::~statistics_accumulator_gpu() {
    if(!isCopy_)
        cudaFree(data_);
}


template<>
statistics_accumulator_gpu<policy::Host>&
statistics_accumulator_gpu<policy::Host>::operator = (const statistics_accumulator_gpu<policy::Host>& other) {
    cudaMemcpy(data_, other.data_, 5*sizeof(uint64_t), cudaMemcpyDefault);
    return *this;
};
template<>
statistics_accumulator_gpu<policy::Host>&
statistics_accumulator_gpu<policy::Host>::operator = (const statistics_accumulator_gpu<policy::Device>& other) {
    cudaMemcpy(data_, other.data_, 5*sizeof(uint64_t), cudaMemcpyDefault);
    return *this;
};


template<>
statistics_accumulator_gpu<policy::Host>&
statistics_accumulator_gpu<policy::Host>::operator += (const statistics_accumulator_gpu<policy::Host>& other) {
    if(*max_() < *(other.max_())) *max_() = *(other.max_());
    *size_() += *(other.size_());
    *sum_()  += *(other.sum_());
    *sum2_() += *(other.sum2_());
    *sum3_() += *(other.sum3_());

    return *this;
}



template<class Value>
__global__
void accumultate_statistics(statistics_accumulator_gpu<policy::Device> stats, Value * values, size_t numValues)
{
    using value_type  = Value;
    using size_type   = uint64_t;
    using result_type = double;

    __shared__  size_type   s_count;
    __shared__  value_type  s_max;
    __shared__  result_type s_sum;
    __shared__  result_type s_sum2;
    __shared__  result_type s_sum3;

    if(threadIdx.x == 0) {
        s_count = 0;
        s_max = 0;
        s_sum = 0;
        s_sum2 = 0;
        s_sum3 = 0;
    }
    __syncthreads();

    size_type   count = 0;
    value_type  max  = 0;
    result_type sum  = 0;
    result_type sum2 = 0;
    result_type sum3 = 0;

    for(int gid = threadIdx.x + blockIdx.x * blockDim.x;
        gid < numValues;
        gid += blockDim.x * gridDim.x)
    {
        auto x = values[gid];
        ++count;
        max = (x > max) ? x : max;
        sum += x;
        auto x2 = x*x;
        sum2 += x2;
        sum3 += x2*x;
    }

    // warp reduce
    for(int stride = 1; stride < 32; stride *= 2) {
        count += __shfl_down_sync(0xFFFFFFFF, count, stride);
        auto x = __shfl_down_sync(0xFFFFFFFF, max, stride);
        max = (x > max) ? x : max;
        sum  += __shfl_down_sync(0xFFFFFFFF, sum, stride);
        sum2 += __shfl_down_sync(0xFFFFFFFF, sum2, stride);
        sum3 += __shfl_down_sync(0xFFFFFFFF, sum3, stride);
    }

    // block update
    if(threadIdx.x % 32 == 0) {
        atomicAdd(&s_count, count);
        atomicMax(&s_max, max);
        atomicAdd(&s_sum, sum);
        atomicAdd(&s_sum2, sum2);
        atomicAdd(&s_sum3, sum3);
    }

    __syncthreads();

    // global update
    if(threadIdx.x == 0) {
        stats.atomic_update(s_count, s_max, s_sum, s_sum2, s_sum3);
    }
}


template<>
template<class Value>
void statistics_accumulator_gpu<policy::Device>::accumulate(Value * values, size_type numValues) {
    accumultate_statistics<<<512, 512>>>(*this, values, numValues);
}

template void statistics_accumulator_gpu<policy::Device>::accumulate(uint64_t *, size_type);

} //namespace mc

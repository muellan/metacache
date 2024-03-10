/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2024 André Müller (muellan@uni-mainz.de)
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

#ifndef MC_SPAN_H_
#define MC_SPAN_H_


#include <type_traits>
#include <cstddef>
#include <vector>
#include <array>


namespace mc {


template<class Value>
class span {
    using no_const_value_type = typename std::remove_const<Value>::type;
public:
    using value_type = Value;
    using size_type = size_t;

    span() : begin_{nullptr}, size_(0) {}

    explicit span(value_type * first, size_type extend) :
        begin_{first},
        size_{extend}
    {}

    explicit span(std::vector<value_type>& v) :
        begin_{v.data()},
        size_{v.size()}
    {}

    explicit span(const std::vector<no_const_value_type>& v) :
        begin_{v.data()},
        size_{v.size()}
    {}

    template<size_t Extend>
    explicit span(std::array<value_type, Extend>& a) :
        begin_{a.data()},
        size_{Extend}
    {}

    value_type * begin() const noexcept { return begin_; }

    value_type * end() const noexcept { return begin_ + size_; }

    size_t size() const noexcept { return size_; }

    bool empty() const noexcept { return size_ == 0; }

          value_type& operator[] (size_type index)       { return begin_[index]; }
    const value_type& operator[] (size_type index) const { return begin_[index]; }

private:
    value_type * const begin_;
    const size_type size_;
};


} // namespace mc


#endif

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

#ifndef MC_SEQUENCE_VIEW_H_
#define MC_SEQUENCE_VIEW_H_


#include <type_traits>
#include <iterator>
#include <vector>


namespace mc {


/*************************************************************************//**
 *
 *
 *
 *****************************************************************************/
template<class RAIterator, class SizeT = std::size_t>
class sequence_view
{
public:
    using iterator   = RAIterator;
    using size_type  = SizeT;
    using value_type = std::decay_t<decltype(*std::declval<iterator>())>;

    //---------------------------------------------------------------
    explicit
    sequence_view(iterator begin, iterator end) noexcept :
       beg_{begin}, end_{end}
    {}


    //---------------------------------------------------------------
    size_type size() const noexcept {
        using std::distance;
        auto d = distance(beg_,end_);
        return d > 0 ? d : 0;
    }

    //---------------------------------------------------------------
    decltype(auto)
    operator [] (size_type i) const noexcept {
        return beg_[i];
    }


    //---------------------------------------------------------------
    bool empty() const noexcept {
        return beg_ == end_;
    }
    //-----------------------------------------------------
    explicit
    operator bool() const noexcept(noexcept(empty())) {
        return !empty();
    }


    //---------------------------------------------------------------
    iterator begin() const noexcept {
        return beg_;
    }
    iterator end() const noexcept {
        return end_;
    }

    template<class Ostream>
    inline friend Ostream&
    operator << (Ostream& os, const sequence_view& s) {
        for(const auto& x : s) os << x;
        return os;
    }

private:
    iterator beg_;
    iterator end_;
};



/*************************************************************************//**
 *
 *
 *****************************************************************************/
template<class RAIterator>
inline sequence_view<RAIterator>
make_view(RAIterator begin, RAIterator end) noexcept
{
    return sequence_view<RAIterator>{begin, end};
}



/*************************************************************************//**
 *
 *
 *
 *****************************************************************************/
template<class OStream>
inline OStream&
operator << (OStream& os, const std::vector<char>& seq)
{
    for(const auto& x : seq)
        os << x;

    return os;
}


} // namespace mc


#endif

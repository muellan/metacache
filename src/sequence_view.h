#ifndef SRCG_SEQUENCE_VIEW_H_
#define SRCG_SEQUENCE_VIEW_H_


#include <type_traits>


namespace mc {


/*****************************************************************************
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
    using value_type = typename std::decay<decltype(*std::declval<iterator>())>::type;

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
    auto
    operator [] (size_type i) const noexcept
        -> decltype(std::declval<iterator>()[i])
    {
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



/*****************************************************************************
 *
 *
 *****************************************************************************/
template<class RAIterator>
inline sequence_view<RAIterator>
make_view(RAIterator begin, RAIterator end) noexcept
{
    return sequence_view<RAIterator>{begin, end};
}



/*****************************************************************************
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

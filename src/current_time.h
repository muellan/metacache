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

#ifndef MC_CURRENT_DATE_TIME_H_
#define MC_CURRENT_DATE_TIME_H_


#include <iosfwd>
#include <string>
#include <time.h>


namespace mc {


/*****************************************************************************
 *
 * @brief timestamp stream modifier
 *
 *****************************************************************************/
template<class Char, class Traits>
inline std::basic_ostream<Char,Traits>&
time_stamp(std::basic_ostream<Char, Traits>& os)
{
    time_t ts = ::time(0);
    tm* time = localtime(&ts);

    os << (time->tm_year + 1900) << '-'
       << ((time->tm_mon+1 < 10)  ? "0" : "") << time->tm_mon+1  << '-'
       << ((time->tm_mday < 10) ? "0" : "") << time->tm_mday << ' '
       << ((time->tm_hour < 10) ? "0" : "") << time->tm_hour << ':'
       << ((time->tm_min < 10)  ? "0" : "") << time->tm_min  << ':'
       << ((time->tm_sec < 10)  ? "0" : "") << time->tm_sec;

    return os;
}



/*****************************************************************************
 *
 * @brief returns the current date / time
 *
 *****************************************************************************/
class current
{
    //---------------------------------------------------------------
    template<typename T>
    static std::string
    to_2digit_string(const T& value) {
        return ((value < 10) ? std::string("0") : std::string("")) +
               std::to_string(value);
    }

public:
    //---------------------------------------------------------------
    //class is not instantiable
    current() = delete;
    ~current() = delete;
    current(const current&) = delete;
    current(current &&) = delete;
    current& operator () (const current&) = delete;
    current& operator () (current&&) = delete;


    //---------------------------------------------------------------
    using int_t    = std::uint_least16_t;
    using string_t = std::string;


    //---------------------------------------------------------------
    class time_type {
    public:
        time_type(): h_(0), m_(0), s_(0) {}

        explicit
        time_type(int_t h, int_t  m, int_t  s) :
            h_(h), m_(m), s_(s)
        {}

        int_t hours()   const noexcept { return h_; }
        int_t minutes() const noexcept { return m_; }
        int_t seccons() const noexcept { return s_; }

        string_t hh_str() const { return to_2digit_string(h_); }
        string_t mm_str() const { return to_2digit_string(m_); }
        string_t ss_str() const { return to_2digit_string(s_); }

        string_t hh_mm_ss_str() const {
            return (hh_str() + ':' + mm_str() + ':' + ss_str());
        }

        string_t hh_mm_str() const { return (hh_str() + ':' + mm_str()); }

        string_t mm_ss_str() const { return (mm_str() + ':' + ss_str()); }

    private:
        int_t h_, m_, s_;
    };


    //---------------------------------------------------------------
    class date_type {
    public:
        date_type(): y_(0), m_(0), d_(0) {}

        explicit
        date_type(int_t y, int_t m, int_t d) :
            y_(y), m_(m), d_(d)
        {}

        int_t year()  const noexcept { return y_; }
        int_t month() const noexcept { return m_; }
        int_t day()   const noexcept { return d_; }

        string_t yy_str() const { return std::to_string(y_).substr(2); }
        string_t mm_str() const { return to_2digit_string(m_); }
        string_t dd_str() const { return to_2digit_string(d_); }

        string_t yyyy_str() const { return std::to_string(y_); }

        string_t mmmm_str() const {
            switch(m_) {
                case  1: return "January";
                case  2: return "February";
                case  3: return "March";
                case  4: return "April";
                case  5: return "May";
                case  6: return "June";
                case  7: return "July";
                case  8: return "August";
                case  9: return "September";
                case 10: return "October";
                case 11: return "November";
                case 12: return "December";
                default: return "?";
            }
        }

        string_t yy_mm_dd_str() const {
            return (yy_str() + '-' + mm_str()+ '-' + dd_str());
        }
        string_t dd_mm_yy_str() const {
            return (dd_str() + '.' + mm_str() + '.' + yy_str());
        }
        string_t mm_dd_yy_str() const {
            return (mm_str() + '/' + dd_str() + '/' + yy_str());
        }

        string_t yyyy_mm_dd_str() const {
            return (yyyy_str() + '-' + mm_str()+ '-' + dd_str());
        }
        string_t dd_mm_yyyy_str() const {
            return (dd_str() + '.' + mm_str() + '.' + yyyy_str());
        }
        string_t mm_dd_yyyy_str() const {
            return (mm_str() + '/' + dd_str() + '/' + yyyy_str());
        }

        string_t yyyy_mmmm_dd_str() const {
            return (yyyy_str() + '-' + mmmm_str()+ '-' + dd_str());
        }
        string_t dd_mmmm_yyyy_str() const {
            return (dd_str() + '.' + mmmm_str() + '.' + yyyy_str());
        }
        string_t mmmm_dd_yyyy_str() const {
            return (mmmm_str() + '/' + dd_str() + '/' + yyyy_str());
        }

    private:
        int_t y_, m_, d_;
    };


    //---------------------------------------------------------------
    static time_type
    time() {
        time_t ts = ::time(0);
        tm* t = localtime(&ts);
        return time_type{ static_cast<int_t>(t->tm_hour),
                          static_cast<int_t>(t->tm_min),
                          static_cast<int_t>(t->tm_sec) };
    }


    //---------------------------------------------------------------
    static date_type
    date() {
        time_t ts = ::time(0);
        tm* t = localtime(&ts);
        return date_type{ static_cast<int_t>(t->tm_year + 1900),
                          static_cast<int_t>(t->tm_mon + 1),
                          static_cast<int_t>(t->tm_mday) };
    }

};


} //namespace mc

#endif

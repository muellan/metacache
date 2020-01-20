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

#ifndef MC_TIMER_H_
#define MC_TIMER_H_


#include <sstream>
#include <ctime>
#include <ratio>
#include <chrono>


namespace mc {


/*************************************************************************//**
 *
 * @brief simple std::chrono based timer class
 *
 *****************************************************************************/
class timer
{
    using basic_duration_t = std::chrono::microseconds;

public:
    //---------------------------------------------------------------
    void
    start() noexcept {
        if(!running_) {
            running_ = true;
            start_ = std::chrono::high_resolution_clock::now();
        }
    }
    //-----------------------------------------------------
    void
    reset() noexcept {
        total_ = basic_duration_t(0);
        running_ = false;
    }
    //-----------------------------------------------------
    void
    restart() noexcept {
        reset();
        start();
    }

    //-----------------------------------------------------
    void
    stop() noexcept {
        if(running_) {
            stop_ = std::chrono::high_resolution_clock::now();
            total_ += std::chrono::duration_cast<basic_duration_t>(stop_-start_);
            running_ = false;
        }
    }


    //---------------------------------------------------------------
    bool
    running() const noexcept {
        return running_;
    }


    //---------------------------------------------------------------
    template<class Unit>
    Unit
    elapsed() const noexcept {
        return std::chrono::duration_cast<Unit>(current());
    }


    //-----------------------------------------------------
    int64_t
    microseconds() const noexcept {
        return elapsed<std::chrono::microseconds>().count();
    }

    int64_t
    milliseconds() const noexcept {
        return elapsed<std::chrono::milliseconds>().count();
    }

    int64_t
    full_seconds() const noexcept {
        return elapsed<std::chrono::seconds>().count();
    }

    int
    full_minutes() const noexcept {
        return elapsed<std::chrono::minutes>().count();
    }

    int
    full_hours() const noexcept {
        return elapsed<std::chrono::hours>().count();
    }


    //-----------------------------------------------------
    double
    seconds() const noexcept {
        return (milliseconds() / 1000.0);
    }

    double
    minutes() const noexcept {
        return (milliseconds() / 60000.0);
    }

    double
    hours() const noexcept {
        return (milliseconds() / 3600000.0);
    }


    //-----------------------------------------------------
    std::string
    hh_mm_ss() const noexcept
    {
        std::ostringstream ss;
        int h = static_cast<int>(full_hours());
        int m = static_cast<int>(full_minutes());
        int s = static_cast<int>(full_seconds());
        if(h < 10) { ss << "0"; } ss << h << ":";
        if(m < 10) { ss << "0"; } ss << m << ":";
        if(s < 10) { ss << "0"; } ss << s;
        return ss.str();
    }


private:
    //-----------------------------------------------------
    basic_duration_t
    current() const noexcept {
        if(!running_) return total_;

        return (
            total_ + (std::chrono::duration_cast<basic_duration_t>(
                std::chrono::high_resolution_clock::now() - start_)) );
    }

    //---------------------------------------------------------------
    std::chrono::high_resolution_clock::time_point start_, stop_;
    basic_duration_t total_ = basic_duration_t(0);
    bool running_ = false;
};


}  // namespace mc


#endif

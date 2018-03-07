/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) Copyright (C) 2016-2018 André Müller (muellan@uni-mainz.de)
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

#ifndef MC_ARGS_PARSER_H_
#define MC_ARGS_PARSER_H_

#include <string>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdint>
#include <iostream>


namespace mc {

namespace detail {

/*************************************************************************//**
 *
 * @brief type conversion helper
 *
 *****************************************************************************/
struct convert_string
{
    template<class T>
    static inline T
    to(const std::string& s) {
        return T(s.c_str());
    }
};


//---------------------------------------------------------
template<> inline
bool
convert_string::to<bool>(const std::string& s) {
    return static_cast<bool>(std::atoi(s.c_str()));
}


//---------------------------------------------------------
template<> inline
unsigned char
convert_string::to<unsigned char>(const std::string& s) {
    return static_cast<unsigned char>(std::atoi(s.c_str()));
}

//---------------------------------------------------------
template<> inline
unsigned short int
convert_string::to<unsigned short int>(const std::string& s) {
    return static_cast<unsigned short int>(std::atoi(s.c_str()));
}
//---------------------------------------------------------
template<> inline
unsigned int
convert_string::to<unsigned int>(const std::string& s) {
    return std::atoi(s.c_str());
}
//---------------------------------------------------------
template<> inline
unsigned long int
convert_string::to<unsigned long int>(const std::string& s) {
    return std::atol(s.c_str());
}
//---------------------------------------------------------
template<> inline
unsigned long long int
convert_string::to<unsigned long long int>(const std::string& s) {
    return std::atol(s.c_str());
}


//---------------------------------------------------------
template<> inline
char
convert_string::to<char>(const std::string& s) {
    return static_cast<char>(std::atoi(s.c_str()));
}
//---------------------------------------------------------
template<> inline
short int
convert_string::to<short int>(const std::string& s) {
    return static_cast<short int>(std::atoi(s.c_str()));
}
//---------------------------------------------------------
template<> inline
int
convert_string::to<int>(const std::string& s) {
    return std::atoi(s.c_str());
}
//---------------------------------------------------------
template<> inline
long int
convert_string::to<long int>(const std::string& s) {
    return std::atol(s.c_str());
}
//---------------------------------------------------------
template<> inline
long long int
convert_string::to<long long int>(const std::string& s) {
    return std::atol(s.c_str());
}


//---------------------------------------------------------
template<> inline
float
convert_string::to<float>(const std::string& s) {
    return std::atof(s.c_str());
}
//---------------------------------------------------------
template<> inline
double
convert_string::to<double>(const std::string& s) {
    return std::atof(s.c_str());
}
//---------------------------------------------------------
template<> inline
long double
convert_string::to<long double>(const std::string& s) {
    return static_cast<long double>(std::atof(s.c_str()));
}

//---------------------------------------------------------
template<> inline
std::string
convert_string::to<std::string>(const std::string& s) {
    return s;
}



/*************************************************************************//**
 *
 * @brief default values for fundamental types
 *
 *****************************************************************************/
template<class T, bool = std::is_fundamental<T>::value>
struct make_default { static constexpr T value = T(); };

template<class T>
struct make_default<T,true> { static constexpr T value = T(0); };

}  // namespace detai





/*************************************************************************//**
 *
 * @brief simple ad-hoc command line argument parser
 *
 * 2008-2017 André Müller
 *
 * The parser distinguishes between non-prefixed and prefixed arguments
 * based on the prefix setting. The default prefix is '-'.
 *
 *
 * CONSTRUCTION
 *     args_parser args(int argc, char** argv)
 *
 *
 * READING PARAMETERS
 *     args.contains("parameter_string")
 *
 *     args.get<type>(non_prefixed_parameter_index)
 *     args.get<type>(non_prefixed_parameter_index, not_provided_value)
 *
 *     args.get<type>("parameter_string", not_provided_value)
 *     args.get<type>("parameter_string", not_provided_value, default_value)
 *
 *     not_provided_value: if the parameter string couldn't be found
 *
 *     default_value:      if the parameter string was found but
 *                         no parameter value was defined
 *
 *   use the member function
 *     arg_prefix(char)
 *   to set the character which signals the begin of a prefixed parameter
 *
 *
 * EXAMPLES
 *
 *   call string:
 *      "myExecutable -a -b 12 -c 20.3 -d filename"
 *
 *   query parameter count:
 *
 *      args.size()                  returns 7
 *      args.prefixed_count()        returns 4   (a, b, c, d)
 *      args.non_prefixed_count()    returns 3   (12, 20.3, filename)
 *
 *
 *   access prefixed arguments (and following value agument to the right):
 *
 *      args.contains("z")        returns false (no such argument provided)
 *      args.get<int>("z", 1)     returns 1
 *      args.get<int>("z", 1, 2)  returns 1
 *
 *      args.contains("a")                 returns true
 *      args.get<string>("a", "xyz")       returns the string "xyz"
 *      args.get<int>("a", 10)             returns the integer 10
 *      args.get<double>("a", 10.0)        returns the double 10.0
 *      args.get<double>("a", 10.0, 88.0)  returns the double 88.0
 *                                         would 10.0 if "a" wasn't found
 *
 *      args.contains("b")                 returns true
 *      args.get<double>("b", 99.0)        returns the double 12.0
 *      args.get<int>("b", 99)             returns the integer 12
 *
 *      args.contains("c")                 returns true
 *      argc.get<double>("c", 0.0)         returns the double 20.3
 *      argc.get<int>("c", 0)              returns the integer 20
 *
 *      args.contains("d")                 returns true
 *      args.get<double>("d", "default")   returns the string "filename"
 *
 *
 *   access non-prefixed arguments (arguments without the prefix token "-"):
 *
 *      args.non_prefixed_count()          returns 2
 *      args.non_prefixed(0)               returns "20.3"
 *      args.non_prefixed(1)               returns "filename"
 *
 *      args.get<double>(0)                returns 20.3
 *      args.get<std::string>(1)           returns "filename"
 *
 *
 *****************************************************************************/
class args_parser
{
    using args_store = std::vector<std::string>;

public:
    //---------------------------------------------------------------
    using size_type = args_store::size_type;
    using string    = args_store::value_type;


    //---------------------------------------------------------------
    explicit
    args_parser(int argc, char** argv):
        args_(),
        prefix_('-'), listDelimiter_(',')
    {
        if(argc > 0) {
            //ignore first arg
            args_.reserve(argc-1);
            for(int i = 1; i < argc; ++i) {
                args_.emplace_back(argv[i]);
            }
        }
    }

    //---------------------------------------------------------------
    explicit
    args_parser(std::vector<std::string> args):
        args_(std::move(args)),
        prefix_('-'), listDelimiter_(',')
    {}


    //---------------------------------------------------------------
    void arg_prefix(char c) {
        prefix_ = c;
    }
    char arg_prefix() const noexcept {
        return prefix_;
    }


    //---------------------------------------------------------------
    void list_delimiter(char c) {
        listDelimiter_ = c;
    }
    char list_delimiter() const noexcept {
        return listDelimiter_;
    }


    //---------------------------------------------------------------
    size_type size() const noexcept {
        return args_.size();
    }
    //-------------------------------------------
    size_type non_prefixed_count() const noexcept {
        size_type n = 0;
        for(const auto& arg : args_)
            if(!arg.empty() && arg[0] != prefix_) ++n;
        return n;
    }
    //-------------------------------------------
    size_type prefixed_count() const noexcept {
        size_type n = 0;
        for(const auto& arg : args_)
            if(!arg.empty() && arg[0] == prefix_) ++n;
        return n;
    }


    //---------------------------------------------------------------
    bool contains(const string& arg) const {
        return parse(arg);
    }

    bool contains(std::initializer_list<string> args) const {
        for(const auto& arg : args) {
            if(contains(arg)) return true;
        }
        return false;
    }


    //---------------------------------------------------------------
    bool contains_singleton(const string& arg) const {
        return parse(arg,false);
    }

    bool contains_singleton(std::initializer_list<string> args) const {
        for(const auto& arg : args) {
            if(contains_singleton(arg)) return true;
        }
        return false;
    }


    //---------------------------------------------------------------
    bool is_prefixed(size_type i) const noexcept {
        return (args_[i][0] == prefix_);
    }
    //-----------------------------------------------------
    bool is_preceded_by_prefixed_arg(size_type i) const noexcept {
        return (args_[i][0] == prefix_);
    }


    //---------------------------------------------------------------
    // get parameter
    //---------------------------------------------------------------
    string
    operator [] (size_type i) const {
        return (i < args_.size()) ? string(args_[i]) : string("");
    }
    //-----------------------------------------------------
    string
    non_prefixed(size_type i, const string& defaultValue = "") const {
        return non_prefixed_str(i, defaultValue);
    }
    //-----------------------------------------------------
    string
    prefixed(size_type i, const string& defaultValue = "") const {
        return prefixed_str(i, defaultValue);
    }


    //---------------------------------------------------------------
    template<class T>
    T get(size_type i,
          const T& notProvidedValue = detail::make_default<T>::value) const
    {
        string p = non_prefixed_str(i);
        if(p != "")
            return detail::convert_string::to<T>(p);
        else
            return notProvidedValue;
    }


    //---------------------------------------------------------------
    /**
     * @param arg               argument name to find
     * @param notProvidedValue  return value if either arg name was not found
     *                          nor any associated value
     */
    template<class T>
    T get(const string& arg,
          const T& notProvidedValue = detail::make_default<T>::value) const
    {
        string p;
        if(parse(arg, p)) {
            if(!p.empty())
                return detail::convert_string::to<T>(p);
            else
                return notProvidedValue;
        }
        return notProvidedValue;
    }


    //---------------------------------------------------------------
     /**
     * @param arg               argument name to find
     * @param notProvidedValue  return value if arg name was not found
     * @param defaultValue      return value if arg name was found,
     *                          but no associated value
     */
    template<class T>
    T get(const string& arg,
          const T& notProvidedValue,
          const T& defaultValue) const
    {
        string p;
        if(parse(arg, p)) {
            if(!p.empty())
                return detail::convert_string::to<T>(p);
            else
                return defaultValue;
        }
        return notProvidedValue;
    }


    //---------------------------------------------------------------
    /**
     * @param args              list of alternative argument names
     * @param notProvidedValue  return value if either no arg name was found
     *                          or nor any associated value
     */
    template<class T>
    T get(std::initializer_list<string> args,
          const T& notProvidedValue = detail::make_default<T>::value) const
    {
        T res = notProvidedValue;
        for(const auto& arg : args) {
            res = get<T>(arg, res);
        }
        return res;
    }


    //---------------------------------------------------------------
    /**
     * @param args              list of alternative argument names
     * @param notProvidedValue  return value if no arg name was found
     * @param defaultValue      return value if arg name was found,
     *                          but no associated value
     */
    template<class T>
    T get(std::initializer_list<string> args,
          const T& notProvidedValue,
          const T& defaultValue) const
    {
        T res = notProvidedValue;
        for(const auto& arg : args) {
            res = get<T>(arg, res, defaultValue);
        }
        return res;
    }


    //---------------------------------------------------------------
    /**
     * @param  arg  argument name to find
     * @return
     */
    template<class T>
    std::vector<T> get_list(const string& arg) const
    {
        string p;
        if(parse(arg, p)) {
            std::vector<T> v;
            size_type pos = 0;
            for(size_type i = 0; i < p.length(); ++i) {
                size_type j = p.find(listDelimiter_,pos);
                if(j == string::npos) {
                    j = p.length();
                    v.push_back(detail::convert_string::to<T>(p.substr(pos,j-pos)));
                    break;
                }
                v.push_back(detail::convert_string::to<T>(p.substr(pos,j-pos)));
                pos = j+1;
            }
            return v;
        }
        return std::vector<T>();
    }


private:
    //---------------------------------------------------------------
    bool
    parse(const string& arg, bool subarg = true) const {
        if(args_.empty()) return false;

        for(const auto& a : args_) {
            if(a[0] == prefix_) {
                string s = "";
                for(size_type j = 1; j < a.size(); ++j) {
                    s += a[j];
                    if(subarg && s == arg) return true;
                }
                if(s == arg) return true;
            }
        }
        return false;
    }


    //---------------------------------------------------------------
    bool
    parse(const string& arg, string& value) const
    {
        if(args_.empty()) return false;

        for(size_type i = 0; i < args_.size(); ++i) {
            if(args_[i][0] == prefix_) {
                string s = "";
                for(size_type j = 1; j < args_[i].size(); ++j) {
                    s += args_[i][j];
                    if(s == arg) {
                        string p(args_[i]);
                        value = p.substr(j+1, args_[i].size()-j-1);
                        for(size_type k = i+1; k < args_.size(); ++k) {
                            if(args_[k][0] == prefix_) break;
                            string ps(args_[k]);
                            if(value.length() > 0) value += " ";
                            value += ps;
                        }
                        return true;
                    }
                }
            }
        }
        return false;
    }

    //---------------------------------------------------------------
    string
    prefixed_str(size_type index, const string& def = "") const
    {
        if(args_.empty()) return def;

        size_type i = 0;
        for(const auto& arg : args_) {
            if(arg[0] == prefix_) {
                if(i == index) return arg;
                ++i;
            }
        }
        return def;
    }
    //---------------------------------------------------------------
    string
    non_prefixed_str(size_type index, const string& def = "") const
    {
        if(args_.empty()) return def;

        size_type j = 0;
        for(const auto& arg : args_) {
            if(arg[0] != prefix_) {
                if(j == index) return arg;
                ++j;
            }
        }
        return def;
    }

    //---------------------------------------------------------------
    args_store args_;
    char prefix_;
    char listDelimiter_;
};



} //namespace mc

#endif


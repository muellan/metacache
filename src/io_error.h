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

#ifndef MC_IO_ERROR_H_
#define MC_IO_ERROR_H_


#include <string>
#include <stdexcept>


namespace mc {


/*************************************************************************//**
 *
 *
 *
 *****************************************************************************/
class io_error :
    public std::runtime_error
{
public:
    explicit
    io_error(const std::string& what):
        std::runtime_error(what)
    {}
};



/*************************************************************************//**
 *
 *
 *
 *****************************************************************************/
class io_format_error :
    public io_error
{
public:
    explicit
    io_format_error(const std::string& what):
        io_error(what)
    {}

};



/*************************************************************************//**
 *
 *
 *
 *****************************************************************************/
class file_io_error :
    public io_error
{
public:
    explicit
    file_io_error(const std::string& what):
        io_error(what)
    {}

    explicit
    file_io_error(const std::string& what, const std::string& filename):
        io_error(what), filename_(filename)
    {}

    const char* filename() const noexcept {
        return filename_.c_str();
    }

private:
    std::string filename_;
};



/*************************************************************************//**
 *
 *
 *
 *****************************************************************************/
class file_access_error :
    public file_io_error
{
public:
    explicit
    file_access_error(const std::string& what):
        file_io_error(what)
    {}

    explicit
    file_access_error(const std::string& what, const std::string& filename):
        file_io_error(what,filename)
    {}
};



/*************************************************************************//**
 *
 *
 *
 *****************************************************************************/
class file_read_error :
    public file_io_error
{
public:
    explicit
    file_read_error(const std::string& what):
        file_io_error(what)
    {}

    explicit
    file_read_error(const std::string& what, const std::string& filename):
        file_io_error(what,filename)
    {}
};



/*************************************************************************//**
 *
 *
 *
 *****************************************************************************/
class file_write_error :
    public file_io_error
{
public:
    explicit
    file_write_error(const std::string& what):
        file_io_error(what)
    {}

    explicit
    file_write_error(const std::string& what, const std::string& filename):
        file_io_error(what,filename)
    {}
};


}  // namespace mc



#endif

/*******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2019 André Müller (muellan@uni-mainz.de)
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
 * FASTA / FASTQ file parsing inspired by https://github.com/attractivechaos/klib
 * by Attractive Chaos <attractor@live.co.uk>
 *
 *****************************************************************************/

#ifndef CHAR_ISTREAM_H
#define CHAR_ISTREAM_H

#ifdef MC_ZLIB
#include <zlib.h>
#else
#include <stdio.h>
#endif

#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include <memory>
#include <string>



namespace mc {



class char_sequence
{
public:
    using size_type = std::size_t;
    using char_type = char;

    using value_type     = char_type;
    using iterator       = value_type*;
    using const_iterator = const value_type*;

private:
    char_type* data_;
    size_type size_;
    size_type capacity_;

public:
    char_sequence() :
        data_{nullptr}, size_{0}, capacity_{0}
    {}

    char_sequence(const char_sequence& other) : char_sequence()
    {
        reserve(other.capacity_);
        size_ = other.size_;
        memcpy(data_, other.data_, size_*sizeof(char_type));
    }

    char_sequence(char_sequence&& other):
        data_{other.data_},
        size_{other.size_}, capacity_{other.capacity_}
    {
        other.data_ = nullptr;
        other.capacity_ = 0;
        other.size_ = 0;
    }

    ~char_sequence() {
        free(data_);
    }

    char_sequence& operator = (const char_sequence& other)
    {
        char_sequence temp{other};
        std::swap(temp.data_, data_);
        size_ = temp.size_;
        capacity_ = temp.capacity_;
        return *this;
    }

    char_sequence& operator = (char_sequence&& other)
    {
        std::swap(data_, other.data_);
        std::swap(capacity_, other.capacity_);
        other.size_ = 0;
        return *this;
    }

    size_type capacity() const noexcept { return capacity_; }

    size_type size() const noexcept { return size_; }

    bool empty() const noexcept { return size_ < 1; }

    void resize(size_type newSize) {
        reserve(newSize);
        size_ = newSize;
    }

    void clear() noexcept { resize(0); }

    void reserve(size_type newCapacity) {
        if(capacity_ < newCapacity) {
            capacity_ = round_up_32(newCapacity);
            data_ = (char_type*)realloc(data_, capacity_*sizeof(char_type));
        }
    }

    void reserve_exactly(size_type newCapacity) {
        if(capacity_ < newCapacity) {
            capacity_ = newCapacity;
            data_ = (char_type*)realloc(data_, capacity_*sizeof(char_type));
        }
    }

    void push_back(char_type c) {
        const auto oldSize = size_;
        resize(oldSize+1);
        data_[oldSize] = c;
    }

    void append(char_type* first, char_type* last) {
        append(first, last - first);
    }

    void append(char_type* first, size_type n) {
        const auto oldSize = size_;
        resize(oldSize + n);
        memcpy(data_ + oldSize, first, n);
    }

          char_type* data()       noexcept { return data_; }
    const char_type* data() const noexcept { return data_; }

    char_type& operator[] (size_type pos)       noexcept { return data_[pos]; }
    char_type  operator[] (size_type pos) const noexcept { return data_[pos]; }

          iterator  begin()       noexcept { return data_; }
    const_iterator  begin() const noexcept { return data_; }
    const_iterator cbegin() const noexcept { return data_; }

          iterator  end()       noexcept { return data_ + size_; }
    const_iterator  end() const noexcept { return data_ + size_; }
    const_iterator cend() const noexcept { return data_ + size_; }

private:
    constexpr size_type round_up_32(size_type x) {
        --x;
        x|=x>>1;
        x|=x>>2;
        x|=x>>4;
        x|=x>>8;
        x|=x>>16;
        ++x;
        return x;
    }
};


class char_istream
{
public:
    enum class Status : char {
        no_err   = 0, // no error
        eof      = 1, // end-of-file
        err      = 2  // error reading stream
    };

    struct filehandle
    {
        filehandle() : file_{NULL} {}
        filehandle(filehandle&& other) :
            file_{other.file_}
        {
            other.file_ = NULL;
        }
        ~filehandle() {
            if(file_) close();
        }

#ifdef MC_ZLIB
        gzFile file_;

        Status open(const char *filename) {
            if(file_) close();
            file_ = gzopen(filename, "r");
            if (file_ == Z_NULL)
                return Status::err;
            return Status::no_err;
        }
        void close() {
            gzclose(file_);
        }
        int read(void *buffer, size_t size) {
            return gzread(file_, buffer, size);
        }
#else
        FILE* file_;

        Status open(const char *filename) {
            if(file_) close();
            file_ = fopen(filename, "r");
            if (file_ == NULL)
                return Status::err;
            return Status::no_err;
        }
        void close() {
            fclose(file_);
        }
        int read(void *buffer, size_t size) {
            size_t readSize = fread(buffer, 1, size, file_);
            if(ferror(file_))
                return -1;
            else
                return readSize;
        }
#endif
    };

    filehandle filehandle_;
    std::unique_ptr<unsigned char[]> buf_;
    int begin_, end_;
    char lastChar_;
    Status status_;

    static constexpr size_t bufsize() noexcept {
        return 16384;
    }

    char_istream() :
        filehandle_{}, buf_{nullptr},
        begin_{0}, end_{0},
        lastChar_{0},
        status_{Status::err}
    {}

    explicit char_istream(const char* filename) : char_istream()
    {
        open(filename);
    }

    char_istream(char_istream&& other) :
        filehandle_{std::move(other.filehandle_)},
        buf_{std::move(other.buf_)},
        begin_{other.begin_},
        end_{other.end_},
        lastChar_{other.lastChar_},
        status_{other.status_}
    {
        other.filehandle_.file_ = NULL;
        other.buf_ = nullptr;
        other.begin_ = 0;
        other.end_ = 0;
        other.lastChar_ = 0;
        other.status_ = Status::err;
    }

    char_istream& operator = (char_istream&& other)
    {
        std::swap(other.filehandle_.file_, filehandle_.file_);
        std::swap(other.buf_, buf_);
        std::swap(other.begin_, begin_);
        std::swap(other.end_, end_);
        std::swap(other.lastChar_, lastChar_);
        std::swap(other.status_, status_);
        return *this;
    }

    void open(const char* filename) {
        status_ = filehandle_.open(filename);
        buf_ = std::make_unique<unsigned char[]>(bufsize());
        begin_ = 0;
        end_ = 0;
        lastChar_ = 0;
    }

    Status status() const noexcept { return status_; }
    bool good() const noexcept { return status_ == Status::no_err; }
    char last_char() const noexcept { return lastChar_; }

private:
    bool buffer_empty() const noexcept { return begin_ >= end_; }

    bool validate_buffer() {
        if(status_ != Status::no_err) return false;

        if (buffer_empty()) {
            begin_ = 0;

            end_ = filehandle_.read(buf_.get(), bufsize());

            if (end_ == 0) { // nothing read
                status_ = Status::eof;
                return false;
            }
            if (end_ == -1) { // read error
                status_ = Status::err;
                return false;
            }
        }

        return true;
    }

public:
    void read_char()
    {
        if(validate_buffer())
            lastChar_ = buf_[(begin_++)];
        else
            lastChar_ = 0;
    }

    void append_line(char_sequence& str) {
        append_line_impl(str);
    }

    void append_line(std::string& str) {
        append_line_impl(str);
    }

private:
    template<class StringT>
    void append_line_impl(StringT& str)
    {
        bool gotany = false;
        for (;;) {
            if(!validate_buffer()) break;

            // find line seperator
            unsigned char *sep = (unsigned char*)memchr(buf_.get() + begin_, '\n', end_ - begin_);
            // if seperator was not found, set lineEnd to buffer end
            int lineEnd = (sep != 0) ? (sep - buf_.get()) : end_;

            // append to str
            str.append((char*)buf_.get() + begin_, lineEnd - begin_);

            gotany = true;

            // advance buffer begin
            begin_ = lineEnd + 1;

            if (lineEnd < end_) {
                break;
            }
            // else: line continues
        }
        if (gotany && str.size() > 0 && str[str.size()-1] == '\r') {
            str.resize(str.size()-1);
        }
    }

public:
    void skip_line() {
        for (;;) {
            if(!validate_buffer()) break;

            // find line seperator
            unsigned char *sep = (unsigned char*)memchr(buf_.get() + begin_, '\n', end_ - begin_);
            // if seperator was not found, set lineEnd to buffer end
            int lineEnd = (sep != 0) ? (sep - buf_.get()) : end_;
            // advance buffer begin
            begin_ = lineEnd + 1;

            if (lineEnd < end_) {
                break;
            }
            // else: line continues
        }
    }
};


} //namespace mc

#endif

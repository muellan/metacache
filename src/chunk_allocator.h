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

#ifndef MC_CHUNK_ALLOCATOR_H_
#define MC_CHUNK_ALLOCATOR_H_

#include <vector>
#include <memory>
#include <algorithm>
#include <iostream>
//#include <mutex>


namespace mc {

/*************************************************************************//**
 *
 *
 *
 *****************************************************************************/
template<class T>
class chunk_allocator
{

    class chunk  {
    public:
        explicit
        chunk(std::size_t size) noexcept :
            mem_{nullptr}, bof_{nullptr}, end_{nullptr}
        {
            try {
                bof_ = new T[size];
                mem_.reset(bof_);
                end_ = bof_ + size;
            } catch (std::exception&) {
                bof_ = nullptr;
            }
        }

        chunk(const chunk&):
            mem_{nullptr}, bof_{nullptr}, end_{nullptr}
        {}

        chunk(chunk&& src) noexcept :
            mem_{std::move(src.mem_)}, bof_{src.bof_}, end_{src.end_}
        {
            src.bof_ = nullptr;
            src.end_ = nullptr;
        }

        chunk& operator = (const chunk& src) = delete;

        chunk& operator = (chunk&& src) noexcept {
            mem_.swap(src.mem_);
            std::swap(src.bof_, bof_);
            std::swap(src.end_, end_);
            return *this;
        }

        T* begin()      const noexcept { return mem_.get(); }
        T* begin_free() const noexcept { return bof_; }
        T* end()        const noexcept { return end_; }

        std::size_t total_size() const noexcept { return end() - begin(); }
        std::size_t used_size()  const noexcept { return begin_free() - begin(); }
        std::size_t free_size()  const noexcept { return end() - begin_free(); }

        bool owns(const T* p) const noexcept { return p >= begin() && p < end(); }

        T* next_buffer(std::size_t n) noexcept {
            if(free_size() < n) return nullptr;
            auto p = bof_;
            bof_ += n;
            return p;
        }

    private:
        std::unique_ptr<T[]> mem_;
        T* bof_;
        T* end_;
    };

public:
    using value_type = T;

    chunk_allocator():
        minChunkSize_(128*1024*1024/sizeof(T)), //128 MiB
        freeSize_(0),
        chunks_{}
    {}

    chunk_allocator(const chunk_allocator& src):
        minChunkSize_(src.minChunkSize_),
        freeSize_(),
        chunks_{}
    {}

    chunk_allocator& operator = (const chunk_allocator& src) {
        minChunkSize_ = src.minChunkSize_;
        return *this;
    }

    chunk_allocator(chunk_allocator&&) = default;
    chunk_allocator& operator = (chunk_allocator&&) = default;


    void min_chunk_size(std::size_t n) {
//        std::lock_guard<std::mutex> lock(mutables_);
        minChunkSize_ = n;
    }
    std::size_t min_chunk_size() const noexcept {
//        std::lock_guard<std::mutex> lock(mutables_);
        return minChunkSize_;
    }

    bool reserve(std::size_t total)
    {
//        std::lock_guard<std::mutex> lock(mutables_);
        if(total > freeSize_) {
            chunks_.emplace_back(total - freeSize_);
            freeSize_ += chunks_.back().free_size();
        }
        return (freeSize_ >= total);
    }

    T* allocate(std::size_t n)
    {
//        std::lock_guard<std::mutex> lock(mutables_);
        //at the moment chunks will only be used,
        //if they have been reserved explicitly
        if(n <= freeSize_) {
            for(auto& c : chunks_) {
                auto p = c.next_buffer(n);
                if(p) {
                    freeSize_ -= n;
                    return p;
                }
            }
        }
        //make new chunk
//        chunks_.emplace_back(std::max(minChunkSize_,n));
//        auto p = chunks_.back().next_buffer(n);
//        if(p) return p;
        //fallback
        try {
            auto p = new T[n];
            return p;
        } catch(std::exception&) {
            return nullptr;
        }
    }

    void deallocate(T* p, std::size_t)
    {
//        std::lock_guard<std::mutex> lock(mutables_);

        //at the moment occupied chunk buffers are not given back
        auto it = std::find_if(begin(chunks_), end(chunks_),
                               [p](const chunk& c){ return c.owns(p); });

        if(it == end(chunks_)) {
            delete[] p;
        }
    }

    chunk_allocator
    select_on_container_copy_construction() const {
        //don't propagate
        return chunk_allocator{};
    }

private:
//    std::mutex mutables_;
    std::size_t minChunkSize_;
    std::size_t freeSize_;
    std::vector<chunk> chunks_;
};




/*************************************************************************//**
 *
 *
 *
 *****************************************************************************/
template<class T>
bool operator == (const chunk_allocator<T>& a, const chunk_allocator<T>& b) {
    return (&a == &b);
}

//-------------------------------------------------------------------
template<class T>
bool operator != (const chunk_allocator<T>& a, const chunk_allocator<T>& b) {
    return !(a == b);
}


} // namespace mc

#endif

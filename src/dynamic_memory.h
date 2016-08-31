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

/*****************************************************************************
 *
 * @file    dynamci_memory.h
 *
 * @brief   a (very much incomplete and quick & dirty) implementation
 *          of C++17's upcoming memory_resource and polymorphic_allocator
 *          classes together with the memory resource types that are needed
 *          for our hash_multimap
 *
 * @details allocators allow us to have different memory allocation strategies
 *          and thus runtume vs. memory consuption tradeoffs
 *          for the hash table build phase and the query phase
 *
 *****************************************************************************/

#ifndef MC_DYNAMIC_MEMORY_H_
#define MC_DYNAMIC_MEMORY_H_

#include <cstddef>
#include <utility>


namespace mc {


/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
class memory_resource
{
public:
    virtual ~memory_resource() = default;

    void*
    allocate(std::size_t bytes,
             std::size_t alignment = alignof(std::max_align_t))
    {
        return do_allocate(bytes, alignment);
    }

    void*
    deallocate(void* p, std::size_t bytes,
               std::size_t alignment = alignof(std::max_align_t))
    {
        return do_deallocate(p, bytes, alignment);
    }

    bool is_equal(const memory_resource& other) const {
        return do_is_equal(other);
    }

    friend bool
    operator == (const memory_resource& a, const memory_resource& b) {
        return a.is_equal(b);
    }
    friend bool
    operator != (const memory_resource& a, const memory_resource& b) {
        return !a.is_equal(b);
    }
protected:
    virtual void* do_allocate(std::size_t bytes, std::size_t alignment) = 0;

    virtual void* do_deallocate(void* p, std::size_t bytes, std::size_t alignment) = 0;

    virtual bool do_is_equal(const memory_resource&) const = 0;
};



/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
class new_delete_resource : public memory_resource
{
public:

protected:
    void* do_allocate(std::size_t bytes, std::size_t) override
    {
        return reinterpret_cast<void*>(new char[bytes]);
    }

    void* do_deallocate(void* p, std::size_t, std::size_t) override
    {
        delete[] reinterpret_cast<char*>(p);
        return nullptr;
    }

    bool do_is_equal(const memory_resource& other) const override
    {
        return std::addressof(*this) == std::addressof(other);
    }
};



/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
template<class T>
class polymorphic_allocator
{
public:
    using value_type = T;

    polymorphic_allocator():
        mr_{}
    {}

    polymorphic_allocator(memory_resource* mr):
        mr_{mr}
    {}

    T* allocate(std::size_t n)
    {
        return reinterpret_cast<T*>(mr_->allocate(n * sizeof(T)));
    }

    void deallocate(T* p, std::size_t n)
    {
        mr_->deallocate(p, n * sizeof(T));
    }

    template<class U, class... Args>
    void construct(U* p, Args&&... args)
    {
        new (reinterpret_cast<void*>(p)) U(std::forward<Args>(args)...);
    }

    template<class U>
    void destroy(U* p) {
        p->~U();
    }

    polymorphic_allocator
    select_on_container_copy_construction() const {
        //don't propagate
        return polymorphic_allocator{};
    }

    memory_resource* resource() const { return mr_; }

    friend bool
    operator == (const polymorphic_allocator& a, const polymorphic_allocator& b) {
        return a.mr_->is_equal(*b.mr_);
    }
    friend bool
    operator != (const polymorphic_allocator& a, const polymorphic_allocator& b) {
        return !(a == b);
    }

private:
    memory_resource* mr_;
};



} // namespace mc

#endif

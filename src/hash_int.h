#ifndef MC_HASHES_H_
#define MC_HASHES_H_


#include <cstdint>


namespace mc {


//-------------------------------------------------------------------
inline std::uint32_t
thomas_mueller_hash(std::uint32_t x) noexcept {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x);
    return x;
}

//makes sure we cant't use the wrong types
inline void thomas_mueller_hash(std::uint64_t) = delete;
inline void thomas_mueller_hash(std::uint16_t) = delete;
inline void thomas_mueller_hash(std::uint8_t) = delete;



//-------------------------------------------------------------------
inline std::uint32_t
nvidia_hash(std::uint32_t x) noexcept {
    x = (x + 0x7ed55d16) + (x << 12);
    x = (x ^ 0xc761c23c) ^ (x >> 19);
    x = (x + 0x165667b1) + (x <<  5);
    x = (x + 0xd3a2646c) ^ (x <<  9);
    x = (x + 0xfd7046c5) + (x <<  3);
    x = (x ^ 0xb55a4f09) ^ (x >> 16);
    return x;
}

//makes sure we cant't use the wrong types
inline void nvidia_hash(std::uint64_t) = delete;
inline void nvidia_hash(std::uint16_t) = delete;
inline void nvidia_hash(std::uint8_t) = delete;



//-------------------------------------------------------------------
inline std::uint32_t
murmur3_int_init(std::uint32_t x) noexcept {
    x ^= x >> 16;
    x *= 0x85ebca6b;
    x ^= x >> 13;
    x *= 0xc2b2ae35;
    x ^= x >> 16;
    return x;
}

//makes sure we cant't use the wrong types
inline void murmur3_int_init(std::uint64_t) = delete;
inline void murmur3_int_init(std::uint16_t) = delete;
inline void murmur3_int_init(std::uint8_t) = delete;



//-------------------------------------------------------------------
inline std::uint64_t
murmur_hash3_finalizer(std::uint64_t v) noexcept {
    v ^= v >> 33;
    v *= 0xff51afd7ed558ccd;
    v ^= v >> 33;
    v *= 0xc4ceb9fe1a85ec53;
    v ^= v >> 33;
    return v;
}

//makes sure we cant't use the wrong types
inline void murmur_hash3_finalizer(std::uint32_t) = delete;
inline void murmur_hash3_finalizer(std::uint16_t) = delete;
inline void murmur_hash3_finalizer(std::uint8_t) = delete;



/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
template<class T>
struct default_hash;

template<>
struct default_hash<std::uint32_t> {
    std::uint32_t operator () (std::uint32_t x) const noexcept {
        return thomas_mueller_hash(x);
    }
};

template<>
struct default_hash<std::uint64_t> {
    std::uint64_t operator () (std::uint64_t x) const noexcept {
        return murmur_hash3_finalizer(x);
    }
};



/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
struct identity_hash
{
    template<class T>
    T&& operator () (T&& x) const noexcept {
        return std::forward<T>(x);
    }
};


} // namespace mc


#endif

#ifndef MC_HASHES_H_
#define MC_HASHES_H_


#include <cstdint>


namespace mc {

namespace { //internal linkage

/*****************************************************************************
 *
 * @brief integer hash: 32 bits -> 32 bits
 *        invented by Thomas Mueller:
 *        http://stackoverflow.com/users/382763/thomas-mueller
 *
 *****************************************************************************/
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



/*****************************************************************************
 *
 * @brief integer hash: 32 bits -> 32 bits
 *
 *****************************************************************************/
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



/*****************************************************************************
 *
 * @brief integer hash: 32 bits -> 32 bits
 *
 *****************************************************************************/
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



/*****************************************************************************
 *
 * @brief integer hash: 64 bits -> 64 bits
 *
 *****************************************************************************/
inline std::uint64_t
murmur_hash3_finalizer(std::uint64_t x) noexcept {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccd;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53;
    x ^= x >> 33;
    return x;
}

//makes sure we cant't use the wrong types
inline void murmur_hash3_finalizer(std::uint32_t) = delete;
inline void murmur_hash3_finalizer(std::uint16_t) = delete;
inline void murmur_hash3_finalizer(std::uint8_t) = delete;



/*****************************************************************************
 *
 * @brief integer hash: 64 bits -> 64 bits
 *        based on splitmix64 by Sebastiano Vigna (vigna@acm.org)
 *
 *****************************************************************************/
inline std::uint64_t
splitmix64_hash(std::uint64_t x) noexcept
{
    x = (x ^ (x >> 30)) * std::uint64_t(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * std::uint64_t(0x94d049bb133111eb);
    x = x ^ (x >> 31);
    return x;
}

//makes sure we cant't use the wrong types
inline void splitmix64_hash(std::uint32_t) = delete;
inline void splitmix64_hash(std::uint16_t) = delete;
inline void splitmix64_hash(std::uint8_t) = delete;



/*****************************************************************************
 *
 * @brief integer "down" hash: 64 bits -> 32 bits
 *
 *****************************************************************************/
inline std::uint32_t
hash_down(std::uint64_t x) noexcept
{
    x = (~x) + (x << 18);
    x = x ^ (x >> 31);
    x = x * 21;
    x = x ^ (x >> 11);
    x = x + (x << 6);
    x = x ^ (x >> 22);
    return std::uint32_t(x);
}

//makes sure we cant't use the wrong types
inline void hash_down(std::uint32_t) = delete;
inline void hash_down(std::uint16_t) = delete;
inline void hash_down(std::uint8_t) = delete;

}  //internal linkage



/*****************************************************************************
 *
 * @brief for testing purposes
 *
 *****************************************************************************/
struct identity_hash
{
    template<class T>
    T&& operator () (T&& x) const noexcept {
        return std::forward<T>(x);
    }
};



/*****************************************************************************
 *
 * @brief same number of output bits as input bits
 *
 *****************************************************************************/
template<class T>
struct same_size_hash;

template<>
struct same_size_hash<std::uint32_t> {
    std::uint32_t operator () (std::uint32_t x) const noexcept {
        return thomas_mueller_hash(x);
    }
};

template<>
struct same_size_hash<std::uint64_t> {
    std::uint64_t operator () (std::uint64_t x) const noexcept {
        return murmur_hash3_finalizer(x);
    }
};



/*****************************************************************************
 *
 * @brief output bits are half of the input bits
 *
 *****************************************************************************/
template<class T>
struct to32bits_hash;

template<>
struct to32bits_hash<std::uint32_t> {
    std::uint32_t operator () (std::uint32_t x) const noexcept {
        return thomas_mueller_hash(x);
    }
};

template<>
struct to32bits_hash<std::uint64_t> {
    std::uint32_t operator () (std::uint64_t x) const noexcept {
        return hash_down(murmur_hash3_finalizer(x));
    }
};


} // namespace mc


#endif

#include "gpu_hashmap.cuh"
#include "hash_dna.h"
#include "hash_int.h"
#include "sketch_database.h"

namespace mc {

    //---------------------------------------------------------------
    template<
        class Key,
        class ValueT,
        class Hash,
        class KeyEqual,
        class BucketSizeT
    >
    gpu_hashmap<Key,ValueT,Hash,KeyEqual,BucketSizeT>::gpu_hashmap()
    :
        numKeys_(0), numValues_(0), maxLoadFactor_(default_max_load_factor()),
        hash_{}, keyEqual_{}
    {
        cudaSetDevice(0);
    }

    //-----------------------------------------------------
    template<
        class Key,
        class ValueT,
        class Hash,
        class KeyEqual,
        class BucketSizeT
    >
    gpu_hashmap<Key,ValueT,Hash,KeyEqual,BucketSizeT>::gpu_hashmap(const key_equal& keyComp)
    :
        numKeys_(0), numValues_(0), maxLoadFactor_(default_max_load_factor()),
        hash_{}, keyEqual_{keyComp}
    {
        cudaSetDevice(0);
    }

    //-----------------------------------------------------
    template<
        class Key,
        class ValueT,
        class Hash,
        class KeyEqual,
        class BucketSizeT
    >
    gpu_hashmap<Key,ValueT,Hash,KeyEqual,BucketSizeT>::gpu_hashmap(const hasher& hash,
                  const key_equal& keyComp)
    :
        numKeys_(0), numValues_(0), maxLoadFactor_(default_max_load_factor()),
        hash_{hash}, keyEqual_{keyComp}
    {
        cudaSetDevice(0);
    }

    template class gpu_hashmap<
            unsigned int,
            // sketch_database<
            //     std::__cxx11::basic_string<
            //         char,
            //         std::char_traits<char>,
            //         std::allocator<char> >,
            //     single_function_unique_min_hasher<
            //         unsigned int, same_size_hash<unsigned int> >,
            //     same_size_hash<unsigned int>,
            //     unsigned short,
            //     unsigned int,
            //     unsigned char
            //     >::target_location,
            uint64_t,
            same_size_hash<unsigned int>,
            std::equal_to<unsigned int>,
            unsigned char
            >;

} // namespace mc


#include "../src/hash_multimap.h"

#include "../src/timer.h"

#include <stdexcept>
#include <random>


namespace mc_test {

using namespace mc;


//-------------------------------------------------------------------
template<class T, class URNG>
std::vector<T>
make_random_keys(int n, URNG& urng)
{
    auto keys = std::vector<T>{};
    keys.resize(n);
    auto keyDistr = std::uniform_int_distribution<T>{};
    std::generate(begin(keys), end(keys), [&]{ return keyDistr(urng); });
    return keys;
}



//-------------------------------------------------------------------
template<class HashMultiMap>
void hash_multimap_retrieval_correctness(HashMultiMap&& hm)
{
    using key_t   = typename HashMultiMap::key_type;
    using value_t = typename HashMultiMap::mapped_type;

    constexpr int n = 1000000;
    const auto maxBucketSize = hm.max_bucket_size() - 1;

    //make keys
    std::mt19937_64 urng;
    auto keys = make_random_keys<key_t>(n, urng);

    if(keys.empty()) {
        throw std::runtime_error{
            "hash_multimap_retrieval_correctness: no keys for testing"};
    }

    //insert
    auto valueDistr = std::uniform_int_distribution<value_t>{};
    for(const auto& k : keys) {
        auto it = hm.insert(k, valueDistr(urng));
        if(it->size() > maxBucketSize) hm.limit(it,maxBucketSize);
    }

    //query
    for(const auto& k : keys) {
        auto it = hm.find(k);
        if(it == hm.end()) {
            throw std::runtime_error{
                "hash_multimap: inserted element was not found!"};
        }
    }

    //erase
    auto m = std::size_t(keys.size() / 4);
    decltype(keys) erased;
    erased.reserve(m);
    for(std::size_t i = 0; i < m; ++i) {
        erased.push_back(keys.back());
        hm.erase(keys.back());
        keys.pop_back();
    }

    //query erased
    for(const auto& k : erased) {
        auto it = hm.find(k);
        if(it != hm.end()) {
            throw std::runtime_error{
                "hash_multimap: erased element could still be found!"};
        }
    }
}



//-------------------------------------------------------------------
template<class HashMultiMap>
void hash_multimap_retrieval_performance(HashMultiMap&& hm, int n)
{
    using key_t   = typename HashMultiMap::key_type;
    using value_t = typename HashMultiMap::mapped_type;

    const auto maxBucketSize = hm.max_bucket_size() - 1;

    //make keys
    std::mt19937_64 urng;
    auto keys = make_random_keys<key_t>(n, urng);

    //insert
    auto valueDistr = std::uniform_int_distribution<value_t>{};
    timer time;
    time.start();
    for(const auto& k : keys) {
        auto it = hm.insert(k, valueDistr(urng));
        if(it->size() > maxBucketSize) hm.limit(it,maxBucketSize);
    }
    time.stop();
    std::cout << "    insertion: " << time.milliseconds() << " ms" << std::flush;

    //query
    time.restart();
    value_t v = 0;
    for(const auto& k : keys) {
        auto it = hm.find(k);
        if(it != hm.end()) v += (*it)[0];
    }
    time.stop();
    std::cout << "           (ignore this: " << v << ")\n"
              << "    queries: " << time.milliseconds() << " ms" << std::endl;
}



//-------------------------------------------------------------------
void hash_multimap_correctness()
{
    hash_multimap_retrieval_correctness(hash_multimap<uint64_t,uint32_t>{});
    hash_multimap_retrieval_correctness(hash_multimap<uint32_t,uint32_t>{});
    hash_multimap_retrieval_correctness(hash_multimap<uint16_t,uint32_t>{});
    hash_multimap_retrieval_correctness(hash_multimap<uint8_t,uint32_t>{});

    hash_multimap_retrieval_correctness(hash_multimap<uint64_t,uint64_t>{});
    hash_multimap_retrieval_correctness(hash_multimap<uint32_t,uint64_t>{});
    hash_multimap_retrieval_correctness(hash_multimap<uint16_t,uint64_t>{});
    hash_multimap_retrieval_correctness(hash_multimap<uint8_t,uint64_t>{});
}



//-------------------------------------------------------------------
void hash_multimap_performance()
{
    constexpr int n = 1000000;

    std::cout << "test performance with " << n << " keys:\n";

    std::cout << "key: 64 bits, values: 32 bits" << std::endl;
    hash_multimap_retrieval_performance(hash_multimap<uint64_t,uint32_t>{},n);
    std::cout << "key: 32 bits, values: 32 bits" << std::endl;
    hash_multimap_retrieval_performance(hash_multimap<uint32_t,uint32_t>{},n);
    std::cout << "key: 16 bits, values: 32 bits" << std::endl;
    hash_multimap_retrieval_performance(hash_multimap<uint16_t,uint32_t>{},n);
    std::cout << "key: 8 bits, values: 32 bits" << std::endl;
    hash_multimap_retrieval_performance(hash_multimap<uint8_t,uint32_t>{},n);

    std::cout << "key: 64 bits, values: 64 bits" << std::endl;
    hash_multimap_retrieval_performance(hash_multimap<uint64_t,uint64_t>{},n);
    std::cout << "key: 32 bits, values: 64 bits" << std::endl;
    hash_multimap_retrieval_performance(hash_multimap<uint32_t,uint64_t>{},n);
    std::cout << "key: 16 bits, values: 64 bits" << std::endl;
    hash_multimap_retrieval_performance(hash_multimap<uint16_t,uint64_t>{},n);
    std::cout << "key: 8 bits, values: 64 bits" << std::endl;
    hash_multimap_retrieval_performance(hash_multimap<uint8_t,uint64_t>{},n);
}


}

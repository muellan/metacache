
#include "../src/hash_multimap.h"

#include "../src/stat_moments.h"
#include "../src/timer.h"

#include <stdexcept>
#include <random>
#include <fstream>


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
template<class HashMultimap, class OStream = std::ostream>
void print_hash_multimap(const HashMultimap& hm, OStream& os = std::cout)
{
    int i = 0;
    for(const auto& b : hm) {
        os << i << ": \t" << b.key() << " [" << int(b.probe_length()) << "] = \t";
        for(const auto& v : b) {
            os << v << " ";
        }
        os << '\n';
        ++i;
    }
    os << std::endl;
}



//-------------------------------------------------------------------
template<class HashMultiMap, class ValueDistr>
void hash_multimap_correctness(HashMultiMap&& hm, int n, ValueDistr valueDistr)
{
    using key_t   = typename std::decay<HashMultiMap>::type::key_type;

//    std::cout << "------------------" << std::endl;

    //make keys
    std::mt19937_64 urng;
//    urng.seed(std::chrono::system_clock::now().time_since_epoch().count());
    auto inkeys = make_random_keys<key_t>(n, urng);

    if(inkeys.empty()) {
        throw std::runtime_error{
            "hash_multimap_retrieval_correctness: abort - no keys for testing!"};
    }

    //insert
//    hm.reserve_keys(n);
//    hm.reserve_values(n);

    std::cout << "insert & query" << std::endl;
    decltype(inkeys) keys;
    keys.reserve(inkeys.size());
    for(const auto& k : inkeys) {
        auto it = hm.insert(k, valueDistr(urng));
        //if insertion failed due to the limited bucket size, ignore key
        if(it != hm.end()) keys.push_back(k);
    }
    inkeys.clear();

//    hm.rehash(1 + (hm.key_count() / hm.max_load_factor()));
//    print_hash_multimap(hm);

    //query
    for(const auto& k : keys) {
        auto it = hm.find(k);
        if(it == hm.end()) {
//            std::ofstream os{"error.log"};
//            os << "queried for " << k << " in: " << std::endl;
//            print_hash_multimap(hm, os);
            throw std::runtime_error{
                "hash_multimap: inserted key was not found!"};
        }
    }

    //rehash to very high load factor and query again
    std::cout << "rehash to very high load factor & query" << std::endl;
    hm.max_load_factor(0.999);
    hm.rehash(1 + (hm.key_count() / hm.max_load_factor()));
    for(const auto& k : keys) {
        auto it = hm.find(k);
        if(it == hm.end()) {
            throw std::runtime_error{
                "hash_multimap: inserted key was not found after rehash (very high load factor)!"};
        }
    }
    //relax a bit
    hm.max_load_factor(0.95);
    hm.rehash(1 + (hm.key_count() / hm.max_load_factor()));
    std::cout << "rehash to high load factor & query" << std::endl;
    for(const auto& k : keys) {
        auto it = hm.find(k);
        if(it == hm.end()) {
            throw std::runtime_error{
                "hash_multimap: inserted key was not found after rehash (high load factor)!"};
        }
    }

//    print_hash_multimap(hm);

    //erase
    std::cout << "erase & query" << std::endl;
    auto m = std::size_t(keys.size() / 4);
    auto erased = decltype(keys){};
    erased.reserve(m);
    for(std::size_t i = 0; i < m; ++i) {
//        std::cout << "  erase " << keys.back() << std::endl;
        auto k = keys.back();
        erased.push_back(k);
        hm.erase(k);
        //remove all occurrences of k in key vector
        keys.erase(std::remove(begin(keys), end(keys), k), end(keys));
    }
    //query erased
    for(const auto& k : erased) {
        auto it = hm.find(k);
        if(it != hm.end()) {
            throw std::runtime_error{
                "hash_multimap: erased key could still be found!"};
        }
    }
    //query non-erased
    for(const auto& k : keys) {
        auto it = hm.find(k);
        if(it == hm.end()) {
            throw std::runtime_error{
                "hash_multimap: inserted key was not found after erasing others!"};
        }
    }

    //copy & query
    {
        std::cout << "copy & query" << std::endl;
        auto hm2(hm);
        for(const auto& k : keys) {
            auto it = hm2.find(k);
            if(it == hm2.end()) {
                throw std::runtime_error{
                    "hash_multimap: key was not found any more after copying!"};
            }
        }
    }

    std::cout << "clear & query" << std::endl;
    hm.clear();
    //query again
    for(const auto& k : keys) {
        auto it = hm.find(k);
        if(it != hm.end()) {
            throw std::runtime_error{
                "hash_multimap: element found after clearing!"};
        }
    }

    //get new keys
    std::cout << "(insert, query, erase) - loop" << std::endl;
    keys = make_random_keys<key_t>(n, urng);
    //insert, query, erase, query
    for(const auto& k : keys) {
        if(hm.insert(k, valueDistr(urng)) != hm.end()) {
        auto it = hm.find(k);
            if(it == hm.end()) {
                throw std::runtime_error{
                    "hash_multimap: key was not found immediatly after insertion!"};
            }
            hm.erase(k);
            it = hm.find(k);
             if(it != hm.end()) {
                throw std::runtime_error{
                     "hash_multimap: key could still be found immediatly after deletion!"};
             }
        }
    }
}



//-------------------------------------------------------------------
template<class HashMultiMap, class ValueDistr>
void hash_multimap_performance(HashMultiMap&& hm, int n, ValueDistr&& valueDistr)
{
    using key_t   = typename HashMultiMap::key_type;
    using value_t = typename HashMultiMap::mapped_type;

    //make keys
    std::mt19937_64 urng;
    auto keys = make_random_keys<key_t>(n, urng);

    //insert
    timer time;
    time.start();
    for(const auto& k : keys) {
        hm.insert(k, valueDistr(urng));
    }
    time.stop();

    //rehash to the worst case with current maximum load factor
    hm.rehash(1 + (hm.key_count() / hm.max_load_factor()));

    std::cout << "    insertion took " << time.milliseconds() << " ms"
              << "\n    keys:    " << hm.key_count()
              << "\n    values:  " << hm.value_count()
              << "\n    buckets: " << hm.bucket_count()
              << "\n    load:    " << (100 * hm.load_factor()) << "%"
              << std::flush;

    variance_accumulator<double> bs;
    variance_accumulator<double> pl;
    for(const auto& b : hm) {
        bs += b.size();
        pl += b.probe_length();
    }
    std::cout << "\n    bucket sizes = " << bs.mean() << " +/- " << bs.stddev();
    std::cout << "\n    probe lengths = " << pl.mean() << " +/- " << pl.stddev();

    //query
    time.restart();
    value_t v = 0;
    for(const auto& k : keys) {
        auto it = hm.find(k);
        if(it != hm.end()) v += (*it)[0];
    }
    time.stop();
    std::cout << "           (ignore this: " << v << ")\n"
              << "    queries took " << time.milliseconds() << " ms "
              << " => " << (keys.size()/time.seconds()) << " queries/s"
              << std::endl;
}



//-------------------------------------------------------------------
void hash_multimap_correctness()
{
    constexpr int n = 100000;

    auto vd32 = std::uniform_int_distribution<uint32_t>{};
    hash_multimap_correctness(hash_multimap<uint32_t,uint32_t>{},n,vd32);
    hash_multimap_correctness(hash_multimap<uint64_t,uint32_t>{},n,vd32);

    auto vd64 = std::uniform_int_distribution<uint64_t>{};
    hash_multimap_correctness(hash_multimap<uint32_t,uint64_t>{},n,vd64);
    hash_multimap_correctness(hash_multimap<uint64_t,uint64_t>{},n,vd64);


    { //reserve space for all half of the values upfront
        hash_multimap<uint32_t,uint32_t> hm{};
        hm.reserve_values(n/2);
        hash_multimap_correctness(hm,n,vd32);
    }
    { //reserve space for all half of the values upfront
        hash_multimap<uint32_t,uint32_t> hm{};
        hm.reserve_keys(n/2);
        hm.reserve_values(n/2);
        hash_multimap_correctness(hm,n,vd32);
    }
    { //reserve space for all values upfront
        hash_multimap<uint32_t,uint32_t> hm{};
        hm.reserve_values(n);
        hash_multimap_correctness(hm,n,vd32);
    }
    { //reserve space for all values and all keys upfront
        hash_multimap<uint32_t,uint32_t> hm{};
        hm.reserve_keys(n);
        hm.reserve_values(n);
        hash_multimap_correctness(hm,n,vd32);
    }
}



//-------------------------------------------------------------------
void hash_multimap_performance()
{
    constexpr int n = 1000000;

    std::cout << "test performance with " << n << " keys:\n";

    auto vd32 = std::uniform_int_distribution<uint32_t>{};
    std::cout << "key: 64 bits, values: 32 bits" << std::endl;
    hash_multimap_performance(hash_multimap<uint64_t,uint32_t>{},n,vd32);
    std::cout << "key: 32 bits, values: 32 bits" << std::endl;
    hash_multimap_performance(hash_multimap<uint32_t,uint32_t>{},n,vd32);

    auto vd64 = std::uniform_int_distribution<uint64_t>{};
    std::cout << "key: 64 bits, values: 64 bits" << std::endl;
    hash_multimap_performance(hash_multimap<uint64_t,uint64_t>{},n,vd64);
    std::cout << "key: 32 bits, values: 64 bits" << std::endl;
    hash_multimap_performance(hash_multimap<uint32_t,uint64_t>{},n,vd64);
}


}

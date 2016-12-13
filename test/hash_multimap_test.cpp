
#include "../src/hash_multimap.h"

#include "../src/stat_moments.h"
#include "../src/timer.h"

#include <stdexcept>
#include <random>
#include <fstream>
#include <vector>
#include <set>


namespace mc_test {

using namespace mc;


//-------------------------------------------------------------------
template<class Key, class Value>
class key_value_pair_filler
{
public:
    using key_type = Key;
    using value_type = Value;
    using key_distr = std::uniform_int_distribution<Key>;
    using value_distr = std::uniform_int_distribution<Value>;
    using multi_distr = std::uniform_int_distribution<std::size_t>;
    using result_type = std::vector<std::pair<Key,Value>>;

    explicit
    key_value_pair_filler(std::size_t seed = 0):
        urng_(seed), keyDistr_{}, valueDistr_{},
        multiDistr_{1,10}
    {
//    urng_.seed(std::chrono::system_clock::now().time_since_epoch().count());
    }

    void key_range(key_type min, key_type max) {
        keyDistr_ = key_distr{min,max};
    }

    void value_range(const value_type& min, const value_type& max) {
        valueDistr_ = value_distr{min,max};
    }

    void values_per_key(std::size_t min, std::size_t max) {
        if(min < 0) min = 1;
        if(max < 0) max = 1;
        if(min > max) std::swap(min,max);
        multiDistr_ = multi_distr{min,max};
    }

    template<class HashMap>
    result_type operator () (std::size_t n, HashMap& hm)
    {
        //insert
    //    hm.reserve_keys(n);
    //    hm.reserve_values(n);

        std::vector<std::pair<Key,Value>> kvpairs;
        kvpairs.reserve(n);
        for(std::size_t i = 0; i < n; ++i) {
            auto k = keyDistr_(urng_);
            auto m = multiDistr_(urng_);
            for(std::size_t j = 0; j < m; ++j) {
                auto v = valueDistr_(urng_);
                auto it = hm.insert(k, v);
                //if insertion failed due to the limited bucket size, ignore (key,value)
                if(it != hm.end()) kvpairs.emplace_back(k,v);
            }
        }

        return kvpairs;
    }


private:
    std::mt19937_64 urng_;
    key_distr keyDistr_;
    value_distr valueDistr_;
    multi_distr multiDistr_;
};



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
template<class HashMultiMap>
void hash_multimap_check_binary_IO(const HashMultiMap& hm)
{
    {
        std::ofstream os {"test.map", std::ios::binary};
        write_binary(os, hm);
    }

    HashMultiMap hm2;
    {
        std::ifstream is {"test.map", std::ios::binary};
        read_binary(is, hm2);
    }

    if(hm.key_count() != hm2.key_count()) {
        std::ofstream os1 {"test1.out"};
        print_hash_multimap(hm, os1);

        std::ofstream os2 {"test2.out"};
        print_hash_multimap(hm2, os2);

        throw std::runtime_error{
            "hash_multimap::key_count() inconsistent after deserialization"};
    }
    if(hm.value_count() != hm2.value_count()) {

        throw std::runtime_error{
            "hash_multimap::value_count() inconsistent after deserialization"};
    }
}



//-------------------------------------------------------------------
template<class HashMultiMap, class K, class V>
void hash_multimap_check_book_keeping(
    HashMultiMap&& hm,
    const std::vector<std::pair<K,V>>& kvpairs,
    const std::string& message = "")
{
    using key_t = typename std::decay<HashMultiMap>::type::key_type;

    if(kvpairs.size() != hm.value_count()) {
        std::cout << kvpairs.size() << " != " << hm.value_count() << std::endl;
        throw std::runtime_error{
            "hash_multimap::value_count() incorrect " + message + "!"};
    }

    //make set of unique keys
    std::set<key_t> ukeys;
    for(const auto& p : kvpairs) ukeys.insert(p.first);

    if(ukeys.size() != hm.key_count()) {
        std::cout << ukeys.size() << " != " << hm.key_count() << std::endl;
        throw std::runtime_error{
            "hash_multimap::key_count() incorrect " + message + "!"};
    }
}



//-------------------------------------------------------------------
template<class HashMultiMap, class K, class V>
void hash_multimap_check_presence(HashMultiMap&& hm,
                                  const std::vector<std::pair<K,V>>& kvpairs,
                                  const std::string& message = "")
{
    hash_multimap_check_book_keeping(hm, kvpairs, message);

    for(const auto& p : kvpairs) {
        auto it = hm.find(p.first);
        if(it == hm.end()) {

            std::ofstream os{"error.log"};
            os << "queried for " << p.first << " in: " << std::endl;
            print_hash_multimap(hm, os);

//            std::cout << "queried for " << p.first << " in: " << std::endl;
//            print_hash_multimap(hm, std::cout);

            throw std::runtime_error{
                "hash_multimap: inserted key was not found " + message};
        }
        //check if value is there
        if(std::find(it->begin(), it->end(), p.second) == it->end()) {
            throw std::runtime_error{
                "hash_multimap: inserted value was not found " + message};
        }
    }
}



//-------------------------------------------------------------------
template<class HashMultiMap, class K>
void hash_multimap_check_absence(HashMultiMap&& hm,
                                 const std::vector<K>& keys,
                                 const std::string& message = "")
{
    for(const auto& k : keys) {
        auto it = hm.find(k);
        if(it != hm.end()) {
//            std::ofstream os{"error.log"};
//            os << "queried for " << k << " in: " << std::endl;
//            print_hash_multimap(hm, os);
            throw std::runtime_error{
                "hash_multimap: unexpected key found " + message};
        }
    }
}



//-------------------------------------------------------------------
template<class HashMultiMap, class KeyValGen>
void hash_multimap_correctness(HashMultiMap&& hm, std::size_t n, KeyValGen&& keyValGen)
{

//    std::cout << "------------------" << std::endl;

//    std::cout << "insert & query" << std::endl;
    auto kvpairs = keyValGen(n, hm);
    hash_multimap_check_presence(hm, kvpairs, "after initial insert");

    //rehash to very high load factor and query again
//    std::cout << "rehash to very high load factor & query" << std::endl;
    hm.max_load_factor(0.999);
    hm.rehash(1 + (hm.key_count() / hm.max_load_factor()));
//    print_hash_multimap(hm);
    hash_multimap_check_presence(hm, kvpairs, "after rehash to very high load factor");

    //relax a bit
//    std::cout << "rehash to high load factor & query" << std::endl;
    hm.max_load_factor(0.95);
    hm.rehash(1 + (hm.key_count() / hm.max_load_factor()));
    hash_multimap_check_presence(hm, kvpairs, "after rehash to high load factor");

    //copy & query
    {
//        std::cout << "copy & query" << std::endl;
        auto hm2(hm);
        hash_multimap_check_presence(hm2, kvpairs, "after copying");
    }

    //erase
/*
//    std::cout << "erase & query" << std::endl;
    using key_t = typename std::decay<HashMultiMap>::type::key_type;
    using val_t = typename std::decay<HashMultiMap>::type::value_type;

    auto m = std::size_t(kvpairs.size() / 3);
    std::vector<key_t> erased;
    erased.reserve(m);
    std::size_t numerased = 0;
    auto kvcand = kvpairs;
    auto numvalsold = hm.value_count();
    for(std::size_t i = 0; i < m; ++i) {
//        std::cout << "  erase " << keys.back() << std::endl;
        auto k = kvcand.back().first;
        kvcand.pop_back();
        auto s = hm.erase(k);
        if(s > 0) {
            numerased += s;
            erased.push_back(k);
            kvpairs.erase(std::remove_if(begin(kvpairs), end(kvpairs),
                           [=](const std::pair<key_t,val_t>& p) {
                                return p.first == k; }),
                    end(kvpairs));
        }
    }

    if(hm.value_count() != (numvalsold - numerased)) {
        throw std::runtime_error{
            "hash_multimap::erase return value inconsistent with value_count" };
    }

    //query non-erased
    hash_multimap_check_presence(hm, kvpairs, "after erasing others");
    //query erased
    hash_multimap_check_absence(hm, erased, ": was erased before");
*/

    hash_multimap_check_binary_IO(hm);
}



//-------------------------------------------------------------------
template<class HashMultiMap, class KeyValGen>
void hash_multimap_performance(HashMultiMap&& hm, std::size_t n, KeyValGen&& keyValGen)
{
    using value_t = typename std::decay<HashMultiMap>::type::value_type;

    //insert
    timer time;
    time.start();
    auto kvpairs = keyValGen(n, hm);
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
    for(const auto& p : kvpairs) {
        auto it = hm.find(p.first);
        if(it != hm.end()) v += (*it)[0];
    }
    time.stop();
    std::cout << "           (ignore this: " << v << ")\n"
              << "    queries took " << time.milliseconds() << " ms  => "
              << ((kvpairs.size()/time.seconds())/1e6) << " Mqueries/s"
              << std::endl;
}



//-------------------------------------------------------------------
void hash_multimap_correctness()
{
    constexpr std::size_t n = 10000;

    std::cout << "test correctness with " << n << " keys" << std::endl;

    auto k32v32 = key_value_pair_filler<uint32_t,uint32_t>{};
    hash_multimap_correctness(hash_multimap<uint32_t,uint32_t>{},n,k32v32);
    //make multiple non-consecutive inserts on same key more likely
    k32v32.key_range(1, 2*n);
    hash_multimap_correctness(hash_multimap<uint32_t,uint32_t>{},n,k32v32);

    auto k64v64 = key_value_pair_filler<uint64_t,uint64_t>{};
    hash_multimap_correctness(hash_multimap<uint64_t,uint64_t>{},n,k64v64);

//    auto k32v64 = key_value_pair_filler<uint32_t,uint64_t>{};
//    hash_multimap_correctness(hash_multimap<uint32_t,uint64_t>{},n,k32v64);
//    auto k64v32 = key_value_pair_filler<uint64_t,uint32_t>{};
//    hash_multimap_correctness(hash_multimap<uint64_t,uint32_t>{},n,k64v32);


    { //reserve space for all half of the values upfront
        hash_multimap<uint32_t,uint32_t> hm{};
        hm.reserve_values(n/2);
        hash_multimap_correctness(hm,n,k32v32);
    }
    { //reserve space for all half of the values upfront
        hash_multimap<uint32_t,uint32_t> hm{};
        hm.reserve_keys(n/2);
        hm.reserve_values(n/2);
        hash_multimap_correctness(hm,n,k32v32);
    }
    { //reserve space for all values upfront
        hash_multimap<uint32_t,uint32_t> hm{};
        hm.reserve_values(n);
        hash_multimap_correctness(hm,n,k32v32);
    }
    { //reserve space for all values and all keys upfront
        hash_multimap<uint32_t,uint32_t> hm{};
        hm.reserve_keys(n);
        hm.reserve_values(n);
        hash_multimap_correctness(hm,n,k32v32);
    }
}



//-------------------------------------------------------------------
void hash_multimap_performance()
{
    constexpr std::size_t n = 1000000;

    std::cout << "test performance with " << n
              << " keys and around 1-10 values per key:\n";

    //key distributions, make sure we get multiple values per key

    std::cout << "key: 32 bits, values: 32 bits" << std::endl;
    auto k32v32 = key_value_pair_filler<uint32_t,uint32_t>{};
    hash_multimap_performance(hash_multimap<uint32_t,uint32_t>{},n,k32v32);
    std::cout << "key: 32 bits, values: 64 bits" << std::endl;
    auto k32v64 = key_value_pair_filler<uint32_t,uint64_t>{};
    hash_multimap_performance(hash_multimap<uint32_t,uint64_t>{},n,k32v64);

    std::cout << "key: 64 bits, values: 32 bits" << std::endl;
    auto k64v32 = key_value_pair_filler<uint64_t,uint32_t>{};
    hash_multimap_performance(hash_multimap<uint64_t,uint32_t>{},n,k64v32);
    std::cout << "key: 64 bits, values: 64 bits" << std::endl;
    auto k64v64 = key_value_pair_filler<uint64_t,uint64_t>{};
    hash_multimap_performance(hash_multimap<uint64_t,uint64_t>{},n,k64v64);
}


}


#include <iostream>

#include "sketch_database_test.h"

namespace mc_test {

void hash_multimap_correctness();
void hash_multimap_performance();

} // namespace mc_test


int main()
{
    using namespace mc_test;

    try {
        hash_multimap_correctness();
        hash_multimap_performance();

//        evaluate_database_lookup_performane("G2_w128_k16_s16.db",
//                                            "test/random_keys.txt");

//        evaluate_database_lookup_performane("G2_w512_k16_s16.db",
//                                            "test/random_keys.txt");

//        evaluate_database_lookup_performane("G2_w8192_k16_s16.db",
//                                            "test/random_keys.txt");

        return 0;
    }
    catch (std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }
}

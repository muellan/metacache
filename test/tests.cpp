
#include <iostream>

#include "performance.h"

int main()
{
    try {
        mc_test::evaluate_database_lookup_performane("G2_w128_k16_s16.db",
                                                     "test/random_keys.txt");

//        mc_test::evaluate_database_lookup_performane("G2_w512_k16_s16.db",
//                                                     "test/random_keys.txt");

//        mc_test::evaluate_database_lookup_performane("G2_w8192_k16_s16.db",
//                                                     "test/random_keys.txt");

        return 0;
    }
    catch (std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }
}

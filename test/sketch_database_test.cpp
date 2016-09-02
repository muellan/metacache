
#include "../src/config.h"


namespace mc_test {

using namespace mc;


//-------------------------------------------------------------------
void evaluate_database_lookup_performane(const std::string& dbfname,
                                         const std::string& keyfname)
{
    //read database
    auto db = make_database<database>(dbfname);

    //read keys from file
    std::cout << "Reading keys from file " << keyfname << " ... " << std::flush;
    std::ifstream is {keyfname};
    std::vector<std::uint32_t> keys;
    while(is.good()) {
        std::uint32_t k;
        if(is >> k) keys.push_back(k);
    }
    std::cout << " done.\n" << "Running test ... " << std::endl;

    //run test
    std::vector<database::reference_pos> values;
    values.reserve(keys.size() * 1024);
    auto time = db.hash_table_query_benchmark(keys, values);

    std::cout
        << "key count:   " << keys.size() << '\n'
        << "refs found:  " << values.size() << '\n'
        << "query time:  " << time.milliseconds() << " ms\n"
        << "query speed: " << (keys.size() / time.seconds()) << " queries/s"
        << std::endl;
}


} //namespace mc_test

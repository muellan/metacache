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

#include <iostream>
#include <fstream>

#include "args_handling.h"
#include "modes.h"
#include "dna_encoding.h"

#include "timer.h"


/*****************************************************************************
 *
 * @brief  select mode and display usage info
 *
 *****************************************************************************/
int main(int argc, char* argv[])
{
    using namespace mc;

    /*
    //database / hash map query time benchmarking
    auto db = make_database<database_type>("refseq_w128_k16_s16.db");

    std::ifstream is {"random_keys.txt"};

    std::vector<std::uint32_t> keys;

    while(is.good()) {
        std::uint32_t k;
        is >> k;
        keys.push_back(k);
    }

    std::vector<database_type::reference_pos> values;
    values.reserve(keys.size() * 1024);

    auto time = db.hash_table_query_benchmark(keys, values);

    std::cout << time.milliseconds() << " ms " << std::endl;
    */


    args_parser args{argc, argv};
    
    bool nomode = false;
    auto modestr = std::string{};
    if(args.non_prefixed_count() > 0) {
        modestr = args.non_prefixed(0);
    }

    try {
        if(modestr == "help") {
            main_mode_help(args);
        }
        else if(modestr == "build") {
            main_mode_build(args);
        }
        else if(modestr == "add") {
            main_mode_build_add(args);
        }
        else if(modestr == "query") {
            main_mode_query(args);
        }
        else if(modestr == "info") {
            main_mode_info(args);
        }
        else if(modestr == "annotate") {
            main_mode_annotate(args);
        }
        else {
            nomode = true;
        }
    }
    catch(std::runtime_error& e) {
        std::cout << "\nABORT: " << e.what() << "!" << std::endl;
    }
    catch(std::invalid_argument& e) {
        std::cout << "\nABORT: " << e.what() << "!" << std::endl;
        nomode = true;
    }

    if(nomode) {
       main_mode_help(args);
    }

}

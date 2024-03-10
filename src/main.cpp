/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2024 André Müller (muellan@uni-mainz.de)
 *                       & Robin Kobus  (kobus@uni-mainz.de)
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


#include "modes.h"

#include <iostream>
#include <stdexcept>


/*************************************************************************//**
 *
 * @brief  selects & executes main mode / shows quick help
 *
 *****************************************************************************/
int main(int argc, char** argv)
{
    using namespace mc;

    std::string modestr;
    if (argc > 1) { modestr = argv[1]; }

    try {

        if (modestr == "build") {
            main_mode_build(make_args_list(argv+2, argv+argc));
        }
        else if (modestr == "modify") {
            main_mode_modify(make_args_list(argv+2, argv+argc));
        }
        else if (modestr == "query") {
            main_mode_query(make_args_list(argv+2, argv+argc));
        }
        else if (modestr == "build+query") {
            main_mode_build_query(make_args_list(argv+2, argv+argc));
        }
        else if (modestr == "merge") {
            main_mode_merge(make_args_list(argv+2, argv+argc));
        }
        else if (modestr == "info") {
            main_mode_info(make_args_list(argv+2, argv+argc));
        }
        else {
            main_mode_help(make_args_list(argv, argv+argc));
        }
    }
    catch(std::runtime_error& e) {
        std::cerr << "\nABORT: " << e.what() << '\n';
    }
    catch(std::invalid_argument& e) {
        std::cerr << "ERROR: Invalid command line arguments!\n\n"
                  << e.what() << "\n";
    }
    catch(std::exception& e) {
        std::cerr << e.what() << "\n\n";
    }

}

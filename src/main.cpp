/*****************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
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
 * @brief  select mode or display quick help
 *
 *****************************************************************************/
int main(int argc, char* argv[])
{
    using namespace mc;

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
        else if(modestr == "modify") {
            main_mode_build_modify(args);
        }
        else if(modestr == "query") {
            main_mode_query(args);
        }
        else if(modestr == "info") {
            main_mode_info(args);
        }
        else if(modestr == "evaluate") {
            main_mode_evaluate(args);
        }
        else if(modestr == "annotate") {
            main_mode_annotate(args);
        }
        else {
            nomode = true;
        }
    }
    catch(std::runtime_error& e) {
        std::cerr << "\nABORT: " << e.what() << "!" << std::endl;
    }
    catch(std::invalid_argument& e) {
        std::cerr << "\nABORT: " << e.what() << "!" << std::endl;
        nomode = true;
    }

    if(nomode) {
       main_mode_help(args);
    }

}

/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2019 André Müller (muellan@uni-mainz.de)
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

#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include "args_handling.h"
#include "filesys_utility.h"


namespace mc {


//-------------------------------------------------------------------
void show_available_help_topics()
{
    auto files = files_in_directory("docs");

    if(files.empty()) {
        std::cerr << "Documentation files are missing!\n"
                  << "These should be in a folder called 'docs'.\n"
                  << "Please download a new copy of MetaCache!"
                  << std::endl;
        return;
    }

    std::cout << "Available topics are:\n";
    for(const auto& name : files) {
        auto i = name.find(".txt");
        auto topic = name.substr(5,i-5);
        auto spc = std::string(10 - topic.size(), ' ');
        auto synopsis = std::string("");

        std::ifstream is{name};
        if(is.good()) {
            //forward to synopsis
            while(is.good() && synopsis.find("SYNOPSIS:") == std::string::npos) {
                getline(is, synopsis);
            }
            getline(is, synopsis);
        }

        std::cout << "    " << topic << spc << synopsis << '\n';
    }
    std::cout.flush();
}


//-------------------------------------------------------------------
void show_help_for_topic(const std::string& topic)
{
    std::string filename = "docs/" + topic + ".txt";

    std::ifstream is {filename};

    if(is.good()) {
        std::string line;
        do {
            getline(is, line);
            std::cout << line << '\n';
        } while(is.good());

        std::cout.flush();
    }
    else {
        std::cerr << "Documentation for topic '" << topic << "' not available.\n";
        show_available_help_topics();
    }
}


//-------------------------------------------------------------------
void main_mode_help(const args_parser& args)
{
    if(args.non_prefixed_count() < 1) {
        show_help_for_topic("quick");
    }
    else if(args.non_prefixed_count() == 1 ) {
       if(args.non_prefixed(0) == "help") {
           std::cerr << "You need to specify a help topic:\n"
                     << "    ./metacache help <topic>\n";

           show_available_help_topics();
       } else {
           show_help_for_topic("quick");
       }
    }
    else {
        show_help_for_topic(args.non_prefixed(1));
    }

}


} // namespace mc

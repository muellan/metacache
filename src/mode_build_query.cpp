/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2021 André Müller (muellan@uni-mainz.de)
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



#include "building.h"
#include "database.h"
#include "options.h"
#include "querying.h"

#include <iostream>



namespace mc {

using std::cout;
using std::cerr;
using std::endl;




/*************************************************************************//**
 *
 * @brief adds targets to datbase and queries
 *
 *****************************************************************************/
void add_to_database_and_query(database& db, build_query_options& opt)
{
    timer time;
    time.start();

    add_to_database(db, opt.build, time);

    time.stop();

    db.initialize_taxonomy_caches();

    adapt_options_to_database(opt.query, db);

    if (!opt.query.infiles.empty()) {
        cerr << "Classifying query sequences.\n";

        process_input_files(db, opt.query);
    }
    else {
        cout << "No input files provided.\n";

        run_interactive_query_mode(db, opt.query);
    }

    if (opt.saveDatabase) {
        write_database(db, opt.build);
    }

    //prevents slow deallocation
    db.clear_without_deallocation();
}


/*************************************************************************//**
 *
 * @brief builds a database from reference input sequences and query it
 *
 *****************************************************************************/
void main_mode_build_query(const cmdline_args& args)
{
    auto opt = get_build_query_options(args);

    if (opt.build.infoLevel != info_level::silent) {
        cout << "Building new database from reference sequences." << endl;
    }

    auto db = database{opt.build.sketching};

    add_to_database_and_query(db, opt);
}


} // namespace mc

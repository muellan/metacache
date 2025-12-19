/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2026 André Müller (github.com/muellan)
 *                       & Robin Kobus  (github.com/funatiq)
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


#include "building.hpp"
#include "database.hpp"
#include "options.hpp"
#include "timer.hpp"

#include <iostream>


namespace mc {

using std::cout;
using std::endl;


//-----------------------------------------------------------------------------
/**
 * @brief adds targets to datbase and writes to disk
 */
void add_to_database_and_save (database& db, const build_options& opt)
{
    timer time;
    time.start();

    add_to_database(db, opt, time);

    const auto btime = time.seconds();

    write_database(db, opt);

    time.stop();

    const auto wtime = time.seconds() - btime;

    if (opt.infoLevel != info_level::silent) {
        cout << "------------------------------------------------\n"
             << "Construction time: " << btime << " s\n"
             << "Writing time:      " << wtime << " s\n"
             << "Total build time:  " << time.seconds() << " s" << endl;
    }

    // prevents slow deallocation
    db.clear_without_deallocation();
}



//-----------------------------------------------------------------------------
/**
 * @brief adds reference sequences to an existing database
 */
void main_mode_modify (const cmdline_args& args)
{
    auto opt = get_modify_options(args);

    cout << "Modify database " << opt.dbfile << endl;

    auto db = make_database(opt.dbfile, opt.dbpart);

    if (opt.infoLevel != info_level::silent && !opt.infiles.empty()) {
        cout << "Adding reference sequences to database..." << endl;
    }

    add_to_database_and_save(db, opt);
}



//-----------------------------------------------------------------------------
/**
 * @brief builds a database from reference input sequences
 */
void main_mode_build (const cmdline_args& args)
{
    auto opt = get_build_options(args);

    if (opt.infoLevel != info_level::silent) {
        cout << "Building new database '" << opt.dbfile
             << "' from reference sequences." << endl;
    }

    auto db = database{opt.sketching};

    add_to_database_and_save(db, opt);
}


} // namespace mc

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


#include "io_error.hpp"
#include "options.hpp"
#include "querying.hpp"

#include <iostream>


namespace mc {

using std::cout;
using std::cerr;


//-----------------------------------------------------------------------------
/**
 * @brief read database and modify db content and sketching scheme according to
 *        command line options
 */
database
read_database (const query_options& opt)
{
    const database_storage_options& dbopt = opt.dbconfig;

    database db;

    if (dbopt.maxLoadFactor > 0.4 && dbopt.maxLoadFactor < 0.99) {
        db.max_load_factor(dbopt.maxLoadFactor);
        if (opt.output.showInfo) {
            cerr << "Using custom hash table load factor of "
                 << dbopt.maxLoadFactor << '\n';
        }
    }

    try {
        const auto infoLvl = opt.output.showInfo 
                    ? info_level::moderate : info_level::silent;

        db.read(opt.dbfile, opt.dbpart, opt.performance.replication,
                database::scope::everything, infoLvl);
    }
    catch (const file_access_error& e) {
        cerr << "FAIL\n";
        throw;
    }

    if (dbopt.removeOverpopulatedFeatures) {
        auto old = db.feature_count();

        auto maxlpf = dbopt.maxLocationsPerFeature - 1;
        if (maxlpf < 0 || maxlpf >= database::max_supported_locations_per_feature())
            maxlpf = database::max_supported_locations_per_feature() - 1;

        maxlpf = std::min(maxlpf, db.max_locations_per_feature() - 1);
        if (maxlpf > 0) { // always keep buckets with size 1
            if (opt.output.showInfo) {
                cerr << "\nRemoving features with more than "
                     << maxlpf << " locations...\n";
            }
            auto rem = db.remove_features_with_more_locations_than(maxlpf);

            if (opt.output.showInfo) {
                cerr << rem << " of " << old << " removed.\n";
            }
        }
        // in case new max is less than the database setting
        db.max_locations_per_feature(dbopt.maxLocationsPerFeature);
    }
    else if (dbopt.maxLocationsPerFeature > 1) {
        db.max_locations_per_feature(dbopt.maxLocationsPerFeature);
        if (opt.output.showInfo) {
            cerr << "Max locations per feature set to "
                 << dbopt.maxLocationsPerFeature << '\n';
        }
    }

    return db;
}




//-----------------------------------------------------------------------------
/**
 * @brief    run query reads against pre-built database
 *           entry point for query mode
 *
 * @details  note that precision (positive predictive value) and
 *           clade exclusion testing is much slower and only intended
 *           for classification performance evaluation
 */
void main_mode_query (const cmdline_args& args)
{
    auto opt = get_query_options(args);

    auto db = read_database(opt);
    adapt_options_to_database(opt, db);

    if (not opt.infiles.empty()) {
        if (opt.output.showInfo) {
            cerr << "Classifying query sequences.\n";
        }
        process_input_files(db, opt);
    }
    else {
        if (opt.output.showInfo) {
            cout << "No input files provided.\n";
        }
        run_interactive_query_mode(db, opt);
    }
}


} // namespace mc

/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2020 André Müller (muellan@uni-mainz.de)
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


#include "io_error.h"
#include "options.h"
#include "querying.h"

#include <iostream>
#include <stdexcept>


namespace mc {

using std::cout;
using std::cerr;
using std::endl;
using std::flush;


/*************************************************************************//**
 *
 * @brief read database and modify db content and sketching scheme according to
 *        command line options
 *
 *****************************************************************************/
database
read_database(const query_options& opt)
{
    const database_storage_options& dbopt = opt.dbconfig;

    database db;

    if(dbopt.maxLoadFactor > 0.4 && dbopt.maxLoadFactor < 0.99) {
        db.max_load_factor(dbopt.maxLoadFactor);
        cerr << "Using custom hash table load factor of "
             << dbopt.maxLoadFactor << '\n';
    }

    try {
        db.read(opt.dbfile, opt.dbpart, opt.performance.replication);
    }
    catch(const file_access_error& e) {
        cerr << "FAIL" << endl;
        throw;
    }

    if(dbopt.removeOverpopulatedFeatures) {
        auto old = db.feature_count();

        auto maxlpf = dbopt.maxLocationsPerFeature - 1;
        if(maxlpf < 0 || maxlpf >= database::max_supported_locations_per_feature())
            maxlpf = database::max_supported_locations_per_feature() - 1;

        maxlpf = std::min(maxlpf, db.max_locations_per_feature() - 1);
        if(maxlpf > 0) { //always keep buckets with size 1
            cerr << "\nRemoving features with more than "
                 << maxlpf << " locations... " << std::flush;

            auto rem = db.remove_features_with_more_locations_than(maxlpf);

            cerr << rem << " of " << old << " removed.\n";
        }
        //in case new max is less than the database setting
        db.max_locations_per_feature(dbopt.maxLocationsPerFeature);
    }
    else if(dbopt.maxLocationsPerFeature > 1) {
        db.max_locations_per_feature(dbopt.maxLocationsPerFeature);
        cerr << "Max locations per feature set to "
             << dbopt.maxLocationsPerFeature << '\n';
    }

    const sketching_options& skopt = opt.sketching;
    //use a different sketching scheme for querying?
    if((skopt.sketchlen > 0) || (skopt.winlen > 0) || (skopt.winstride > 0)) {
        auto s = db.query_sketcher();
        if(skopt.sketchlen > 0) s.sketch_size(skopt.sketchlen);
        if(skopt.winlen > 0)    s.window_size(skopt.winlen);
        // if no custom window stride requested => set to w-k+1
        auto winstride = skopt.winstride;
        if(winstride < 0) winstride = s.window_size() - s.kmer_size() + 1;
        s.window_stride(winstride);

        cerr << "custom query sketching settings:"
             << " -winlen "    << s.window_size()
             << " -winstride " << s.window_stride()
             << " -sketchlen " << s.sketch_size() << '\n';

        db.query_sketcher(std::move(s));
    }

    return db;
}



/*************************************************************************//**
 *
 * @brief    run query reads against pre-built database
 *           entry point for query mode
 *
 * @details  note that precision (positive predictive value) and
 *           clade exclusion testing is much slower and only intended
 *           for classification performance evaluation
 *
 *****************************************************************************/
void main_mode_query(const cmdline_args& args)
{
    auto opt = get_query_options(args);

    auto db = read_database(opt);

    if(!opt.infiles.empty()) {
        cerr << "Classifying query sequences.\n";

        adapt_options_to_database(opt.classify, db);
        process_input_files(db, opt);
    }
    else {
        cout << "No input files provided.\n";

        run_interactive_query_mode(db, opt);
    }
}


} // namespace mc

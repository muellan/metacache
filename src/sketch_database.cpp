/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2019 André Müller (muellan@uni-mainz.de)
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

#include "sketch_database.h"


namespace mc {


// ----------------------------------------------------------------------------
database
make_database(const std::string& filename, database::scope what, info_level info)
{
    if(filename.empty()) throw file_access_error{"No database name given"};

    database db;

    const bool showInfo = info != info_level::silent;

    if(showInfo) {
        std::cerr << "Reading database from file '"
                  << filename << "' ... " << std::flush;
    }
    try {
        db.read(filename, what);
        if(showInfo) std::cerr << "done." << std::endl;
    }
    catch(const file_access_error& e) {
        std::cerr << "FAIL" << std::endl;
        throw file_access_error{"Could not read database file '" + filename + "'"};
    }

    return db;
}




} // namespace mc

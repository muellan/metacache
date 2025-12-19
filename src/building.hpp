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

#ifndef MC_BUILDING_HPP_
#define MC_BUILDING_HPP_


#include "database.hpp"
#include "options.hpp"
#include "timer.hpp"


namespace mc {


//-----------------------------------------------------------------------------
/**
 * @brief writes database to disk
 */
void write_database (const database& db, const build_options& opt);



//-----------------------------------------------------------------------------
/**
 * @brief prepares datbase for build, adds targets
 */
void add_to_database (database& db, const build_options& opt, timer& time);


} // namespace mc

#endif

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

#ifndef MC_CLASSIFICATION_HPP_
#define MC_CLASSIFICATION_HPP_


#include "config.hpp"
#include "classification_statistics.hpp"
#include "database.hpp"
#include "timer.hpp"

#include <vector>
#include <string>
#include <iostream>


namespace mc {

// / @brief forward declarations
struct query_options;


//-----------------------------------------------------------------------------
/**
 * @brief Compare taxa by rank in descending order; root > ... > species.
 *        If ranks are equal, compare using sequence ids.
 */
struct rank_higher {
    bool operator() (const taxon* lhs, const taxon* rhs) const noexcept {
        if (lhs->rank() > rhs->rank()) return true;
        if (lhs->rank() < rhs->rank()) return false;
        return lhs->id() < rhs->id();
    }
};

using taxon_count_map = std::map<const taxon*, double, rank_higher>;




//-----------------------------------------------------------------------------
/**
 * @brief classification result target
 */
struct classification_results
{
    explicit
    classification_results (
        std::ostream& readOutputTgt   = std::cout,
        std::ostream& targetOutputTgt = std::cout,
        std::ostream& taxonOutputTgt  = std::cout,
        std::ostream& statusTarget    = std::cerr)
    :
        perReadOut(readOutputTgt),
        perTargetOut(targetOutputTgt),
        perTaxonOut(taxonOutputTgt),
        status(statusTarget)
    {}

    void flush_all_streams() {
        perReadOut.flush();
        perTargetOut.flush();
        perTaxonOut.flush();
        status.flush();
    }

    std::ostream& perReadOut;
    std::ostream& perTargetOut;
    std::ostream& perTaxonOut;
    std::ostream& status;
    timer time;
    classification_statistics statistics;
    taxon_count_map taxCounts; // global taxon -> read count
};




//-----------------------------------------------------------------------------
/**
 * @brief try to map each read from the input files to a taxon
 *        according to the query options
 */
void map_queries_to_targets (
    const std::vector<std::string>& inputFilenames,
    const database&, const query_options&,
    classification_results&);




//-----------------------------------------------------------------------------
/**
 * @brief needed for 'merge' mode: try to map candidates to a taxon
 *        according to the query options
 */
void map_candidates_to_targets (
    std::vector<std::string>&&,
    const std::vector<classification_candidates>&,
    const database&, const query_options&,
    classification_results&);


} // namespace mc

#endif

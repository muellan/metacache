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

#ifndef MC_CLASSIFICATION_H_
#define MC_CLASSIFICATION_H_

#include <vector>
#include <string>
#include <iostream>

#include "classification_statistics.h"
#include "timer.h"
#include "config.h"
#include "query_options.h"
#include "querying.h"


namespace mc {

/// @brief forward declarations
struct query_options;


/*************************************************************************//**
 *
 * @brief classification result target
 *
 *****************************************************************************/
struct classification_results
{
    explicit
    classification_results(std::ostream& readOutputTgt   = std::cout,
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
};


/*************************************************************************//**
 *
 * @brief try to map each read from the input files to a taxon
 *        according to the query options
 *
 *****************************************************************************/
void map_queries_to_targets(
    const std::vector<std::string>& inputFilenames,
    const database&, const query_options&,
    classification_results&);


/*************************************************************************//**
 *
 * @brief needed for 'merge' mode: try to map candidates to a taxon
 *        according to the query options
 *
 *****************************************************************************/
void map_candidates_to_targets(
    const std::vector<std::string>& queryHeaders,
    const std::vector<classification_candidates>&,
    const database&, const query_options&,
    classification_results&);


/*************************************************************************//**
 *
 * @brief Compare taxa by rank in descending order; root > ... > species.
 *        If ranks are equal, compare using sequence ids.
 *
 *****************************************************************************/
struct rank_higher {
    bool operator() (const taxon* lhs, const taxon* rhs) const noexcept {
        if(lhs->rank() > rhs->rank()) return true;
        if(lhs->rank() < rhs->rank()) return false;
        return lhs->id() < rhs->id();
    }
};

using taxon_count_map = std::map<const taxon*, float, rank_higher>;
// using taxon_count_map = std::unordered_map<const taxon*, query_id>;


} // namespace mc

#endif

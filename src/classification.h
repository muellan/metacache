/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2018 André Müller (muellan@uni-mainz.de)
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
    classification_results(std::ostream& mapOutputTgt = std::cout,
                           std::ostream& auxOutputTgt = std::cout,
                           std::ostream& statusTarget = std::cerr)
    :
        mapout(mapOutputTgt),
        auxout(auxOutputTgt),
        status(statusTarget)
    {}

    void flush_all_streams() {
        mapout.flush();
        auxout.flush();
        status.flush();
    }

    std::ostream& mapout;
    std::ostream& auxout;
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
 * @brief sorting taxa by rank
 *
 *****************************************************************************/
struct sortTaxaByRank {
    //sort by rank starting with root, then by id
    bool operator() (const taxon* lhs, const taxon* rhs) {
        if(lhs->rank() > rhs->rank()) return true;
        if(lhs->rank() < rhs->rank()) return false;
        return lhs->id() < rhs->id();
    }
};
using taxon_count_map = std::map<const taxon*, float, sortTaxaByRank>;
// using taxon_count_map = std::unordered_map<const taxon*, query_id>;


} // namespace mc

#endif

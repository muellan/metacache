/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2018 André Müller (muellan@uni-mainz.de)
 *                       & Robin Kobus  (rkobus@uni-mainz.de)
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

#ifndef MC_PRINT_RESULTS_H_
#define MC_PRINT_RESULTS_H_

#include "config.h"
#include "classification_statistics.h"
#include "matches_per_target.h"
#include "query_options.h"


namespace mc {


/*************************************************************************//**
 *
 * @brief prints classification parameters
 *
 *****************************************************************************/
void show_query_parameters(std::ostream&, const query_options&);


/*************************************************************************//**
 *
 * @brief prints taxon information according to output options
 *
 *****************************************************************************/
void show_taxon(std::ostream&,
                         const database&,
                         const classification_output_options&,
                         const taxon* classified);


/*************************************************************************//**
 *
 * @brief prints header for taxon information
 *
 *****************************************************************************/
void show_taxon_header(std::ostream&,
                       const classification_output_options&,
                       const std::string& prefix = "");


/*************************************************************************//**
 *
 * @brief prints top classification candidates
 *
 *****************************************************************************/
void show_matches(std::ostream&,
                  const database&,
                  const classification_candidates&,
                  taxon_rank lowest = taxon_rank::Sequence);


/*************************************************************************//**
 *
 * @brief prints target.window hit matches
 *
 *****************************************************************************/
void show_matches(std::ostream&,
                  const database&,
                  const matches_per_location&,
                  taxon_rank lowest = taxon_rank::Sequence);


/*************************************************************************//**
 *
 * @brief prints target.window hit statistics from database
 *
 *****************************************************************************/
void show_candidate_ranges(std::ostream&,
                           const database&,
                           const classification_candidates&);


/*************************************************************************//**
 *
 * @brief prints a list of all match locations for each classification target
 *
 *****************************************************************************/
void show_matches_per_targets(std::ostream& os,
                              const database&,
                              const matches_per_target&,
                              const classification_output_options&);


/*************************************************************************//**
 *
 * @brief prints the number of matches for each classification target
 *
 *****************************************************************************/
void show_num_matches_per_targets(std::ostream& os,
                                  const database&,
                                  const matches_per_target&,
                                  const classification_output_options&);


/*************************************************************************//**
 *
 * @brief prints the number of matches for each classification target
 *
 *****************************************************************************/
void show_features_of_targets(std::ostream& os,
                                  const database&,
                                  const matches_per_target&,
                                  const classification_output_options&);


/*************************************************************************//**
 *
 * @brief prints per-rank classification statistics
 *
 *****************************************************************************/
void show_taxon_statistics(std::ostream&,
                                    const classification_statistics&,
                                    const std::string& prefix = "");


} // namespace mc


#endif

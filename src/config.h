/*****************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016 André Müller (muellan@uni-mainz.de)
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

/*****************************************************************************
 *
 * @file contains important type aliases
 *
 *****************************************************************************/

#ifndef MC_CONFIG_H_
#define MC_CONFIG_H_


#include "sketch_database.h"
#include "hash_dna.h"


namespace mc {


/********************************************************************
 * @brief limits max. k (= number of letters) per k-mer
 */
#ifdef MC_KMER_TYPE
    using kmer_type = MC_KMER_TYPE ;
#else
    using kmer_type = std::uint32_t;
//    using kmer_type = std::uint64_t;
#endif


/********************************************************************
 * @brief limits number of targets (reference genomes) per database
 */
#ifdef MC_TARGET_ID_TYPE
    using target_id = MC_TARGET_ID_TYPE ;
#else
    using target_id = std::uint16_t;
#endif


/********************************************************************
 * @brief limits max. number of windows per target (reference genome)
 */
#ifdef MC_WINDOW_ID_TYPE
    using window_id = MC_WINDOW_ID_TYPE ;
#else
    using window_id = std::uint32_t;
#endif


/********************************************************************
 * @brief limits max. number of locations per feature
 */
#ifdef MC_LOCATION_LIST_SIZE_TYPE
    using loclist_size_t = MC_LOCATION_LIST_SIZE_TYPE ;
#else
    using loclist_size_t = std::uint8_t;
#endif


/**************************************************************************
 * @brief nucleotide sequence storage type
 */
using sequence = std::string;


/**************************************************************************
 * @brief controls how nucleotide sequences are transformed into 'features'
 */
using sketcher = single_function_min_hasher<kmer_type>;
//using sketcher = minimizer_hasher<kmer_type>;
//using sketcher = kmer_statistics_hasher;
//using sketcher = entropy_hasher<kmer_type>;


/**************************************************************************
 * @brief hash function for database hash multi-map
 */
using feature_hash = std::hash<typename sketcher::feature_type>;


/**************************************************************************
 * @brief stores a multimap that maps features to locations
 *        within targets (= reference genomes)
 */
using database = sketch_database<sequence,sketcher,feature_hash,
                                 target_id,window_id,loclist_size_t>;


/**************************************************************************
 * @brief controls how a classification is derived from a location hit list;
 *        default is a top 2 voting scheme;
 *        will use majority voting scheme if MC_VOTE_TOP > 2
 */
#ifndef MC_VOTE_TOP
    #define MC_VOTE_TOP 2
#endif

using top_matches_in_contiguous_window_range
        = matches_in_contiguous_window_range_top< MC_VOTE_TOP ,target_id>;


} // namespace mc


#endif

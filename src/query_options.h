/*****************************************************************************
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

#ifndef MC_QUERY_OPTIONS_H_
#define MC_QUERY_OPTIONS_H_

#include <thread>

#include "taxonomy.h"


namespace mc {


using taxon_rank = taxonomy::rank;


/*****************************************************************************
 * @brief pairing of queries
 *****************************************************************************/
enum class pairing_mode : unsigned char {
    none, files, sequences
};


/*****************************************************************************
 * @brief alignment mode
 *****************************************************************************/
enum class align_mode : unsigned char {
    none, semi_global
};


/*****************************************************************************
 * @brief how to show mapping
 *****************************************************************************/
enum class map_view_mode : unsigned char {
    none, mapped_only, all
};


/*****************************************************************************
 *
 * @brief how taxon formatting will be done
 *
 *****************************************************************************/
enum class taxon_print_mode : unsigned char {
    rank_id, rank_name, rank_name_id,
    id, name, name_id
};



/*************************************************************************//**
 *
 * @brief database configuration options
 *
 *****************************************************************************/
struct database_query_options
{
    // database tuning options
    bool removeOverpopulatedFeatures = false;
    int maxLocationsPerFeature = -1; //< 0 : use value from database
    float maxLoadFactor = -1;        //< 0 : use database default

    // query sketching options
    int sketchlen = -1;  //< 0 : use value from database
    int winlen    = -1;  //< 0 : use value from database
    int winstride = -1;  //< 0 : use value from database
};



/*************************************************************************//**
 *
 * @brief database configuration options
 *
 *****************************************************************************/
struct query_processing_options
{
    pairing_mode pairing = pairing_mode::none;

    int numThreads = std::thread::hardware_concurrency();

    std::size_t batchSize = 4096 * std::thread::hardware_concurrency();

    //limits number of reads per sequence source (file)
    std::uint_least64_t queryLimit = std::numeric_limits<std::uint_least64_t>::max();
};



/*************************************************************************//**
 *
 * @brief classification options
 *
 *****************************************************************************/
struct classification_options
{
    //ranks/taxa to classify on
    taxon_rank lowestRank  = taxon_rank::Sequence;
    taxon_rank highestRank = taxon_rank::Domain;

    std::uint16_t hitsMin  = 0;  //< 1 : deduced from database parameters
    float hitsDiffFraction = 0.9;
    //maximum range in sequence that read (pair) is expected to be in
    std::size_t insertSizeMax = 0;
};



/*************************************************************************//**
 *
 * @brief ground truth based testing
 *
 *****************************************************************************/
struct evaluation_options
{
    //test precision (ground truth must be available)
    bool precision = false;
    bool taxonCoverage = false;
    //ground truth rank to exclude (for clade exclusion test)
    taxon_rank excludeRank = taxon_rank::none;

    //show known taxon (or complete lineage if 'showLineage' on)
    bool determineGroundTruth = false;
};



/*************************************************************************//**
 *
 * @brief classification output options
 *
 *****************************************************************************/
struct classification_output_options
{
    using taxon_rank = taxonomy::rank;

    //how to show classification (read mappings), if 'none', only summary will be shown
    map_view_mode mapViewMode = map_view_mode::all;

    bool showQueryIds = false;
    //show top candidate sequences and their associated k-mer hash hit count
    bool showTopHits = false;
    //show all k-mer-hash hits in database for each given read
    bool showAllHits = false;
    //show candidate position(s) in reference sequence(s)
    bool showLocations = false;
    //show all ranks that a sequence could be classified on
    bool showLineage = false;
    //what to show of a taxon
    taxon_print_mode showTaxaAs = taxon_print_mode::rank_name;
    bool separateTaxaInfo = false;
    //show ground thruth if available
    bool showGroundTruth = false;
    //make statistics of semi-global alignment scores of queries against
    //target candidate(s)
    bool showAlignment = false;
    //show list of target -> hit mappings
    bool showHitsPerTargetList = false;
    //show error messages?
    bool showErrors = true;

    //ranks/taxa to show
    taxon_rank lowestRank  = taxon_rank::Sequence;
    taxon_rank highestRank = taxon_rank::Domain;

    //show top candidate sequences and their associated k-mer hash hit count
    //prefix for each non-mapping line
    std::string comment = "# ";
    //separates individual mapping fields
    std::string separator = "\t|\t";
};



/*****************************************************************************
 *
 * @brief query options
 *
 *****************************************************************************/
struct query_options
{
    query_processing_options process;
    classification_options classify;
    evaluation_options evaluate;
    classification_output_options output;

    //make a separate output file for each input file
    bool splitFiles = false;
    //output filename
    std::string outfile;
    //show database properties
    bool showDBproperties = false;
    bool showQueryParams = true;
    bool showSummary = true;
};



/*************************************************************************//**
 *
 * @brief forward declarations
 *
 *****************************************************************************/
class args_parser;


/*****************************************************************************
 *
 * @brief command line args -> query options
 *
 *****************************************************************************/
query_options
get_query_options(const args_parser&,
                  const std::vector<std::string>& inputFiles,
                  const query_options& defaults = query_options{});


 /*****************************************************************************
 *
 * @brief command line args -> database query configuration options
 *
 *****************************************************************************/
database_query_options
get_database_query_options(
    const args_parser&,
    const database_query_options& defaults = database_query_options{});


} // namespace mc


#endif

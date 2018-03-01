#ifndef MC_QUERY_OPTIONS_H_
#define MC_QUERY_OPTIONS_H_

#include <thread>

#include "taxonomy.h"
#include "print_info.h"


namespace mc {


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
 * @brief
 *****************************************************************************/
enum class target_view_mode : unsigned char {
    none
};


/*************************************************************************//**
 *
 * @brief database configuration options
 *
 *****************************************************************************/
struct database_query_options
{
    // query sketching options
    int sketchlen = -1;  //< 0 : use value from database
    int winlen    = -1;  //< 0 : use value from database
    int winstride = -1;  //< 0 : use value from database

    // database tuning options
    int maxLocationsPerFeature = -1; //< 0 : use value from database
    bool removeOverpopulatedFeatures = false;
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
    using taxon_rank = taxonomy::rank;

    //ranks/taxa to classify on and to show
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
struct classification_testing_options
{
    //test precision (ground truth must be available)
    bool precision = false;
    bool taxonCoverage = false;
    //ground truth rank to exclude (for clade exclusion test)
    taxon_rank excludedRank = taxon_rank::none;

    //show known taxon (or complete lineage if 'showLineage' on)
    bool showGroundTruth = false;

    //make statistics of semi-global alignment scores of queries against
    //target candidate(s)
    bool alignmentScores = false;
    bool makeAlignments = false;
};



/*************************************************************************//**
 *
 * @brief classification output options
 *
 *****************************************************************************/
struct classification_output_options
{
    //how to show classification (read mappings), if 'none', only summary will be shown
    map_view_mode mapViewMode = map_view_mode::all;

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
    //show alignment if available
    bool showAlignment = false;
    //show list of target -> hit mappings
    bool showHitsPerTargetList = false;
    //show error messages?
    bool showErrors = true;

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
    classification_testing_options test;
    classification_output_options output;

    //make a separate output file for each input file
    bool splitFiles = false;
    //output filename
    std::string outfile;
    //show database properties
    bool showDBproperties = false;
    bool showClassificationParams = true;
    bool showResults = true;
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

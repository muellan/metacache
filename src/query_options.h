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
 *
 * @brief query options
 *
 *****************************************************************************/
struct query_options
{
    using string = std::string;
    using taxon_rank = taxonomy::rank;

    //-----------------------------------------------------
    //output options & formatting
    //-----------------------------------------------------
    //prefix for each non-mapping line
    string comment = "# ";
    //separates individual mapping fields
    string outSeparator = "\t|\t";
    //show database properties
    bool showDBproperties = false;
    //make a separate output file for each input file
    bool splitOutput = false;
    //show top candidate sequences and their associated k-mer hash hit count
    bool showTopHits = false;
    //show all k-mer-hash hits in database for each given read
    bool showAllHits = false;
    //how to show classification (read mappings), if 'none', only summary will be shown
    map_view_mode mapViewMode = map_view_mode::all;
    //show candidate position(s) in reference sequence(s)
    bool showLocations = false;
    //show known taxon (or complete lineage if 'showLineage' on)
    bool showGroundTruth = false;
    //show all ranks that a sequence could be classified on
    bool showLineage = false;
    //what to show of a taxon
    taxon_print_mode showTaxaAs = taxon_print_mode::rank_name;
    bool separateTaxaInfo = false;
    bool showAlignment = false;
    //show error messages?
    bool showErrors = true;

    //-----------------------------------------------------
    //classification options
    //-----------------------------------------------------
    pairing_mode pairing = pairing_mode::none;
    //ranks/taxa to classify on and to show
    taxon_rank lowestRank  = taxon_rank::Sequence;
    taxon_rank highestRank = taxon_rank::Domain;
    //ground truth rank to exclude (for clade exclusion test)
    taxon_rank excludedRank = taxon_rank::none;
    //
    std::uint16_t hitsMin  = 0;  //< 1 : deduced from database parameters
    float hitsDiffFraction = 0.9;
    //maximum range in sequence that read (pair) is expected to be in
    std::size_t insertSizeMax = 0;
    //limit number of reads per file (for quick tests)
    std::uint_least64_t queryLimit =
        std::numeric_limits<std::uint_least64_t>::max();

    //-----------------------------------------------------
    //analysis options
    //-----------------------------------------------------
    //test precision (ground truth must be available)
    bool testPrecision = false;
    bool testCoverage = false;
    bool testAlignment = false;
    //additional file with query -> ground truth mapping
    string sequ2taxonPreFile;

    //-----------------------------------------------------
    //query sampling scheme
    //-----------------------------------------------------
    int sketchlen = -1;  //< 0 : use value from database
    int winlen = -1;     //< 0 : use value from database
    int winstride = -1;  //< 0 : use value from database

    //-----------------------------------------------------
    // tuning parameters
    //-----------------------------------------------------
    int maxLocationsPerFeature = -1; //< 0 : use value from database
    int numThreads = std::thread::hardware_concurrency();
    int batchSize = 0; // automatic
    bool removeOverpopulatedFeatures = false;


    //-----------------------------------------------------
    //filenames
    //-----------------------------------------------------
    string dbfile;
    std::vector<string> infiles;
    string outfile;
};

} // namespace mc


#endif

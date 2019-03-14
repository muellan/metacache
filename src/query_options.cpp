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

#include "args_handling.h"
#include "query_options.h"


namespace mc {

using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::flush;


/*************************************************************************//**
 *
 * @brief command line args -> database query configuration options
 *
 *****************************************************************************/
database_query_options
get_database_query_options(const args_parser& args,
                           const database_query_options& defaults)
{
    database_query_options opt;

    opt.maxLocationsPerFeature = args.get<int>({"max-locations-per-feature",
                                                "max_locations_per_feature"},
                                                defaults.maxLocationsPerFeature);

    opt.removeOverpopulatedFeatures = defaults.removeOverpopulatedFeatures ||
        args.contains({"remove_overpopulated_features" ,
                       "remove-overpopulated-features"});

    //query sketching
    opt.sketchlen = args.get<int>("sketchlen", defaults.sketchlen);
    opt.winlen    = args.get<int>("winlen", defaults.winlen);
    opt.winstride = args.get<int>("winstride",
                                  opt.winlen > 0 ? opt.winlen : defaults.winstride);

    opt.maxLoadFactor = args.get<float>({"max-load-fac", "max_load_fac",
                                         "maxloadfac"},
                                        defaults.maxLoadFactor);

    return opt;
}



/*************************************************************************//**
 *
 * @brief command line args -> query options
 *
 *****************************************************************************/
query_processing_options
get_query_processing_options(const args_parser& args,
                             const vector<string>& infiles,
                             const query_processing_options& defaults)
{
    query_processing_options opt;

    //pairing
    if(args.contains({"paired_files", "paired-files",
                      "pair_files", "pair-files", "pairfiles"}))
    {
        if(infiles.size() > 1) {
            opt.pairing = pairing_mode::files;
        }
    }
    else if(args.contains({"paired_sequences", "paired-sequences",
                           "pair_sequences", "pair-sequences",
                           "pair_sequ", "pair-sequ",
                           "pair_seq", "pair-seq",
                           "pairsequ", "pairseq", "paired"}) )
    {
        opt.pairing = pairing_mode::sequences;
    }

    auto numThreads = args.get<int>("threads", defaults.numThreads);
    if(numThreads >= 1) opt.numThreads = numThreads;

    auto batchSize = args.get<int>("batch-size", defaults.batchSize);
    if(batchSize >= 1) opt.batchSize = batchSize;

    auto numSeq = args.get<int>("per-thread-sequential-queries",
                                defaults.perThreadSequentialQueries);

    if(numSeq >= 1) {
        opt.perThreadSequentialQueries = numSeq;
    }
    else {
        //adapt to number of threads
        if(numThreads < 64) {
            opt.perThreadSequentialQueries = 4;
        } else {
            opt.perThreadSequentialQueries = 8;
        }
    }

    opt.queryLimit = args.get<std::int_least64_t>({"query-limit"},
                                          defaults.queryLimit);

    if(opt.queryLimit < 0) opt.queryLimit = 0;

    return opt;
}



/*************************************************************************//**
 *
 * @brief command line args -> classification options
 *
 *****************************************************************************/
classification_options
get_classification_options(const args_parser& args,
                           const classification_options& defaults)
{
    classification_options opt;

    //classification ranks
    opt.lowestRank  = defaults.lowestRank;
    auto lowestRank = args.get<string>("lowest", "");
    if(!lowestRank.empty()) {
        auto r = taxonomy::rank_from_name(lowestRank);
        if(r < taxon_rank::root) opt.lowestRank = r;
    }

    opt.highestRank = defaults.highestRank;
    auto highestRank = args.get<string>("highest", "");
    if(!highestRank.empty()) {
        auto r = taxonomy::rank_from_name(highestRank);
        if(r < taxon_rank::root) opt.highestRank = r;
    }

    if(opt.lowestRank  > opt.highestRank) {
        opt.lowestRank  = opt.highestRank;
    }
    if(opt.highestRank < opt.lowestRank) {
        opt.highestRank = opt.lowestRank;
    }

    opt.insertSizeMax = args.get<std::size_t>({"insertsize", "insert-size"},
                                                defaults.insertSizeMax);

    opt.hitsMin  = args.get<std::uint16_t>("hitmin", defaults.hitsMin);
    opt.hitsDiffFraction = args.get<float>("hitdiff", defaults.hitsDiffFraction);
    //interprest numbers > 1 as percentage
    if(opt.hitsDiffFraction > 1) opt.hitsDiffFraction *= 0.01;

    opt.maxNumCandidatesPerQuery = args.get<int>({
        "max-cand", "maxcand", "max_cand",
        "max-candidates", "maxcandidates", "max_candidates" },
        defaults.maxNumCandidatesPerQuery);

    if(opt.maxNumCandidatesPerQuery < 1) {
        opt.maxNumCandidatesPerQuery = std::numeric_limits<std::size_t>::max();
    }

    return opt;
}



/*************************************************************************//**
 *
 * @brief command line args -> classification testing options
 *
 *****************************************************************************/
evaluation_options
get_evaluation_options(const args_parser& args,
                       const evaluation_options& defaults)
{
    evaluation_options opt;

    opt.determineGroundTruth = defaults.determineGroundTruth ||
                          args.contains({"ground-truth", "ground_truth",
                                         "groundtruth"});

    //ground truth based analysis options
    opt.taxonCoverage = defaults.taxonCoverage || args.contains("taxon-coverage");

    opt.precision = defaults.precision || opt.taxonCoverage || args.contains("precision");

    opt.excludeRank = defaults.excludeRank;
    auto excludedRank = args.get<string>("exclude", "");
    if(!excludedRank.empty()) {
        auto r = taxonomy::rank_from_name(excludedRank);
        if(r < taxon_rank::root) opt.excludeRank = r;
    }

    return opt;
}



/*************************************************************************//**
 *
 * @brief command line args -> output options
 *
 *****************************************************************************/
classification_output_options
get_classification_output_options(const args_parser& args,
                                  const classification_options& classify,
                                  const evaluation_options& evaluate,
                                  const classification_output_options& defaults)
{
    classification_output_options opt;

    //output tokens
    opt.format.comment = args.get<string>("comment", defaults.format.comment);

    opt.format.column = args.get<string>("separator", defaults.format.column);

    bool separateTaxInfo = args.contains(
        {"separate-cols", "separatecols", "separate_cols",
         "separate-columns", "separatecolumns", "separate_columns"});

    if(separateTaxInfo) {
        opt.collapseUnclassified = false;
        opt.format.taxSeparator = opt.format.column;
        opt.format.rankSuffix   = opt.format.column;
        opt.format.taxidPrefix  = opt.format.column;
        opt.format.taxidSuffix  = "";
    }


    //output formatting
    opt.showQueryIds = args.contains({"queryids", "query-ids", "query_ids",
                                      "queryid", "query-id", "query_id"});

    opt.lowestRank  = classify.lowestRank;
    opt.highestRank = classify.highestRank;

    opt.showLineage = defaults.showLineage || args.contains("lineage");

    opt.showLocations = defaults.showLocations || args.contains("locations");

    opt.showTopHits = defaults.showTopHits || opt.showLocations ||
                      args.contains({"tophits", "top-hits"});

    opt.showAllHits = defaults.showAllHits ||
                      args.contains({"allhits", "all-hits"});

    bool showRanks = !args.contains({"omit-ranks", "omitranks", "omit_ranks"});

    if(args.contains({"taxidsonly","taxids-only","taxids_only",
                      "taxidonly", "taxid-only", "taxid_only"}))
    {
        opt.showTaxaAs = showRanks ? taxon_print_mode::rank_id
                                   : taxon_print_mode::id;
    }
    else if(args.contains({"taxids", "taxid"})) {
        opt.showTaxaAs = showRanks ? taxon_print_mode::rank_name_id
                                   : taxon_print_mode::name_id;
    }
    else {
        opt.showTaxaAs = showRanks ? taxon_print_mode::rank_name
                                   : taxon_print_mode::name;
    }

    opt.mapViewMode = defaults.mapViewMode;
    if(args.contains({"nomap","no-map","noshowmap","nomapping","nomappings"})) {
        opt.mapViewMode = map_view_mode::none;
    }
    else if(args.contains({"mapped-only", "mapped_only", "mappedonly"})) {
        opt.mapViewMode = map_view_mode::mapped_only;
    }
    else if(args.contains({"showallmap","show-all-map","show-all-mappings"})) {
        opt.mapViewMode = map_view_mode::all;
    }

    //showing hits changes the mapping mode!
    if(opt.mapViewMode == map_view_mode::none && opt.showTopHits) {
        opt.mapViewMode = map_view_mode::mapped_only;
    }
    else if(opt.showAllHits) opt.mapViewMode = map_view_mode::all;

    opt.showHitsPerTargetList = defaults.showHitsPerTargetList ||
        args.contains({"hits-per-seq", "hitsperseq", "hits_per_seq",
                       "hits-per-sequence", "hitspersequence", "hits_per_sequence"});

    opt.targetsFile = args.get<string>(
                      {"hits-per-seq", "hitsperseq", "hits_per_seq",
                       "hits-per-sequence", "hitspersequence", "hits_per_sequence"},
                       defaults.targetsFile);

    opt.showQueryIds = opt.showQueryIds || opt.showHitsPerTargetList;

    opt.showTaxAbundances = defaults.showTaxAbundances || args.contains("abundances");

    opt.abundanceFile = args.get<string>("abundances", defaults.abundanceFile);

    opt.showAbundanceEstimatesOnRank = defaults.showAbundanceEstimatesOnRank;
    auto estimationRank = args.get<string>({"abundance-per", "abundances-per",
                                            "abundance_per", "abundances_per"}, "");
    if(!estimationRank.empty()) {
        auto r = taxonomy::rank_from_name(estimationRank);
        if(r < taxon_rank::root) opt.showAbundanceEstimatesOnRank = r;
    }

    if(opt.showTaxAbundances || (opt.showAbundanceEstimatesOnRank != taxon_rank::none))
        opt.makeTaxCounts = true;

    opt.showErrors = defaults.showErrors && !args.contains({"-noerr","-noerrors"});

    opt.showGroundTruth = evaluate.determineGroundTruth;

    opt.showAlignment = defaults.showAlignment ||
        args.contains({"align", "alignment", "showalignment",
                       "showalign", "show-align", "show_align"});

    opt.showDBproperties = defaults.showDBproperties || args.contains("verbose");

    opt.showQueryParams = (defaults.showQueryParams || args.contains("verbose"))
                        && !args.contains({"no-query-params", "noqueryparams",
                                                          "no_query_params"});

    opt.showSummary = (defaults.showSummary || args.contains("verbose")) &&
                      !(args.contains({"no-summary", "nosummary", "no_summary"}));

    //output files
    opt.splitFiles = defaults.splitFiles ||
                     args.contains({"splitout","split-out"});

    opt.queryMappingsFile = args.get<string>("out", defaults.queryMappingsFile);
    if(opt.queryMappingsFile.empty()) {
        //use string after "splitout" as output filename prefix
        opt.queryMappingsFile = args.get<string>({"splitout","split-out"},
                                                 defaults.queryMappingsFile);
    }

    if(opt.targetsFile == opt.queryMappingsFile) opt.targetsFile.clear();

    if(opt.abundanceFile == opt.queryMappingsFile) opt.abundanceFile.clear();

    return opt;
}



/*************************************************************************//**
 *
 * @brief command line args -> query options
 *
 *****************************************************************************/
query_options
get_query_options(const args_parser& args,
                  const vector<string>& infiles,
                  const query_options& defaults)
{
    query_options opt;

    opt.process  = get_query_processing_options(args, infiles, defaults.process);
    opt.classify = get_classification_options(args, defaults.classify);
    opt.evaluate = get_evaluation_options(args, defaults.evaluate);

    opt.output   = get_classification_output_options(args,
                                                     opt.classify,
                                                     opt.evaluate,
                                                     defaults.output);

    return opt;
}


} // namespace mc

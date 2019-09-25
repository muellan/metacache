/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2019 André Müller (muellan@uni-mainz.de)
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

#include <cstdint>
#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <thread>
#include <chrono>

#include "timer.h"
#include "args_handling.h"
#include "args_parser.h"
#include "filesys_utility.h"
#include "cmdline_utility.h"
#include "io_error.h"
#include "io_options.h"
#include "config.h"
#include "sequence_io.h"
#include "taxonomy_io.h"


namespace mc {

using std::string;
using std::cout;
using std::cerr;
using std::flush;
using std::endl;



/*************************************************************************//**
 *
 * @brief database creation parameters
 *
 *****************************************************************************/
struct build_options
{
    int kmerlen = 16;
    int sketchlen = 16;
    int winlen = 128;
    int winstride = 113;

    float maxLoadFactor = -1;           //< 0 : use database default
    int maxLocationsPerFeature = database::max_supported_locations_per_feature();
    bool removeOverpopulatedFeatures = true;

    taxon_rank removeAmbigFeaturesOnRank = taxon_rank::none;
    int maxTaxaPerFeature = 1;

    taxonomy_options taxonomy;
    bool resetParents = false;

    info_level infoLevel = info_level::moderate;

    string dbfile;
    std::vector<string> infiles;
};



/*************************************************************************//**
 *
 * @brief command line args -> database creation parameters
 *
 *****************************************************************************/
build_options
get_build_options(const args_parser& args, const build_options& defaults = {})
{
    build_options opt;

    opt.dbfile = database_name(args);

    opt.infiles = sequence_filenames(args);

    if(args.contains("silent")) {
        opt.infoLevel = info_level::silent;
    } else if(args.contains("verbose")) {
        opt.infoLevel = info_level::verbose;
    }

    opt.kmerlen   = args.get<int>({"kmerlen"}, defaults.kmerlen);
    opt.sketchlen = args.get<int>({"sketchlen"}, defaults.sketchlen);
    opt.winlen    = args.get<int>({"winlen"}, defaults.winlen);
    opt.winstride = args.get<int>({"winstride"}, opt.winlen - opt.kmerlen + 1);

    opt.maxLoadFactor = args.get<float>({"max-load-fac", "max_load_fac",
                                         "maxloadfac"},
                                        defaults.maxLoadFactor);

    opt.maxLocationsPerFeature = args.get<int>({"max-locations-per-feature",
                                                "max_locations_per_feature" },
                                                 defaults.maxLocationsPerFeature);

    opt.removeOverpopulatedFeatures = args.contains({"remove-overpopulated-features",
                                                     "remove_overpopulated_features" });

    opt.removeAmbigFeaturesOnRank = taxonomy::rank_from_name(
        args.get<string>({"remove-ambig-features",
                               "remove_ambig_features"},
                    taxonomy::rank_name(defaults.removeAmbigFeaturesOnRank)));

    opt.maxTaxaPerFeature = args.get<int>({"max-ambig-per-feature",
                                           "max_ambig_per_feature"},
                                            defaults.maxTaxaPerFeature);

    opt.resetParents = args.contains({"reset-parents", "reset_parents"});

    opt.taxonomy = get_taxonomy_options(args);

    return opt;
}



/*************************************************************************//**
 *
 * @brief extract db creation parameters
 *
 *****************************************************************************/
build_options
get_build_options_from_db(const database& db)
{
    build_options opt;

    opt.kmerlen   = db.target_sketcher().kmer_size();
    opt.sketchlen = db.target_sketcher().sketch_size();
    opt.winlen    = db.target_sketcher().window_size();
    opt.winstride = db.target_sketcher().window_stride();

    opt.maxLoadFactor = db.max_load_factor();

    opt.maxLocationsPerFeature = db.max_locations_per_feature();

    return opt;
}



/*************************************************************************//**
 *
 * @brief Alternative way of providing taxonomic mappings.
 *        The input file must be a text file with each line in the format:
 *        accession  accession.version taxid gi
 *        (like the NCBI's *.accession2version files)
 *
 *****************************************************************************/
void rank_targets_with_mapping_file(database& db,
                                    std::set<const taxon*>& targetTaxa,
                                    const string& mappingFile)
{
    if(targetTaxa.empty()) return;

    std::ifstream is {mappingFile};
    if(!is.good()) return;

    const auto fsize = file_size(mappingFile);
    bool showProgress = fsize > 100000000;
    //accession2taxid files have 40-450M lines
    //update progress indicator every 1M lines
    size_t step = 0;
    size_t statStep = 1UL << 20;

    cout << "Try to map sequences to taxa using '" << mappingFile
         << "' (" << std::max(std::streamoff(1),
                              fsize/(1024*1024)) << " MB)" << endl;

    if(showProgress) show_progress_indicator(cout, 0);

    string acc;
    string accver;
    std::uint64_t taxid;
    string gi;

    //skip header
    getline(is, acc);
    acc.clear();

    while(is >> acc >> accver >> taxid >> gi) {
        //target in database?
        //accession.version is the default
        const taxon* tax = db.taxon_with_name(accver);

        if(!tax) {
            tax = db.taxon_with_similar_name(acc);
            if(!tax) tax = db.taxon_with_name(gi);
        }

        //if in database then set parent
        if(tax) {
            auto i = targetTaxa.find(tax);
            if(i != targetTaxa.end()) {
                db.reset_parent(*tax, taxid);
                targetTaxa.erase(i);
                if(targetTaxa.empty()) break;
            }
        }

        if(showProgress && !(++step % statStep)) {
            auto pos = is.tellg();
            show_progress_indicator(cout, pos / float(fsize));
        }
    }

    if(showProgress) clear_current_line(cout);
}



/*************************************************************************//**
 *
 * @return target taxa that have no parent taxon assigned
 *
 *****************************************************************************/
std::set<const taxon*>
unranked_targets(const database& db)
{
    auto res = std::set<const taxon*>{};

    for(const auto& tax : db.target_taxa()) {
        if(!tax.has_parent()) res.insert(&tax);
    }

    return res;
}



/*************************************************************************//**
 *
 * @return all target taxa
 *
 *****************************************************************************/
std::set<const taxon*>
all_targets(const database& db)
{
    auto res = std::set<const taxon*>{};

    for(const auto& tax : db.target_taxa()) {
        res.insert(&tax);
    }

    return res;
}



/*************************************************************************//**
 *
 * @brief try to assign parent taxa to target taxa using mapping files
 *
 *****************************************************************************/
void try_to_rank_unranked_targets(database& db, const build_options& opt)
{
    std::set<const taxon*> unranked;
    if(opt.resetParents)
        unranked = all_targets(db);
    else
        unranked = unranked_targets(db);

    if(!unranked.empty()) {
        if(opt.infoLevel != info_level::silent) {
            cout << unranked.size()
                 << " targets are unranked." << endl;
        }

        for(const auto& file : opt.taxonomy.mappingPostFiles) {
            rank_targets_with_mapping_file(db, unranked, file);
            if(unranked.empty()) break;
        }
    }

    unranked = unranked_targets(db);
    if(opt.infoLevel != info_level::silent) {
        if(unranked.empty()) {
            cout << "All targets are ranked." << endl;
        }
        else {
            cout << unranked.size()
                 << " targets remain unranked." << endl;
        }
    }
}



/*************************************************************************//**
 *
 * @brief look up taxon id based on an identifier (accession number etc.)
 *
 *****************************************************************************/
taxon_id find_taxon_id(
    const std::map<string,taxon_id>& name2tax,
    const string& name)
{
    if(name2tax.empty()) return taxonomy::none_id();
    if(name.empty()) return taxonomy::none_id();

    //try to find exact match
    auto i = name2tax.find(name);
    if(i != name2tax.end()) return i->second;

    //find nearest match
    i = name2tax.upper_bound(name);
    if(i == name2tax.end()) return taxonomy::none_id();

    //if nearest match contains 'name' as prefix -> good enough
    //e.g. accession vs. accession.version
    if(i->first.compare(0,name.size(),name) != 0) return taxonomy::none_id();
    return i->second;
}



/*************************************************************************//**
 *
 * @brief adds reference sequences from *several* files to database
 *
 *****************************************************************************/
void add_targets_to_database(database& db,
    const std::vector<string>& infiles,
    const std::map<string,taxon_id>& sequ2taxid,
    info_level infoLvl = info_level::moderate)
{
    int n = infiles.size();
    int i = 0;

    //read sequences
    for(const auto& filename : infiles) {
        if(infoLvl == info_level::verbose) {
            cout << "  " << filename << " ... " << flush;
        } else if(infoLvl != info_level::silent) {
            show_progress_indicator(cout, i/float(n));
        }

        try {
            auto fileId = extract_accession_string(filename);
            taxon_id fileTaxId = find_taxon_id(sequ2taxid, fileId);

            auto reader = make_sequence_reader(filename);
            while(reader->has_next()) {

                auto sequ = reader->next();
                if(!sequ.data.empty()) {
                    auto seqId  = extract_accession_string(sequ.header);

                    //make sure sequence id is not empty,
                    //use entire header if neccessary
                    if(seqId.empty()) seqId = sequ.header;

                    taxon_id parentTaxId = fileTaxId;

                    if(parentTaxId == taxonomy::none_id())
                        parentTaxId = find_taxon_id(sequ2taxid, seqId);

                    if(parentTaxId == taxonomy::none_id())
                        parentTaxId = extract_taxon_id(sequ.header);

                    if(infoLvl == info_level::verbose) {
                        cout << "[" << seqId;
                        if(parentTaxId > 0) cout << ":" << parentTaxId;
                        cout << "] ";
                    }

                    //try to add to database
                    bool added = db.add_target(
                        sequ.data, seqId, parentTaxId,
                        database::file_source{filename, reader->index()});

                    if(infoLvl == info_level::verbose && !added) {
                        cout << seqId << " not added to database" << endl;
                    }
                }
            }

            if(infoLvl == info_level::verbose) {
                cout << "done." << endl;
            }
        }
        catch(database::target_limit_exceeded_error&) {
            cout << endl;
            cerr << "Reached maximum number of targets per database ("
                 << db.max_target_count() << ").\n"
                 << "See 'README.md' on how to compile MetaCache with "
                 << "support for databases with more reference targets.\n"
                 << endl;
            break;
        }
        catch(std::exception& e) {
            if(infoLvl == info_level::verbose) {
                cout << "FAIL: " << e.what() << endl;
            }
        }

        ++i;
    }

    try {
        //last batch
        db.wait_until_add_target_complete();
    }
    catch(database::target_limit_exceeded_error&) {
        cout << endl;
        cerr << "Reached maximum number of targets per database ("
                << db.max_target_count() << ").\n"
                << "See 'README.md' on how to compile MetaCache with "
                << "support for databases with more reference targets.\n"
                << endl;
    }
    catch(std::exception& e) {
        if(infoLvl == info_level::verbose) {
            cout << "FAIL: " << e.what() << endl;
        }
    }
}



/*************************************************************************//**
 *
 * @brief prepares datbase for build
 *
 *****************************************************************************/
void prepare_database(database& db, const build_options& opt)
{
    if(opt.maxLocationsPerFeature > 0) {
        db.max_locations_per_feature(opt.maxLocationsPerFeature);
    }

    if(opt.maxLoadFactor > 0.4 && opt.maxLoadFactor < 0.99) {
        db.max_load_factor(opt.maxLoadFactor);
        cerr << "Using custom hash table load factor of "
             << opt.maxLoadFactor << endl;
    }

    if(!opt.taxonomy.path.empty()) {
        db.reset_taxa_above_sequence_level(
            make_taxonomic_hierarchy(opt.taxonomy.nodesFile,
                                     opt.taxonomy.namesFile,
                                     opt.taxonomy.mergeFile,
                                     opt.infoLevel) );

        if(opt.infoLevel != info_level::silent) {
            cout << "Taxonomy applied to database." << endl;
        }
    }

    if(db.non_target_taxon_count() < 1 && opt.infoLevel != info_level::silent) {
        cout << "The datbase doesn't contain a taxonomic hierarchy yet.\n"
                  << "You can add one or update later via:\n"
                  << "   ./metacache add <database> -taxonomy <directory>"
                  << endl;
    }

    if(opt.removeAmbigFeaturesOnRank != taxon_rank::none &&
       opt.infoLevel != info_level::silent)
    {
        if(db.non_target_taxon_count() > 1) {
            cout << "Ambiguous features on rank "
                      << taxonomy::rank_name(opt.removeAmbigFeaturesOnRank)
                      << " will be removed afterwards.\n";
        } else {
            cout << "Could not determine amiguous features "
                      << "due to missing taxonomic information.\n";
        }
    }
}



/*************************************************************************//**
 *
 * @brief database features post-processing
 *
 *****************************************************************************/
void post_process_features(database& db, const build_options& opt)
{
    const bool notSilent = opt.infoLevel != info_level::silent;

    if(opt.removeOverpopulatedFeatures) {
        auto old = db.feature_count();
        auto maxlpf = db.max_locations_per_feature() - 1;
        if(maxlpf > 0) { //always keep buckets with size 1
            if(notSilent) {
                cout << "\nRemoving features with more than "
                     << maxlpf << " locations... " << flush;
            }
            auto rem = db.remove_features_with_more_locations_than(maxlpf);

            if(notSilent) {
                cout << rem << " of " << old << " removed." << endl;
                if(rem != old) print_content_properties(db);
            }
        }
    }

    if(opt.removeAmbigFeaturesOnRank != taxon_rank::none &&
        db.non_target_taxon_count() > 1)
    {
        if(notSilent) {
            cout << "\nRemoving ambiguous features on rank "
                 << taxonomy::rank_name(opt.removeAmbigFeaturesOnRank)
                 << "... " << flush;
        }

        auto old = db.feature_count();
        auto rem = db.remove_ambiguous_features(opt.removeAmbigFeaturesOnRank,
                                                opt.maxTaxaPerFeature);

        if(notSilent) {
            cout << rem << " of " << old << "." << endl;
            if(rem != old) print_content_properties(db);
        }

    }
}



/*************************************************************************//**
 *
 * @brief prepares datbase for build, adds targets and writes database to disk
 *
 *****************************************************************************/
void add_to_database(database& db, const build_options& opt)
{
    prepare_database(db, opt);

    const bool notSilent = opt.infoLevel != info_level::silent;
    if(notSilent) print_static_properties(db);

    timer time;
    time.start();

    if(!opt.infiles.empty()) {
        if(notSilent) cout << "Processing reference sequences." << endl;

        const auto initNumTargets = db.target_count();

        auto taxonMap = make_sequence_to_taxon_id_map(
                            opt.taxonomy.mappingPreFilesLocal,
                            opt.taxonomy.mappingPreFilesGlobal,
                            opt.infiles, opt.infoLevel);

        add_targets_to_database(db, opt.infiles, taxonMap, opt.infoLevel);

        if(notSilent) {
            clear_current_line(cout);
            cout << "Added "
                 << (db.target_count() - initNumTargets) << " reference sequences "
                 << "in " << time.seconds() << " s" << endl;

            print_content_properties(db);
        }
    }

    try_to_rank_unranked_targets(db, opt);

    post_process_features(db, opt);

    if(notSilent) {
        cout << "Writing database to file '" << opt.dbfile << "' ... " << flush;
    }
    try {
        db.write(opt.dbfile);
        if(notSilent) cout << "done." << endl;
    }
    catch(const file_access_error&) {
        if(notSilent) cout << "FAIL" << endl;
        cerr << "Could not write database file!" << endl;
    }

    time.stop();

    if(notSilent) {
        cout << "Total build time: " << time.seconds() << " s" << endl;
    }

    //prevents slow deallocation
    db.clear_without_deallocation();
}



/*************************************************************************//**
 *
 * @brief adds reference sequences to an existing database
 *
 *****************************************************************************/
void main_mode_build_modify(const args_parser& args)
{
    auto dbfile = database_name(args);

    if(dbfile.empty()) {
        throw std::invalid_argument{"No database filename provided."};
    }

    cout << "Modify database " << dbfile << endl;

    auto db = make_database<database>(dbfile);

    auto opt = get_build_options(args, get_build_options_from_db(db));

    if(opt.infoLevel != info_level::silent && !opt.infiles.empty()) {
        cout << "Adding reference sequences to database..." << endl;
    }

    add_to_database(db, opt);
}



/*************************************************************************//**
 *
 * @brief builds a database from reference input sequences
 *
 *****************************************************************************/
void main_mode_build(const args_parser& args)
{
    auto opt = get_build_options(args);

    if(opt.infoLevel != info_level::silent) {
        cout << "Building new database from reference sequences." << endl;
    }

    if(opt.dbfile.empty()) {
        throw std::invalid_argument{"No database filename provided."};
    }

    if(opt.infiles.empty()) {
        throw std::invalid_argument{
            "Nothing to do - no reference sequences provided."};
    }

    //configure sketching scheme
    auto sketcher = database::sketcher{};
    sketcher.kmer_size(opt.kmerlen);
    sketcher.sketch_size(opt.sketchlen);
    sketcher.window_size(opt.winlen);
    sketcher.window_stride(opt.winstride);

    auto db = database{sketcher};

    add_to_database(db, opt);
}


} // namespace mc

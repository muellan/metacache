/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2024 André Müller (muellan@uni-mainz.de)
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


#include "building.h"

#include "batch_processing.h"
#include "cmdline_utility.h"
#include "database.h"
#include "filesys_utility.h"
#include "io_error.h"
#include "io_options.h"
#include "options.h"
#include "printing.h"
#include "sequence_io.h"
#include "taxonomy_io.h"
#include "timer.h"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <unordered_set>
#include <utility>
#include <vector>


namespace mc {

using std::string;
using std::cout;
using std::cerr;
using std::flush;
using std::endl;

using taxon_pointer_set = std::unordered_set<const taxon*>;

/*************************************************************************//**
 *
 * @brief Alternative way of providing taxonomic mappings.
 *        The input file must be a text file with each line in the format:
 *        accession  accession.version taxid gi
 *        (like the NCBI's *.accession2version files)
 *
 *****************************************************************************/
void rank_targets_with_mapping_file(taxonomy_cache& taxonomy,
                                    taxon_pointer_set& targetTaxa,
                                    const string& mappingFile,
                                    info_level infoLvl)
{
    const bool showInfo = infoLvl != info_level::silent;

    if (targetTaxa.empty()) return;

    std::ifstream is {mappingFile};
    if (!is.good()) return;

    const auto fsize = file_size(mappingFile);

    if (showInfo) {
        cout << "Try to map sequences to taxa using '" << mappingFile
             << "' (" << std::max(std::streamoff(1),
                                 fsize/(1024*1024)) << " MB)" << endl;
    }

    bool showProgress = showInfo && fsize > 100000000;
    //accession2taxid files have 40-450M lines
    //update progress indicator every 1M lines
    size_t step = 0;
    size_t statStep = 1UL << 20;
    if (showProgress) show_progress_indicator(cout, 0);

    string acc;
    string accver;
    std::uint64_t taxid;
    string gi;

    //skip header
    getline(is, acc);
    acc.clear();

    while (is >> acc >> accver >> taxid >> gi) {
        //target in database?
        //accession.version is the default
        const taxon* tax = taxonomy.taxon_with_name(accver);

        if (!tax) {
            tax = taxonomy.taxon_with_similar_name(acc);
            if (!tax) tax = taxonomy.taxon_with_name(gi);
        }

        //if in database then set parent
        if (tax) {
            auto i = targetTaxa.find(tax);
            if (i != targetTaxa.end()) {
                taxonomy.reset_target_parent(*tax, taxid);
                targetTaxa.erase(i);
                if (targetTaxa.empty()) break;
            }
        }

        if (showProgress && !(++step % statStep)) {
            auto pos = is.tellg();
            show_progress_indicator(cout, pos / float(fsize));
        }
    }

    if (showProgress) clear_current_line(cout);
}



/*************************************************************************//**
 *
 * @return target taxa that have no parent taxon assigned
 *
 *****************************************************************************/
taxon_pointer_set
unranked_targets(const taxonomy_cache& taxonomy)
{
    auto res = taxon_pointer_set{};

    for (const auto& tax : taxonomy.target_taxa()) {
        if (!tax.has_parent()) res.insert(&tax);
    }

    return res;
}



/*************************************************************************//**
 *
 * @return all target taxa
 *
 *****************************************************************************/
taxon_pointer_set
all_targets(const taxonomy_cache& taxonomy)
{
    auto res = taxon_pointer_set{};

    for (const auto& tax : taxonomy.target_taxa()) {
        res.insert(&tax);
    }

    return res;
}



/*************************************************************************//**
 *
 * @brief try to assign parent taxa to target taxa using mapping files
 *
 *****************************************************************************/
void try_to_rank_unranked_targets(taxonomy_cache& taxonomy, const build_options& opt)
{
    taxon_pointer_set unranked;
    if (opt.resetParents)
        unranked = all_targets(taxonomy);
    else
        unranked = unranked_targets(taxonomy);

    if (!unranked.empty()) {
        if (opt.infoLevel != info_level::silent) {
            cout << unranked.size()
                 << " targets are unranked." << endl;
        }

        for (const auto& file : opt.taxonomy.mappingPostFiles) {
            rank_targets_with_mapping_file(taxonomy, unranked, file, opt.infoLevel);
            if (unranked.empty()) break;
        }
    }

    unranked = unranked_targets(taxonomy);
    if (opt.infoLevel != info_level::silent) {
        if (unranked.empty()) {
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
    if (name2tax.empty()) return taxonomy::none_id();
    if (name.empty()) return taxonomy::none_id();

    //try to find exact match
    auto i = name2tax.find(name);
    if (i != name2tax.end()) return i->second;

    //find nearest match
    i = name2tax.upper_bound(name);
    if (i == name2tax.end()) return taxonomy::none_id();

    //if nearest match contains 'name' as prefix -> good enough
    //e.g. accession vs. accession.version
    if (i->first.compare(0,name.size(),name) != 0) return taxonomy::none_id();
    return i->second;
}



// ---------------------------------------------------------------------------
struct input_sequence {
    sequence_reader::header_type header;
    sequence_reader::data_type data;
    database::file_source fileSource;
    taxon_id fileTaxId = 0;
};

using input_batch = std::vector<input_sequence>;



/*************************************************************************//**
 *
 * @brief add batch of reference targets to database
 *
 * @return false, if database not ready / insertion error
 *
 *****************************************************************************/
bool add_targets_to_database(
    database& db,
    int dbPart,
    input_batch& batch,
    const std::map<string,taxon_id>& sequ2taxid,
    info_level infoLvl)
{
    for (auto& seq : batch) {
        if (!db.check_load_factor(dbPart)) return false;
        if (db.add_target_failed(dbPart)) return false;

        if (!seq.data.empty()) {
            auto seqId = extract_accession_string(
                             seq.header, sequence_id_type::any);

            // make sure sequence id is not empty,
            // use entire header if neccessary
            if (seqId.empty()) seqId = seq.header;

            taxon_id parentTaxId = seq.fileTaxId;

            if (parentTaxId == taxonomy::none_id())
                parentTaxId = find_taxon_id(sequ2taxid, seqId);

            if (parentTaxId == taxonomy::none_id())
                parentTaxId = extract_taxon_id(seq.header);

            // try to add to database
            bool added = db.add_target(
                dbPart, seq.data, seqId, parentTaxId, seq.fileSource);

            if (infoLvl == info_level::verbose) {
                cout << "        [" << seqId;
                if (parentTaxId > 0) cout << ":" << parentTaxId;
                cout << "]  " << seq.data.size() << " bp";
                if (added) {
                    cout << endl;
                } else {
                    cout << "  --  not added to database!" << endl;
                }
            }

            seq.data.clear();
        }
    }

    return true;
}



/*************************************************************************//**
 *
 * @brief adds reference sequences from *several* files to database
 *
 *****************************************************************************/
void add_targets_to_database(database& db,
    const std::vector<string>& infiles,
    const std::map<string,taxon_id>& sequ2taxid,
    info_level infoLvl)
{
    // make executor that runs database insertion (concurrently) in batches
    // IMPORTANT: do not use more than one worker thread!
    batch_processing_options<input_sequence> execOpt;
    execOpt.batch_size(8);
    execOpt.queue_size(8);
#ifndef GPU_MODE
    execOpt.concurrency(db.num_parts(), db.num_parts());
#else
    execOpt.concurrency(8, db.num_parts());
#endif

    execOpt.on_error([&] (std::exception& e) {
        if (dynamic_cast<database::target_limit_exceeded_error*>(&e)) {
            cout << endl;
            cerr << "! Reached maximum number of targets per database ("
                 << db.max_target_count() << ").\n"
                 << "! See 'README.md' on how to compile MetaCache with "
                 << "support for databases with more reference targets.\n\n";
        }
        else if (infoLvl == info_level::verbose) {
            cout << "FAIL: " << e.what() << endl;
        }
    });

    execOpt.on_work_done([&] (int id) { db.wait_until_add_target_complete(id); });

    batch_executor<input_sequence> executor { execOpt,
        [&] (int id, auto& batch) {
            return add_targets_to_database(db, id, batch, sequ2taxid, infoLvl);
        }};

    // spawn threads to read sequences
    std::vector<std::future<void>> producers;
    concurrent_progress readingProgress{};
    const size_t numFiles = infiles.size();
    readingProgress.total = numFiles;
    std::mutex outputMtx;

    for (int producerId = 0; producerId < execOpt.num_producers(); ++producerId) {
        producers.emplace_back(std::async(std::launch::async, [&, producerId] {
            auto fileId = readingProgress.counter++;
            while (fileId < numFiles && executor.valid()) {
                const auto& filename = infiles[fileId];
                if (infoLvl == info_level::verbose) {
                    std::lock_guard<std::mutex> lock(outputMtx);
                    cout << "  (" << fileId << '/' << numFiles << ") "
                         << filename << endl;
                }

                try {
                    const auto fileAccession = extract_accession_string(
                                            filename, sequence_id_type::acc_ver);

                    const taxon_id fileTaxId = find_taxon_id(sequ2taxid, fileAccession);

                    if (infoLvl == info_level::verbose) {
                        std::lock_guard<std::mutex> lock(outputMtx);
                        cout << "      accession '" << fileAccession
                             << "' -> taxid " << fileTaxId << endl;
                    }

                    sequence_reader reader{filename};

                    while (reader.has_next() && executor.valid()) {
                        // get (ref to) next input sequence storage and fill it
                        auto& seq = executor.next_item(producerId);
                        seq.fileSource.filename = filename;
                        seq.fileSource.index = reader.index();
                        seq.fileTaxId = fileTaxId;
                        reader.next_header_and_data(seq.header, seq.data);
                    }
                }
                catch(std::exception& e) {
                    if (infoLvl == info_level::verbose) {
                        std::lock_guard<std::mutex> lock(outputMtx);
                        cout << "FAIL: " << e.what() << endl;
                    }
                }

                fileId = readingProgress.counter++;
            }

            executor.finalize_producer(producerId);
        }));
    }

    if (infoLvl == info_level::moderate) {
        show_progress_until_ready(cerr, readingProgress, producers);
    }
    else {
        for (auto& producer : producers) {
            if (producer.valid()) producer.get();
        }
    }

    float progress = readingProgress.progress();
    if (progress < 1.0f) {
        cout << "WARNING: Could only only process " << 100*progress << "% of input files." << endl;
    }

    db.taxo_cache().mark_cached_lineages_outdated();
}



/*************************************************************************//**
 *
 * @brief prepares datbase for build
 *
 *****************************************************************************/
void prepare_database(database& db, const build_options& opt)
{
    const auto dbconf = opt.dbconfig;
    if (dbconf.maxLocationsPerFeature > 0) {
        db.max_locations_per_feature(dbconf.maxLocationsPerFeature);
        cerr << "Max locations per feature set to "
             << int(db.max_locations_per_feature()) << '\n';
    }

    if (dbconf.maxLoadFactor > 0.4 && dbconf.maxLoadFactor < 0.99) {
        db.max_load_factor(dbconf.maxLoadFactor);
        cerr << "Using custom hash table load factor of "
             << db.max_load_factor() << '\n';
    }

    if (!opt.taxonomy.path.empty()) {
        db.taxo_cache().reset_taxa_above_sequence_level(
            make_taxonomic_hierarchy(opt.taxonomy.nodesFile,
                                     opt.taxonomy.namesFile,
                                     opt.taxonomy.mergeFile,
                                     opt.infoLevel) );

        if (opt.infoLevel != info_level::silent) {
            cout << "Taxonomy applied to database." << endl;
        }
    }

    if (db.taxo_cache().non_target_taxon_count() < 1 && opt.infoLevel != info_level::silent) {
        cout << "The datbase doesn't contain a taxonomic hierarchy yet.\n"
             << "You can add one or update later via:\n"
             << "   metacache modify <database> -taxonomy <directory>"
             << endl;
    }

    if (dbconf.removeAmbigFeaturesOnRank != taxon_rank::none &&
       opt.infoLevel != info_level::silent)
    {
        if (db.taxo_cache().non_target_taxon_count() > 1) {
            cout << "Ambiguous features on rank "
                 << taxonomy::rank_name(dbconf.removeAmbigFeaturesOnRank)
                 << " will be removed afterwards.\n";
        } else {
            cout << "Could not determine amiguous features "
                 << "due to missing taxonomic information.\n";
        }
    }

    db.initialize_hash_table(opt.numDbParts);
}



/*************************************************************************//**
 *
 * @brief database features post-processing
 *
 *****************************************************************************/
void post_process_features(database& db, const build_options& opt)
{
    const bool notSilent = opt.infoLevel != info_level::silent;

    const auto& dbconf = opt.dbconfig;

    if (dbconf.removeOverpopulatedFeatures) {
        auto old = db.feature_count();
        auto maxlpf = db.max_locations_per_feature() - 1;
        if (maxlpf > 0) { //always keep buckets with size 1
            if (notSilent) {
                cout << "\nRemoving features with more than "
                     << maxlpf << " locations... " << flush;
            }
            auto rem = db.remove_features_with_more_locations_than(maxlpf);

            if (notSilent) {
                cout << rem << " of " << old << " removed." << endl;
                if (rem != old) print_content_properties(db);
            }
        }
    }

    if (dbconf.removeAmbigFeaturesOnRank != taxon_rank::none &&
        db.taxo_cache().non_target_taxon_count() > 1)
    {
        if (notSilent) {
            cout << "\nRemoving ambiguous features on rank "
                 << taxonomy::rank_name(dbconf.removeAmbigFeaturesOnRank)
                 << "... " << flush;
        }

        auto old = db.feature_count();
        auto rem = db.remove_ambiguous_features(dbconf.removeAmbigFeaturesOnRank,
                                                dbconf.maxTaxaPerFeature);

        if (notSilent) {
            cout << rem << " of " << old << "." << endl;
            if (rem != old) print_content_properties(db);
        }

    }
}



/*************************************************************************//**
 *
 * @brief writes database to disk
 *
 *****************************************************************************/
void write_database(const database& db, const build_options& opt)
{
    const bool notSilent = opt.infoLevel != info_level::silent;

    if (notSilent) {
        cout << "Writing database to file ... " << flush;
    }
    try {
        db.write(opt.dbfile);
        if (notSilent) cout << "done." << endl;
    }
    catch(const file_access_error&) {
        if (notSilent) cout << "FAIL" << endl;
        cerr << "Could not write database file!\n";
    }
}



/*************************************************************************//**
 *
 * @brief prepares datbase for build, adds targets
 *
 *****************************************************************************/
void add_to_database(database& db, const build_options& opt, timer& time)
{
    prepare_database(db, opt);

    const bool notSilent = opt.infoLevel != info_level::silent;
    if (notSilent) print_static_properties(db);

    if (!opt.infiles.empty()) {
        const auto initNumTargets = db.target_count();

        auto taxonMap = make_sequence_to_taxon_id_map(
                            opt.taxonomy.mappingPreFilesLocal,
                            opt.taxonomy.mappingPreFilesGlobal,
                            opt.infiles, opt.infoLevel);

        if (notSilent) cout << "Processing reference sequences." << endl;

        add_targets_to_database(db, opt.infiles, taxonMap, opt.infoLevel);

        if (notSilent) {
            clear_current_line(cout);
            cout << "Added "
                 << (db.target_count() - initNumTargets) << " reference sequences "
                 << "in " << time.seconds() << " s" << endl;

            print_content_properties(db);
        }
    }

    try_to_rank_unranked_targets(db.taxo_cache(), opt);

    post_process_features(db, opt);
}


} // namespace mc

/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2026 André Müller (github.com/muellan)
 *                       & Robin Kobus  (github.com/funatiq)
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


#include "building.hpp"

#include "batch_processing.hpp"
#include "cmdline_utility.hpp"
#include "database.hpp"
#include "filesys_utility.hpp"
#include "io_error.hpp"
#include "io_options.hpp"
#include "options.hpp"
#include "printing.hpp"
#include "sequence_io.hpp"
#include "taxonomy_io.hpp"
#include "timer.hpp"

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


constexpr std::size_t
pow_of_2_floor (std::size_t n) noexcept
{
    // already power of 2 -> exit
    if (n == 0) return 0;
    if ((n & (n - 1)) == 0) return n;
    // no -> find next lower
    --n;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    if (sizeof(std::size_t) > 4) { n |= n >> 32; }
    return (n + 1) >> 1;
}




//-----------------------------------------------------------------------------
/**
 * @brief Alternative way of providing taxonomic mappings.
 *        The input file must be a text file with each line in the format:
 *        accession  accession.version taxid gi
 *        (like the NCBI's *.accession2version files)
 */
void rank_targets_with_mapping_file (
    taxonomy_cache& taxonomy,
    taxon_pointer_set& targetTaxa,
    const string& mappingFile,
    info_level infoLvl)
{
    const bool showInfo = infoLvl != info_level::silent;

    if (targetTaxa.empty()) return;

    std::ifstream is {mappingFile};
    if (not is.good()) return;

    const auto fsize = file_size(mappingFile);

    if (showInfo) {
        cout << "Try to map sequences to taxa using '" << mappingFile
             << "' (" << std::max(std::streamoff(1),
                                 fsize/(1024*1024)) << " MB)" << endl;
    }

    bool showProgress = showInfo && fsize > 100000000;
    // accession2taxid files have 40-450M lines
    // update progress indicator every 1M lines
    size_t step = 0;
    size_t statStep = 1UL << 20;
    if (showProgress) show_progress_indicator(cout, 0);

    string acc;
    string accver;
    std::uint64_t taxid;
    string gi;

    // skip header
    getline(is, acc);
    acc.clear();

    while (is >> acc >> accver >> taxid >> gi) {
        // target in database?
        // accession.version is the default
        const taxon* tax = taxonomy.taxon_with_name(accver);

        if (not tax) {
            tax = taxonomy.taxon_with_similar_name(acc);
            if (not tax) tax = taxonomy.taxon_with_name(gi);
        }

        // if in database then set parent
        if (tax) {
            auto i = targetTaxa.find(tax);
            if (i != targetTaxa.end()) {
                taxonomy.reset_target_parent(*tax, taxid);
                targetTaxa.erase(i);
                if (targetTaxa.empty()) break;
            }
        }

        if (showProgress && not (++step % statStep)) {
            auto pos = is.tellg();
            show_progress_indicator(cout, pos / float(fsize));
        }
    }

    if (showProgress) clear_current_line(cout);
}




//-----------------------------------------------------------------------------
/**
 * @return target taxa that have no parent taxon assigned
 */
taxon_pointer_set
unranked_targets (const taxonomy_cache& taxonomy)
{
    auto res = taxon_pointer_set{};

    for (const auto& tax : taxonomy.target_taxa()) {
        if (not tax.has_parent()) res.insert(&tax);
    }

    return res;
}




//-----------------------------------------------------------------------------
/**
 * @return all target taxa
 */
taxon_pointer_set
all_targets (const taxonomy_cache& taxonomy)
{
    auto res = taxon_pointer_set{};

    for (const auto& tax : taxonomy.target_taxa()) {
        res.insert(&tax);
    }

    return res;
}




//-----------------------------------------------------------------------------
/**
 * @brief try to assign parent taxa to target taxa using mapping files
 */
void try_to_rank_unranked_targets (taxonomy_cache& taxonomy, const build_options& opt)
{
    const bool showInfo = opt.infoLevel != info_level::silent;

    taxon_pointer_set unranked;
    if (opt.resetParents) {
        unranked = all_targets(taxonomy);
    } else {
        unranked = unranked_targets(taxonomy);
    }

    if (not unranked.empty()) {
        if (showInfo) {
            cout << unranked.size()
                 << " targets are unranked (no taxon was assigned)." << endl;
        }

        for (const auto& file : opt.taxonomy.mappingPostFiles) {
            rank_targets_with_mapping_file(taxonomy, unranked, file, opt.infoLevel);
            if (unranked.empty()) break;
        }
    }

    unranked = unranked_targets(taxonomy);
    if (showInfo) {
        if (unranked.empty()) {
            cout << "All targets are ranked (have a taxon assigned)." << endl;
        }
        else {
            cout << unranked.size()
                 << " targets remain unranked (no taxon was assigned)." << endl;
        }
    }
}




//-----------------------------------------------------------------------------
/**
 * @brief look up taxon id based on an identifier (accession number etc.)
 */
taxon_id find_taxon_id (
    const std::map<string,taxon_id>& name2tax,
    const string& name)
{
    if (name2tax.empty()) return taxonomy::none_id();
    if (name.empty()) return taxonomy::none_id();

    // try to find exact match
    auto i = name2tax.find(name);
    if (i != name2tax.end()) return i->second;

    // find nearest match
    i = name2tax.upper_bound(name);
    if (i == name2tax.end()) return taxonomy::none_id();

    // if nearest match contains 'name' as prefix -> good enough
    // e.g. accession vs. accession.version
    if (i->first.compare(0,name.size(),name) != 0) return taxonomy::none_id();
    return i->second;
}



//-----------------------------------------------------------------------------
struct input_sequence 
{
    sequence_reader::header_type header;
    sequence_reader::data_type data;
    database::file_source fileSource;
    taxon_id fileTaxId = 0;
};

using input_batch = std::vector<input_sequence>;




//-----------------------------------------------------------------------------
/**
 * @brief add batch of reference targets to one database part
 *
 * @return false, if database not ready / insertion error
 */
bool add_target_batch_to_database (
    database& db,
    int dbPart,
    input_batch& batch,
    const std::map<string,taxon_id>& sequ2taxid,
    sequence_id_type seqIdType,
    info_level infoLvl)
{
    for (auto& seq : batch) {
        if (not db.check_load_factor(dbPart)) return false;
        if (db.add_target_failed(dbPart)) return false;

        if (not seq.data.empty()) {
            auto seqId = extract_accession_string(seq.header, seqIdType);

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
                cerr << "    P" << dbPart << "  [" << seqId;
                if (parentTaxId > 0) cerr << ":" << parentTaxId;
                cerr << "]  " << seq.data.size() << " bp";
                if (added) {
                    cerr << '\n';
                } else {
                    cerr << "  --  not added to database!\n";
                }
            }

            seq.data.clear();
        }
    }

    return true;
}




//-----------------------------------------------------------------------------
/**
 * @brief adds reference sequences from *several* files to database
 */
void add_targets_to_database (database& db,
    const std::vector<string>& infiles,
    const std::map<string,taxon_id>& sequ2taxid,
    sequence_id_type seqIdType,
    info_level infoLvl)
{
    // make executor that runs database insertion (concurrently) in batches
    // IMPORTANT: do not use more than one consumer thread per database part!
    batch_processing_options<input_sequence> execOpt;

#ifdef GPU_MODE
    const auto numProducers = 8;
#else
    const auto numProducers = std::max(1, std::min(4, int(db.part_count()/2)));
#endif
    // as many consumers as database parts (CPU: virtual, GPU: #GPUs) 
    const auto numConsumers = db.part_count();
    execOpt.batch_size(8);
    execOpt.queue_size(numProducers);
    execOpt.concurrency(numProducers, numConsumers);

    execOpt.on_error([&] (std::exception& e) {
        if (dynamic_cast<database::target_limit_exceeded_error*>(&e)) {
            cout << endl;
            cerr << "! Reached maximum number of targets per database ("
                 << db.max_target_count() << ").\n"
                 << "! See 'README.md' on how to compile MetaCache with "
                 << "support for databases with more reference targets.\n\n";
        }
        else if (infoLvl == info_level::verbose) {
            cerr << "FAIL: " << e.what() << '\n';
        }
    });

    execOpt.on_work_done([&] (int id) { db.wait_until_add_target_complete(id); });

    const size_t numFiles = infiles.size();

    // each consumer inserts sequences into its associated database part 
    batch_executor<input_sequence> executor { execOpt,
        [&] (int id, auto& inSeqBatch) {
            return add_target_batch_to_database(
                        db, id, inSeqBatch, sequ2taxid, seqIdType, infoLvl);
        }};

    // spawn producer threads to read sequences
    std::vector<std::future<void>> producers;
    concurrent_progress readProgress {};
    readProgress.total = numFiles;

    std::mutex outputMtx;

    // each producer reads sequences from input files
    // and places them into the batch executor's buffers (executor.next_item)
    for (int producerId = 0; producerId < execOpt.producer_count(); ++producerId) {
        producers.emplace_back(std::async(std::launch::async, [&, producerId] {
            auto fileId = readProgress.counter++;
            while (fileId < numFiles && executor.valid()) {
                const auto& filename = infiles[fileId];
                if (infoLvl == info_level::verbose) {
                    std::lock_guard<std::mutex> lock(outputMtx);
                    cerr << "  (" << fileId << '/' << numFiles << ") "
                         << filename << endl;
                }

                try {
                    auto fileAccession = extract_accession_string(filename, seqIdType);
                    taxon_id fileTaxId = find_taxon_id(sequ2taxid, fileAccession);

                    if (fileTaxId == taxonomy::none_id() &&
                        seqIdType == sequence_id_type::smart)
                    {
                        fileAccession = extract_accession_string(filename, sequence_id_type::filename);
                        fileTaxId = find_taxon_id(sequ2taxid, fileAccession);
                    }

                    if (infoLvl == info_level::verbose) {
                        std::lock_guard<std::mutex> lock(outputMtx);
                        cerr << "      accession '" << fileAccession
                             << "' -> taxid " << fileTaxId << endl;
                    }

                    sequence_reader reader {filename};

                    while (reader.has_next() && executor.valid()) {
                        // get (ref to) next input sequence storage and fill it
                        auto& seq = executor.next_item(producerId);
                        seq.fileSource.filename = filename;
                        seq.fileSource.index = reader.index();
                        seq.fileTaxId = fileTaxId;
                        reader.next_header_and_data(seq.header, seq.data);
                    }
                }
                catch (std::exception& e) {
                    if (infoLvl == info_level::verbose) {
                        std::lock_guard<std::mutex> lock(outputMtx);
                        cerr << "FAIL: " << e.what() << '\n';
                    }
                }

                fileId = readProgress.counter++;
            }

            executor.finalize_producer(producerId);
        }));
    }

    // wait until operations are finished
    if (infoLvl == info_level::moderate) {
        // show_progress_until_ready(cerr, readProgress, producers);
        show_progress_until_ready(cerr, readProgress, producers);
    }
    else {
        for (auto& producer : producers) {
            if (producer.valid()) producer.get();
        }
    }

    float progress = readProgress.progress();
    if (progress < 1.0f) {
        cout << "WARNING: Could only only process " << 100*progress << "% of input files." << endl;
    }

    db.taxa().mark_cached_lineages_outdated();
}




//-----------------------------------------------------------------------------
/**
 * @brief prepares datbase for build
 */
void prepare_database (database& db, const build_options& opt)
{
    const auto dbconf = opt.dbconfig;
    const bool showInfo = opt.infoLevel != info_level::silent;

    if (dbconf.maxLocationsPerFeature > 0) {
        db.max_locations_per_feature(dbconf.maxLocationsPerFeature);
        if (showInfo) {
            cout << "Max locations per feature set to "
                 << int(db.max_locations_per_feature()) << endl;
        }
    }

    if (dbconf.maxLoadFactor > 0.4 && dbconf.maxLoadFactor < 0.99) {
        db.max_load_factor(dbconf.maxLoadFactor);
        if (showInfo) {
            cout << "Using custom hash table load factor of "
                 << db.max_load_factor() << endl;
        }
    }

    if (not opt.taxonomy.path.empty()) {
        db.taxa().reset_taxa_above_sequence_level(
            make_taxonomic_hierarchy(opt.taxonomy.nodesFile,
                                     opt.taxonomy.namesFile,
                                     opt.taxonomy.mergeFile,
                                     opt.infoLevel) );

        if (showInfo) { cout << "Taxonomy applied to database." << endl; }
    }

    if (db.taxa().non_target_taxon_count() < 1 && showInfo) {
        cout << "The datbase doesn't contain a taxonomic hierarchy yet.\n"
             << "You can add one or update later via:\n"
             << "   metacache modify <database> -taxonomy <directory>"
             << endl;
    }

    if (dbconf.removeAmbigFeaturesOnRank != taxon_rank::none && showInfo) {
        if (db.taxa().non_target_taxon_count() > 1) {
            cout << "Ambiguous features on rank "
                 << taxonomy::rank_name(dbconf.removeAmbigFeaturesOnRank)
                 << " will be removed afterwards." << endl;
        } else {
            cout << "Could not determine amiguous features "
                 << "due to missing taxonomic information." << endl;
        }
    }
}




//-----------------------------------------------------------------------------
/**
 * @brief database features post-processing
 */
void post_process_features (database& db, const build_options& opt)
{
    const bool showInfo = opt.infoLevel != info_level::silent;

    const auto& dbconf = opt.dbconfig;

    if (dbconf.removeOverpopulatedFeatures) {
        auto old = db.feature_count();
        auto maxlpf = db.max_locations_per_feature() - 1;
        if (maxlpf > 0) { // always keep buckets with size 1
            if (showInfo) {
                cout << "\nRemoving features with more than "
                     << maxlpf << " locations... " << flush;
            }
            auto rem = db.remove_features_with_more_locations_than(maxlpf);

            if (showInfo) {
                cout << rem << " of " << old << " removed." << endl;
                if (rem != old) print_content_properties(db);
            }
        }
    }

    if (dbconf.removeAmbigFeaturesOnRank != taxon_rank::none &&
        db.taxa().non_target_taxon_count() > 1)
    {
        if (showInfo) {
            cout << "\nRemoving ambiguous features on rank "
                 << taxonomy::rank_name(dbconf.removeAmbigFeaturesOnRank)
                 << "... " << flush;
        }

        auto old = db.feature_count();
        auto rem = db.remove_ambiguous_features(dbconf.removeAmbigFeaturesOnRank,
                                                dbconf.maxTaxaPerFeature);

        if (showInfo) {
            cout << rem << " of " << old << "." << endl;
            if (rem != old) print_content_properties(db);
        }

    }
}




//-----------------------------------------------------------------------------
/**
 * @brief writes database to disk
 */
void write_database (const database& db, const build_options& opt)
{
    const bool showInfo = opt.infoLevel != info_level::silent;

    if (showInfo) {
        cout << "------------------------------------------------\n"
             << "Writing database to file ... " << endl;
    }
    try {
        db.write(opt.dbfile, opt.infoLevel);
        if (showInfo) cout << "Completed database writing." << endl;
    }
    catch (const file_access_error&) {
        if (showInfo) cout << "FAIL" << endl;
        cerr << "Could not write database file!\n";
    }
}




//-----------------------------------------------------------------------------
/**
 * @brief prepares datbase for build, adds targets
 */
void add_to_database (database& db, const build_options& opt, timer& time)
{
    const bool showInfo = opt.infoLevel != info_level::silent;

    prepare_database(db, opt);

    // number of (virtual) parts to use for initial hash table construction
#ifdef GPU_MODE
    const std::size_t numBuildParts = opt.numDbParts;
#else
    std::size_t numBuildParts = 
        std::min(opt.infiles.size(),
                 std::max(static_cast<std::size_t>(opt.numDbParts),
                          static_cast<std::size_t>(std::max(1u,opt.numThreads/2))));
    // prefer powers of 2 for monolithic databases
    if (opt.numDbParts == 1 && numBuildParts > 2) {
        numBuildParts = pow_of_2_floor(numBuildParts);
    }
#endif

    db.initialize_parts(numBuildParts);

    if (showInfo) print_static_properties(db);

    if (not opt.infiles.empty()) {
        const auto initNumTargets = db.target_count();

        auto taxonMap = make_sequence_to_taxon_id_map(
                            opt.taxonomy.mappingPreFilesLocal,
                            opt.taxonomy.mappingPreFilesGlobal,
                            opt.infiles, opt.infoLevel);

        if (showInfo) {
            cout << "Sequence ID extraction method: "
                 << to_string(opt.sequenceIdType)
                 << "\nProcessing reference sequences." << endl;
        }

        add_targets_to_database(db, opt.infiles, taxonMap,
                                opt.sequenceIdType, opt.infoLevel);

#ifndef GPU_MODE
        // merge virtual DB parts down to actually requested number of parts
        if (db.part_count() > 1 && db.part_count() > opt.numDbParts) {
            if (showInfo) { cout << "Consolidating Hash Tables. " << endl; }
            db.remove_empty_parts();
            db.merge_reduce_max_parts_max_bytes(opt.numDbParts, opt.maxPartBytes,
                                                opt.numThreads, opt.infoLevel);
        }
#endif

        if (showInfo) {
            clear_current_line(cout);
            cout << "Added "
                 << (db.target_count() - initNumTargets) << " reference sequences "
                 << "in " << time.seconds() << " s" << endl;

            print_content_properties(db);
        }
    }

    try_to_rank_unranked_targets(db.taxa(), opt);

    post_process_features(db, opt);
}


} // namespace mc

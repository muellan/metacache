/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2021 André Müller (muellan@uni-mainz.de)
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


#include "database.h"

#include <future>


namespace mc {


// ----------------------------------------------------------------------------
bool database::add_target(part_id dbPart,
                          const sequence& seq, taxon_name sid,
                          taxon_id parentTaxid,
                          file_source source)
{
    // reached hard limit for number of targets
    if(targetCount_.load() >= max_target_count()) {
        throw target_limit_exceeded_error{};
    }
    const std::uint64_t targetCount = targetCount_++;
    // reached hard limit for number of targets
    // in case of multi concurrent increments
    if(targetCount >= max_target_count()) {
        targetCount_--;
        featureStore_.wait_until_add_target_complete(dbPart, targetSketchingOptions_);
        throw target_limit_exceeded_error{};
    }

    if(seq.empty()) return false;

    // don't allow non-unique sequence ids
    if(taxonomyCache_.contains_name(sid)) return false;

    const auto targetId = target_id(targetCount);
    const auto taxid = taxon_id_of_target(targetId);

    // sketch sequence -> insert features
    source.windows = featureStore_.add_target(dbPart, seq, targetId, targetSketchingOptions_);

    if(parentTaxid < 1) parentTaxid = 0;
    const taxon* newtax = nullptr;
    // insert sequence metadata as a new taxon
    newtax = taxonomyCache_.emplace_target_taxon(
        taxid, parentTaxid, sid, std::move(source));

    // allows lookup via sequence id (e.g. NCBI accession number)
    taxonomyCache_.insert_name(std::move(sid), newtax);

    return true;
}



// ----------------------------------------------------------------------------
part_id database::read_meta(const std::string& filename, std::future<void>& taxonomyReaderThread)
{
    std::ifstream is{filename, std::ios::in | std::ios::binary};

    if(!is.good()) {
        throw file_access_error{"Could not read database file '" + filename + "'"};
    }

    //database version info
    using std::uint64_t;
    using std::uint8_t;
    uint64_t dbVer = 0;
    read_binary(is, dbVer);

    if(uint64_t( MC_DB_VERSION ) != dbVer) {
        throw file_read_error{
            "Database " + filename + " (version " + std::to_string(dbVer) + ")"
            + " is incompatible\nwith this version of MetaCache"
            + " (uses version " + std::to_string(MC_DB_VERSION) + ")" };
    }

    //data type widths
    uint8_t featureSize = 0; read_binary(is, featureSize);
    uint8_t targetSize = 0;  read_binary(is, targetSize);
    uint8_t windowSize = 0;  read_binary(is, windowSize);
    uint8_t bucketSize = 0;  read_binary(is, bucketSize);
    uint8_t partSize = 0;    read_binary(is, partSize);
    uint8_t taxidSize = 0;   read_binary(is, taxidSize);
    uint8_t numTaxRanks = 0; read_binary(is, numTaxRanks);

    if( (sizeof(feature) != featureSize) ||
        (sizeof(target_id) != targetSize) ||
        (sizeof(window_id) != windowSize) ||
        (sizeof(bucket_size_type) != bucketSize) ||
        (sizeof(part_id) != partSize))
    {
        throw file_read_error{
            "Database " + filename +
            " is incompatible with this variant of MetaCache" +
            " due to different data type sizes"};
    }

    if( (sizeof(taxon_id) != taxidSize) ||
        (taxonomy::num_ranks != numTaxRanks) )
    {
        throw file_read_error{
            "Database " + filename +
            " is incompatible with this variant of MetaCache" +
            " due to different taxonomy data types"};
    }

    clear();

    //sketching parameters
    read_binary(is, targetSketchingOptions_);
    //TODO: remove this and change db version
    read_binary(is, targetSketchingOptions_);

    //target insertion parameters
    uint64_t maxLocationsPerFeature = 0;
    read_binary(is, maxLocationsPerFeature);
    max_locations_per_feature(maxLocationsPerFeature);

    target_id targetCount = 0;
    read_binary(is, targetCount);
    targetCount_ = targetCount;

    part_id numParts = 0;
    read_binary(is, numParts);

    // read taxon metadata in separate thread
    taxonomyReaderThread = std::async(std::launch::async,
        [&, is = std::move(is)]() mutable {read_binary(is, taxonomyCache_);});

    return numParts;
}


// ----------------------------------------------------------------------------
void database::read_cache(const std::string& filename, part_id partId,
                          concurrent_progress& readingProgress)
{
    std::ifstream is{filename, std::ios::in | std::ios::binary};

    if(!is.good()) {
        throw file_access_error{"Could not read database file '" + filename + "'"};
    }

    //hash table
    read_binary(is, featureStore_, partId, readingProgress);
}


//-------------------------------------------------------------------
void database::read(const std::string& filename, int singlePartId,
                    unsigned replication,
                    scope what)
{
    std::cerr << "Reading database metadata ...\n";

    std::future<void> taxonomyReaderThread;
    part_id numParts = read_meta(filename+".meta", taxonomyReaderThread);

    if(singlePartId >= 0) {
        if(part_id(singlePartId) >= numParts)
            throw std::runtime_error{
                "Database part "+std::to_string(singlePartId)+" is not available. "
                "Database has only "+std::to_string(numParts)+" parts."};

        numParts = 1;
    }

    // read caches in separate threads
    std::vector<std::future<void>> cacheReaderThreads;
    concurrent_progress readingProgress{};

    if(what != scope::metadata_only) {
        featureStore_.prepare_for_query_hash_tables(numParts, replication);

        cacheReaderThreads.reserve(numParts * replication);

        for(unsigned r = 0; r < replication; ++r) {
            if(singlePartId >= 0) {
                cacheReaderThreads.emplace_back(std::async(std::launch::async, [&, r]() {
                    read_cache(filename+".cache"+std::to_string(singlePartId), r, readingProgress);
                }));
            }
            else {
                for(part_id partId = 0; partId < numParts; ++partId) {
                    cacheReaderThreads.emplace_back(std::async(std::launch::async, [&, r, partId]() {
                        read_cache(filename+".cache"+std::to_string(partId), r*numParts+partId, readingProgress);
                    }));
                }
            }
        }
    }

    taxonomyReaderThread.get();
    initialize_taxonomy_caches();

    if(what != scope::metadata_only) {
        std::cerr << "Reading " << numParts << " database part(s) ...\n";

        show_progress_until_ready(std::cerr, readingProgress, cacheReaderThreads);
    }

    std::cerr << "Completed database reading.\n";
}



//-------------------------------------------------------------------
void database::write_meta(const std::string& filename) const
{
    std::cerr << "Writing database metadata to file '" << filename << "' ... ";

    using std::uint64_t;
    using std::uint8_t;

    std::ofstream os{filename, std::ios::out | std::ios::binary};

    if(!os.good()) {
        throw file_access_error{"can't open file " + filename};
    }

    //database version info
    write_binary(os, uint64_t( MC_DB_VERSION ));

    //data type widths
    write_binary(os, uint8_t(sizeof(feature)));
    write_binary(os, uint8_t(sizeof(target_id)));
    write_binary(os, uint8_t(sizeof(window_id)));
    write_binary(os, uint8_t(sizeof(bucket_size_type)));
    write_binary(os, uint8_t(sizeof(part_id)));
    write_binary(os, uint8_t(sizeof(taxon_id)));
    write_binary(os, uint8_t(taxonomy::num_ranks));

    //sketching parameters
    write_binary(os, targetSketchingOptions_);
    //TODO: remove this and change db version
    write_binary(os, targetSketchingOptions_);

    //target insertion parameters
    write_binary(os, uint64_t(max_locations_per_feature()));

    write_binary(os, target_id(targetCount_));

    write_binary(os, num_parts());

    //taxon & target metadata
    write_binary(os, taxonomyCache_);

    std::cerr << "done." << std::endl;
}


//-------------------------------------------------------------------
void database::write_cache(const std::string& filename, part_id partId) const
{
    std::cerr << "Writing database part to file '" << filename << "' ... ";

    using std::uint64_t;
    using std::uint8_t;

    std::ofstream os{filename, std::ios::out | std::ios::binary};

    if(!os.good()) {
        throw file_access_error{"can't open file " + filename};
    }

    //hash table
    write_binary(os, featureStore_, partId);

    std::cerr << "done." << std::endl;
}


//-------------------------------------------------------------------
void database::write(const std::string& filename) const
{
    write_meta(filename+".meta");

    for(part_id partId = 0; partId < num_parts(); ++partId)
        write_cache(filename+".cache"+std::to_string(partId), partId);
}


// ----------------------------------------------------------------------------
void database::clear() {
    taxonomyCache_.clear();
    featureStore_.clear();
}


// ----------------------------------------------------------------------------
/**
 * @brief very dangerous! clears feature map without memory deallocation
 */
void database::clear_without_deallocation() {
    taxonomyCache_.clear();
    featureStore_.clear_without_deallocation();
}




// ----------------------------------------------------------------------------
database
make_database(const std::string& filename, int dbPart, database::scope what, info_level info)
{
    if(filename.empty()) throw file_access_error{"No database name given"};

    database db;

    const bool showInfo = info != info_level::silent;

    if(showInfo) {
        std::cerr << "Reading database from file '"
                  << filename << "' ... " << std::flush;
    }
    try {
        unsigned replication = 1;
        db.read(filename, dbPart, replication, what);
        if(showInfo) std::cerr << "done." << std::endl;
    }
    catch(const file_access_error& e) {
        std::cerr << "FAIL" << std::endl;
        throw file_access_error{"Could not read database file '" + filename + "'"};
    }

    return db;
}


} // namespace mc

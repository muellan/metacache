/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2020 André Müller (muellan@uni-mainz.de)
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


namespace mc {


// ----------------------------------------------------------------------------
bool database::add_target(part_id dbPart,
                          const sequence& seq, taxon_name sid,
                          taxon_id parentTaxid,
                          file_source source)
{
    using std::begin;
    using std::end;

    // reached hard limit for number of targets
    if(targetCount_.load() >= max_target_count()) {
        throw target_limit_exceeded_error{};
    }
    const std::uint64_t targetCount = targetCount_++;
    // reached hard limit for number of targets
    // in case of multi concurrent increments
    if(targetCount >= max_target_count()) {
        throw target_limit_exceeded_error{};
    }

    if(seq.empty()) return false;

    // don't allow non-unique sequence ids
    {
        std::lock_guard<std::mutex> lock(name2taxMtx);
        if(name2tax_.find(sid) != name2tax_.end()) return false;
    }

    const auto targetId = target_id(targetCount);
    const auto taxid = taxon_id_of_target(targetId);

    // sketch sequence -> insert features
    source.windows = featureStore_.add_target(dbPart, seq, targetId, targetSketcher_);

    if(parentTaxid < 1) parentTaxid = 0;
    const taxon* newtax = nullptr;
    // insert sequence metadata as a new taxon
    {
        std::lock_guard<std::mutex> lock(taxaMtx);
        auto nit = taxa_.emplace(
            taxid, parentTaxid, sid,
            taxon_rank::Sequence, std::move(source));

        // should never happen
        if(nit == taxa_.end()) {
            throw std::runtime_error{"target taxon could not be created"};
        }

        newtax = &(*nit);
    }

    // allows lookup via sequence id (e.g. NCBI accession number)
    {
        std::lock_guard<std::mutex> lock(name2taxMtx);
        name2tax_.insert({std::move(sid), newtax});
    }

    return true;
}



// ----------------------------------------------------------------------------
void database::read_single(const std::string& filename, part_id partId, scope what)
{
    std::cerr << "Reading database from file '" << filename << "' ... ";

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

    //data type info
    {
        //data type widths
        uint8_t featureSize = 0; read_binary(is, featureSize);
        uint8_t targetSize = 0;  read_binary(is, targetSize);
        uint8_t windowSize = 0;  read_binary(is, windowSize);
        uint8_t bucketSize = 0;  read_binary(is, bucketSize);
        uint8_t taxidSize = 0;   read_binary(is, taxidSize);
        uint8_t numTaxRanks = 0; read_binary(is, numTaxRanks);

        if( (sizeof(feature) != featureSize) ||
            (sizeof(target_id) != targetSize) ||
            (sizeof(bucket_size_type) != bucketSize) ||
            (sizeof(window_id) != windowSize) )
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
    }

    if(what != scope::hashtable_only) {
        clear();

        //sketching parameters
        read_binary(is, targetSketcher_);
        read_binary(is, querySketcher_);

        //target insertion parameters
        uint64_t maxLocationsPerFeature = 0;
        read_binary(is, maxLocationsPerFeature);
        max_locations_per_feature(maxLocationsPerFeature);

        //taxon metadata
        read_binary(is, taxa_);

        target_id targetCount = 0;
        read_binary(is, targetCount);
        if(targetCount < 1) return;

        targetCount_ = targetCount;

        //update target id -> target taxon lookup table
        targets_.reserve(targetCount);
        for(decltype(targetCount) i = 0 ; i < targetCount; ++i) {
            targets_.push_back(taxa_[taxon_id_of_target(i)]);
        }

        //sequence id lookup
        for(const auto& t : taxa_) {
            if(t.rank() == taxon_rank::Sequence) {
                name2tax_.insert({t.name(), &t});
            }
        }

        mark_cached_lineages_outdated();
        update_cached_lineages(taxon_rank::Sequence);
    }
    else {
        //skip metadata

        //sketching parameters
        sketcher skecherDummy;
        read_binary(is, skecherDummy);
        read_binary(is, skecherDummy);

        //target insertion parameters
        uint64_t dummy;
        read_binary(is, dummy);

        //taxon metadata
        taxonomy dummyTaxa;
        read_binary(is, dummyTaxa);

        target_id targetCount = 0;
        read_binary(is, targetCount);
        if(targetCount < 1) return;
    }

    if(what != scope::metadata_only) {
        //hash table
        read_binary(is, featureStore_, partId);

#ifdef GPU_MODE
        featureStore_.copy_target_lineages_to_gpu(targetLineages_.lineages(), partId);
#endif
    }

    std::cerr << "done." << std::endl;
}


//-------------------------------------------------------------------
void database::read(const std::string& filename, part_id numParts, scope what)
{
#ifndef GPU_MODE
    part_id partId = 0;
    read_single(filename, partId, what);
#else
    if(numParts > num_parts()) {
        numParts = num_parts();
    }

    if(numParts == 1) {
        part_id partId = 0;
        read_single(filename, partId, what);
    }
    else {
        part_id partId = 0;
        read_single(filename+std::to_string(partId), partId, what);

        if(what != scope::metadata_only) {
            for(part_id partId = 1; partId < numParts; ++partId) {
                read_single(filename+std::to_string(partId), partId, scope::hashtable_only);
            }
        }

        featureStore_.enable_all_peer_access(numParts);
    }
#endif
}



//-------------------------------------------------------------------
void database::write_single(const std::string& filename, part_id partId) const
{
    std::cerr << "Writing database part to file '" << filename << "' ... ";

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
    write_binary(os, uint8_t(sizeof(taxon_id)));
    write_binary(os, uint8_t(taxonomy::num_ranks));

    //sketching parameters
    write_binary(os, targetSketcher_);
    write_binary(os, querySketcher_);

    //target insertion parameters
    write_binary(os, uint64_t(max_locations_per_feature()));

    //taxon & target metadata
    write_binary(os, taxa_);
    write_binary(os, target_id(targetCount_));

    //hash table
    write_binary(os, featureStore_, partId);

    std::cerr << "done." << std::endl;
}


//-------------------------------------------------------------------
void database::write(const std::string& filename) const
{
#ifndef GPU_MODE
    part_id partId = 0;
    write_single(filename, partId);
#else
    if(featureStore_.num_parts() == 1) {
        part_id partId = 0;
        write_single(filename, partId);
    }
    else {
        for(part_id partId = 0; partId < featureStore_.num_parts(); ++partId)
            write_single(filename+std::to_string(partId), partId);
    }
#endif
}


// ----------------------------------------------------------------------------
void database::clear() {
    targets_.clear();
    ranksCache_.clear();
    targetLineages_.clear();
    name2tax_.clear();
    featureStore_.clear();
}


// ----------------------------------------------------------------------------
/**
 * @brief very dangerous! clears feature map without memory deallocation
 */
void database::clear_without_deallocation() {
    targets_.clear();
    ranksCache_.clear();
    targetLineages_.clear();
    name2tax_.clear();
    featureStore_.clear_without_deallocation();
}




// ----------------------------------------------------------------------------
database
make_database(const std::string& filename, database::scope what, info_level info)
{
    if(filename.empty()) throw file_access_error{"No database name given"};

    database db;

    const bool showInfo = info != info_level::silent;

    if(showInfo) {
        std::cerr << "Reading database from file '"
                  << filename << "' ... " << std::flush;
    }
    try {
        unsigned numParts = 1;
        db.read(filename, numParts, what);
        if(showInfo) std::cerr << "done." << std::endl;
    }
    catch(const file_access_error& e) {
        std::cerr << "FAIL" << std::endl;
        throw file_access_error{"Could not read database file '" + filename + "'"};
    }

    return db;
}


} // namespace mc

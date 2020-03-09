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
bool database::add_target(const sequence& seq, taxon_name sid,
                          taxon_id parentTaxid,
                          file_source source)
{
    using std::begin;
    using std::end;

    //reached hard limit for number of targets
    if(targets_.size() >= max_target_count()) {
        throw target_limit_exceeded_error{};
    }

    if(seq.empty()) return false;

    //don't allow non-unique sequence ids
    if(name2tax_.find(sid) != name2tax_.end()) return false;

    const auto targetCount = target_id(targets_.size());
    const auto taxid = taxon_id_of_target(targetCount);

    //sketch sequence -> insert features
    source.windows = add_all_window_sketches(seq, targetCount);

    //insert sequence metadata as a new taxon
    if(parentTaxid < 1) parentTaxid = 0;
    auto nit = taxa_.emplace(
        taxid, parentTaxid, sid,
        taxon_rank::Sequence, std::move(source));

    //should never happen
    if(nit == taxa_.end()) {
        throw std::runtime_error{"target taxon could not be created"};
    }

    //allows lookup via sequence id (e.g. NCBI accession number)
    const taxon* newtax = &(*nit);
    name2tax_.insert({std::move(sid), newtax});

    //target id -> taxon lookup table
    targets_.push_back(newtax);

    targetLineages_.mark_outdated();

    return true;
}



// ----------------------------------------------------------------------------
void database::read(const std::string& filename, scope what)

{
    std::ifstream is{filename, std::ios::in | std::ios::binary};

    if(!is.good()) {
        throw file_access_error{"can't open file " + filename};
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

    clear();

    //sketching parameters
    read_binary(is, targetSketcher_);
    read_binary(is, querySketcher_);

    //target insertion parameters
    read_binary(is, maxLocsPerFeature_);

    //taxon metadata
    read_binary(is, taxa_);

    target_id targetCount = 0;
    read_binary(is, targetCount);
    if(targetCount < 1) return;

    //update target id -> target taxon lookup table
    targets_.reserve(targetCount);
    for(decltype(targetCount) i = 0 ; i < targetCount; ++i) {
        targets_.push_back(taxa_[taxon_id_of_target(i)]);
    }

    //sequence id lookup
    name2tax_.clear();
    for(const auto& t : taxa_) {
        if(t.rank() == taxon_rank::Sequence) {
            name2tax_.insert({t.name(), &t});
        }
    }

    mark_cached_lineages_outdated();
    update_cached_lineages(taxon_rank::Sequence);

    if(what == scope::metadata_only) return;

    //hash table
    read_binary(is, features_);
}



// ----------------------------------------------------------------------------
void database::write(const std::string& filename) const
{
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
    write_binary(os, maxLocsPerFeature_);

    //taxon & target metadata
    write_binary(os, taxa_);
    write_binary(os, target_id(targets_.size()));

    //hash table
    write_binary(os, features_);
}



// ----------------------------------------------------------------------------
void database::max_locations_per_feature(bucket_size_type n)
{
    if(n < 1) n = 1;
    if(n >= max_supported_locations_per_feature()) {
        n = max_supported_locations_per_feature();
    }
    else if(n < maxLocsPerFeature_) {
        for(auto i = features_.begin(), e = features_.end(); i != e; ++i) {
            if(i->size() > n) features_.shrink(i, n);
        }
    }
    maxLocsPerFeature_ = n;
}



// ----------------------------------------------------------------------------
database::feature_count_type
database::remove_features_with_more_locations_than(bucket_size_type n)
{
    //note that features are not really removed, because the hashmap
    //does not support erasing keys; instead all values belonging to
    //the key are cleared and the key is kept without values
    feature_count_type rem = 0;
    for(auto i = features_.begin(), e = features_.end(); i != e; ++i) {
        if(i->size() > n) {
            features_.clear(i);
            ++rem;
        }
    }
    return rem;
}



// ----------------------------------------------------------------------------
database::feature_count_type
database::remove_ambiguous_features(taxon_rank r, bucket_size_type maxambig)
{
    feature_count_type rem = 0;

    if(taxa_.empty()) {
        throw std::runtime_error{"no taxonomy available!"};
    }

    if(maxambig == 0) maxambig = 1;

    if(r == taxon_rank::Sequence) {
        for(auto i = features_.begin(), e = features_.end(); i != e; ++i) {
            if(!i->empty()) {
                std::set<target_id> targets;
                for(auto loc : *i) {
                    targets.insert(loc.tgt);
                    if(targets.size() > maxambig) {
                        features_.clear(i);
                        ++rem;
                        break;
                    }
                }
            }
        }
    }
    else {
        for(auto i = features_.begin(), e = features_.end(); i != e; ++i) {
            if(!i->empty()) {
                std::set<const taxon*> taxa;
                for(auto loc : *i) {
                    taxa.insert(targetLineages_[loc.tgt][int(r)]);
                    if(taxa.size() > maxambig) {
                        features_.clear(i);
                        ++rem;
                        break;
                    }
                }
            }
        }
    }
    return rem;
}



// ----------------------------------------------------------------------------
void database::clear() {
    ranksCache_.clear();
    targetLineages_.clear();
    name2tax_.clear();
    features_.clear();
}


// ----------------------------------------------------------------------------
/**
 * @brief very dangerous! clears feature map without memory deallocation
 */
void database::clear_without_deallocation() {
    ranksCache_.clear();
    targetLineages_.clear();
    name2tax_.clear();
    features_.clear_without_deallocation();
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
        db.read(filename, what);
        if(showInfo) std::cerr << "done." << std::endl;
    }
    catch(const file_access_error& e) {
        std::cerr << "FAIL" << std::endl;
        throw file_access_error{"Could not read database file '" + filename + "'"};
    }

    return db;
}


} // namespace mc

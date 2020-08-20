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

#ifndef MC_SKETCH_DATABASE_H_
#define MC_SKETCH_DATABASE_H_


#include "candidate_generation.h"
#include "config.h"
#include "io_error.h"
#include "io_options.h"
#include "taxonomy.h"
#include "typename.h"
#include "version.h"

#ifndef GPU_MODE
    #include "host_hashmap.h"
#else
    #include "gpu_hashmap.cuh"
    #include "query_batch.cuh"
#endif

#include "../dep/hpc_helpers/include/cuda_helpers.cuh"

#include <algorithm>
#include <atomic>
#include <climits>
#include <cstdint>
#include <fstream>
#include <future>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>


namespace mc {

/*************************************************************************//**
 *
 * @brief  maps 'features' (e.g. hash values obtained by min-hashing)
 *         to 'locations' = positions in targets/reference sequences
 *
 *         not copyable, but movable
 *
 * @details terminology
 *  target          reference sequence whose sketches are stored in the DB
 *
 *  taxon           one node in a taxonomic hierarchy (sequence metadata is
 *                  also stored in taxa)
 *
 *  taxon_id        numeric taxon identifier
 *  taxon_name      alphanumeric sequence identifier (e.g. an NCBI accession)
 *
 *  query           sequence (usually short reads) that shall be matched against
 *                  the reference targets
 *
 *  window_id       window index (starting with 0) within a target
 *
 *  location        (window_id, target id) = "window within a target sequence"
 *
 *  full_lineage    path from root -> lowest taxon (vector of const taxon*);
 *                  may have arbitrary length
 *
 *  ranked_lineage  main rank taxon ids (array of const taxon*);
 *                  has fixed length; each index always refers to the same rank
 *
 *  sketcher        function object type, that maps reference sequence
 *                      windows (= sequence interval) to
 *                      sketches (= collection of features of the same type)
 *
 *  feature         single part of a sketch
 *
 *  feature_hash    hash function for (feature -> location) map
 *
 *  target id        internal target (reference sequence) identification

 *  window id        target window identification
 *
 *  loclist_size_t   bucket (location list) size tracking type
 *
 *****************************************************************************/
class database
{
public:
    //---------------------------------------------------------------
    // from global config
    using sequence         = mc::sequence;
    using sketcher         = mc::sketcher;
    using feature_hash     = mc::feature_hash;
    using target_id        = mc::target_id;
    using window_id        = mc::window_id;
    using bucket_size_type = mc::loclist_size_t;
    //-----------------------------------------------------
    // from taxonomy
    using taxon_id       = taxonomy::taxon_id;
    using taxon          = taxonomy::taxon;
    using taxon_rank     = taxonomy::rank;
    using taxon_name     = taxonomy::taxon_name;
    using file_source    = taxonomy::file_source;

    //-----------------------------------------------------
    using match_count_type = std::uint16_t;

    //---------------------------------------------------------------
    enum class scope { everything, metadata_only, hashtable_only };


    //-----------------------------------------------------
    class target_limit_exceeded_error : public std::runtime_error {
    public:
        target_limit_exceeded_error():
            std::runtime_error{"target count limit exceeded"}
        {}
    };


public:
    //-----------------------------------------------------
    /** @brief internal location representation = (window index, target index)
     *         these are stored in the in-memory database and on disk
     */
#ifndef GPU_MODE
    //avoid padding bits
    #pragma pack(push, 1)
#endif
    struct location
    {
        window_id win;
        target_id tgt;

        HOSTDEVICEQUALIFIER
        friend bool
        operator == (const location& a, const location& b) noexcept {
            return (a.tgt == b.tgt) && (a.win == b.win);
        }

        friend bool
        operator < (const location& a, const location& b) noexcept {
            if(a.tgt < b.tgt) return true;
            if(a.tgt > b.tgt) return false;
            return (a.win < b.win);
        }
    };
#ifndef GPU_MODE
    //avoid padding bits
    #pragma pack(pop)
#endif

    //-----------------------------------------------------
    using sketch  = typename sketcher::sketch_type;  //range of features
    using feature = typename sketch::value_type;


private:
    //use negative numbers for sequence level taxon ids
    static constexpr taxon_id
    taxon_id_of_target(target_id id) noexcept {
        return ranked_lineages_of_targets::taxon_id_of_target(id);
    }


    //-----------------------------------------------------
    /// @brief "heart of the database": maps features to target locations
#ifndef GPU_MODE
    using feature_store = host_hashmap<location>;
#else
    using feature_store = gpu_hashmap<feature, location>; //key, value
#endif

public:
    //---------------------------------------------------------------
    using feature_count_type = typename feature_store::feature_count_type;

#ifndef GPU_MODE
    using matches_sorter     = typename feature_store::matches_sorter;
    using match_locations    = typename feature_store::match_locations;
#else
    using match_locations    = std::vector<location>;
#endif


    //---------------------------------------------------------------
    explicit
    database(sketcher targetSketcher = sketcher{}) :
        database{targetSketcher, targetSketcher}
    {}
    //-----------------------------------------------------
    explicit
    database(sketcher targetSketcher, sketcher querySketcher) :
        targetSketcher_{std::move(targetSketcher)},
        querySketcher_{std::move(querySketcher)},
        targetCount_{0},
        featureStore_{},
        taxonomyCache_{}
    {}

    database(const database&) = delete;
    database(database&& other) :
        targetSketcher_{std::move(other.targetSketcher_)},
        querySketcher_{std::move(other.querySketcher_)},
        targetCount_{other.targetCount_.load()},
        featureStore_{std::move(other.featureStore_)},
        taxonomyCache_{std::move(other.taxonomyCache_)}
    {}

    database& operator = (const database&) = delete;
    database& operator = (database&&)      = default;


    //---------------------------------------------------------------
    /**
     * @return const ref to the object that transforms target sequence
     *         snippets into a collection of features
     */
    const sketcher&
    target_sketcher() const noexcept {
        return targetSketcher_;
    }
    //-----------------------------------------------------
    /**
     * @return const ref to the object that transforms query sequence
     *         snippets into a collection of features
     */
    const sketcher&
    query_sketcher() const noexcept {
        return querySketcher_;
    }
    //-----------------------------------------------------
    /**
     * @brief sets the object that transforms query sequence
     *        snippets into a collection of features
     */
    void query_sketcher(const sketcher& s) {
        querySketcher_ = s;
    }
    void query_sketcher(sketcher&& s) {
        querySketcher_ = std::move(s);
    }


    //---------------------------------------------------------------
    void max_locations_per_feature(bucket_size_type n) {
        return featureStore_.max_locations_per_feature(n);
    }
    //-----------------------------------------------------
    bucket_size_type
    max_locations_per_feature() const noexcept {
        return featureStore_.max_locations_per_feature();
    }
    //-----------------------------------------------------
    static bucket_size_type
    max_supported_locations_per_feature() noexcept {
        return feature_store::max_supported_locations_per_feature();
    }

    //-----------------------------------------------------
    feature_count_type
    remove_features_with_more_locations_than(bucket_size_type n) {
        return featureStore_.remove_features_with_more_locations_than(n);
    }


    //---------------------------------------------------------------
    /**
     * @brief  removes features that have more than 'maxambig' different
     *         taxa on a certain taxonomic rank
     *         e.g. remove features that are present in more than 4 phyla
     *
     * @return number of features (hash table buckets) that were removed
     */
    feature_count_type
    remove_ambiguous_features(taxon_rank r, bucket_size_type maxambig) {
        if(taxonomyCache_.empty()) {
            throw std::runtime_error{"no taxonomy available!"};
        }

        return featureStore_.remove_ambiguous_features(r, maxambig, taxonomyCache_);
    }


    //---------------------------------------------------------------
    part_id num_parts() const noexcept {
        return featureStore_.num_parts();
    }


    //---------------------------------------------------------------
    void initialize_hash_table(part_id numParts) {
        featureStore_.initialize_build_hash_tables(numParts);
    }

    //---------------------------------------------------------------
    void initialize_taxonomy_caches() {
        taxonomyCache_.initialize_caches(targetCount_);

#ifdef GPU_MODE
        featureStore_.copy_target_lineages_to_gpus(taxonomyCache_.target_lineages());
#endif
    }


    //---------------------------------------------------------------
    bool add_target(
        part_id dbPart,
        const sequence& seq, taxon_name sid,
        taxon_id parentTaxid = 0,
        file_source source = file_source{});


    //---------------------------------------------------------------
    void wait_until_add_target_complete(part_id gpuId) {
        featureStore_.wait_until_add_target_complete(gpuId, targetSketcher_);
    }


    //---------------------------------------------------------------
    bool add_target_failed() {
        return !featureStore_.valid();
    }


    //---------------------------------------------------------------
    std::uint64_t
    target_count() const noexcept {
        return targetCount_;
    }
    static constexpr std::uint64_t
    max_target_count() noexcept {
        return std::numeric_limits<target_id>::max();
    }
    static constexpr std::uint64_t
    max_windows_per_target() noexcept {
        return std::numeric_limits<window_id>::max();
    }

    //-----------------------------------------------------
    bool empty() const noexcept {
        return featureStore_.empty();
    }


    //---------------------------------------------------------------
    void clear();

    /**
     * @brief very dangerous! clears feature map without memory deallocation
     */
    void clear_without_deallocation();


    //---------------------------------------------------------------
    const taxonomy_cache&
    taxo_cache() const {
        return taxonomyCache_;
    }
    //-----------------------------------------------------
    taxonomy_cache&
    taxo_cache() {
        return taxonomyCache_;
    }


    //---------------------------------------------------------------
#ifndef GPU_MODE
    classification_candidates
    query_host(const sequence& query1, const sequence& query2,
               const candidate_generation_rules& rules,
               matches_sorter& sorter) const
    {
        return featureStore_.query_host(query1, query2, querySketcher_, taxonomyCache_, rules, sorter);
    }
#else
    void
    query_gpu_async(query_batch<location>& queryBatch,
                    part_id hostId,
                    bool copyAllHits,
                    taxon_rank lowestRank) const
    {
        featureStore_.query_async(
            queryBatch,
            hostId,
            query_sketcher(),
            copyAllHits,
            lowestRank);
    }
#endif


    //---------------------------------------------------------------
    void max_load_factor(float lf) {
        featureStore_.max_load_factor(lf);
    }
    //-----------------------------------------------------
    float max_load_factor() const noexcept {
        return featureStore_.max_load_factor();
    }


private:
    /****************************************************************
     * @brief   read database from binary file
     * @details Note that the map is not just de-serialized but
     *          rebuilt by inserting individual keys and values
     *          This should make DB files more robust against changes in the
     *          internal mapping structure.
     ****************************************************************/
    part_id read_meta(const std::string& filename, std::future<void>& taxonomyReaderThread);
    void read_cache(const std::string& filename, part_id partId);

public:
    /****************************************************************
     * @brief   read all database parts from binary files
     ****************************************************************/
    void read(const std::string& filename, int singlePartId,
              unsigned replication,
              scope what = scope::everything);


private:
    /****************************************************************
     * @brief   write database to binary file
     ****************************************************************/
    void write_meta(const std::string& filename) const;
    void write_cache(const std::string& filename, part_id partId) const;

public:
    /****************************************************************
     * @brief   write all database parts to binary files
     ****************************************************************/
    void write(const std::string& filename) const;


    //---------------------------------------------------------------
    std::uint64_t bucket_count() const noexcept {
        return featureStore_.bucket_count();
    }
    //---------------------------------------------------------------
    std::uint64_t feature_count() const noexcept {
        return featureStore_.key_count();
    }
    //---------------------------------------------------------------
    std::uint64_t dead_feature_count() const noexcept {
        return featureStore_.dead_feature_count();
    }
    //---------------------------------------------------------------
    std::uint64_t location_count() const noexcept {
        return featureStore_.value_count();
    }


    //---------------------------------------------------------------
    auto location_list_size_statistics() const {
        return featureStore_.location_list_size_statistics();
    }


    //---------------------------------------------------------------
    void print_feature_map(std::ostream& os) const {
        featureStore_.print_feature_map(os);
    }


    //---------------------------------------------------------------
    void print_feature_counts(std::ostream& os) const {
        featureStore_.print_feature_counts(os);
    }


private:
    //---------------------------------------------------------------
    sketcher targetSketcher_;
    sketcher querySketcher_;
    std::atomic<std::uint64_t> targetCount_;
    mutable feature_store featureStore_;
    taxonomy_cache taxonomyCache_;
};




/*************************************************************************//**
 *
 * @brief pull some types into global namespace
 *
 *****************************************************************************/
using match_locations = database::match_locations;
using location        = database::location;
using feature         = database::feature;




/*************************************************************************//**
 *
 * @brief reads database from file
 *
 *****************************************************************************/
database
make_database(const std::string& filename, int dbPart,
              database::scope = database::scope::everything,
              info_level = info_level::moderate);



} // namespace mc

#endif

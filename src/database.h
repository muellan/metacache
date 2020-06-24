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

#include <climits>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <algorithm>
#include <type_traits>
#include <cstdint>
#include <string>
#include <limits>
#include <memory>
#include <future>
#include <chrono>

#include "version.h"
#include "config.h"
#include "io_error.h"
#include "io_options.h"
#include "stat_combined.cuh"
#include "taxonomy.h"
#include "hash_multimap.h"
#include "gpu_hashmap.cuh"
#include "dna_encoding.h"
#include "typename.h"
#include "query_batch.cuh"

#include "../dep/cudahelpers/cuda_helpers.cuh"

#include "batch_processing.h"

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
    using ranked_lineage = taxonomy::ranked_lineage;
    using full_lineage   = taxonomy::full_lineage;
    using file_source    = taxonomy::file_source;
    using taxon_iterator = taxonomy::const_iterator;
    using taxon_range    = taxonomy::const_range;

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

    const taxon* taxon_of_target(target_id id) const {return targets_[id]; }

public:
    //-----------------------------------------------------
    /** @brief internal location representation = (window index, target index)
     *         these are stored in the in-memory database and on disk
     */
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

    using match_locations = std::vector<location>;


    //-----------------------------------------------------
    using sketch  = typename sketcher::sketch_type;  //range of features
    using feature = typename sketch::value_type;


private:
    //use negative numbers for sequence level taxon ids
    static constexpr taxon_id
    taxon_id_of_target(target_id id) noexcept { return -taxon_id(id)-1; }


    //-----------------------------------------------------
    /// @brief "heart of the database": maps features to target locations
    using feature_store = gpu_hashmap<feature, location>; //key, value

public:
    //---------------------------------------------------------------
    using feature_count_type = typename feature_store::size_type;

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
        targets_{},
        taxa_{},
        ranksCache_{taxa_, taxon_rank::Sequence},
        targetLineages_{taxa_},
        name2tax_{}
    {}

    database(const database&) = delete;
    database(database&& other) :
        targetSketcher_{std::move(other.targetSketcher_)},
        querySketcher_{std::move(other.querySketcher_)},
        targetCount_{other.targetCount_.load()},
        featureStore_{std::move(other.featureStore_)},
        targets_{std::move(other.targets_)},
        taxa_{std::move(other.taxa_)},
        ranksCache_{std::move(other.ranksCache_)},
        targetLineages_{std::move(other.targetLineages_)},
        name2tax_{std::move(other.name2tax_)}
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
        if(taxa_.empty()) {
            throw std::runtime_error{"no taxonomy available!"};
        }

        return featureStore_.remove_ambiguous_features(r, maxambig);
    }


    //---------------------------------------------------------------
    gpu_id num_gpus() const noexcept {
        return featureStore_.num_gpus();
    }


    //---------------------------------------------------------------
    void initialize_hash_table(gpu_id numGPUs) {
        featureStore_.initialize_build_hash_tables(numGPUs);
    }


    //---------------------------------------------------------------
    bool add_target(
        gpu_id dbPart,
        const sequence& seq, taxon_name sid,
        taxon_id parentTaxid = 0,
        file_source source = file_source{});


    //---------------------------------------------------------------
    void wait_until_add_target_complete(gpu_id gpuId) {
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
    static constexpr taxon_id no_taxon_id() noexcept {
        return taxonomy::none_id();
    }


    //---------------------------------------------------------------
    const taxon*
    taxon_with_id(taxon_id id) const noexcept {
        return taxa_[id];
    }
    //-----------------------------------------------------
    /**
     * @brief will only find sequence-level taxon names == sequence id
     */
    const taxon*
    taxon_with_name(const taxon_name& name) const noexcept {
        if(name.empty()) return nullptr;
        auto i = name2tax_.find(name);
        if(i == name2tax_.end()) return nullptr;
        return i->second;
    }
    //-----------------------------------------------------
    /**
     * @brief will find sequence-level taxon names with different versions
     */
    const taxon*
    taxon_with_similar_name(const taxon_name& name) const noexcept {
        if(name.empty()) return nullptr;
        auto i = name2tax_.upper_bound(name);
        if(i == name2tax_.end()) return nullptr;
        const auto s = name.size();
        if(0 != i->first.compare(0,s,name)) return nullptr;
        return i->second;
    }


    //---------------------------------------------------------------
    void reset_taxa_above_sequence_level(taxonomy&& tax) {
        taxa_.erase_above(taxon_rank::Sequence);
        for(auto& t : tax) {
            taxa_.insert_or_replace(std::move(t));
        }
        //re-initialize ranks cache
        mark_cached_lineages_outdated();
        update_cached_lineages(taxon_rank::Sequence);
    }


    //---------------------------------------------------------------
    std::uint64_t
    non_target_taxon_count() const noexcept {
        return taxa_.size() - targetCount_;
    }
    //-----------------------------------------------------
    taxon_range taxa() const {
        return taxa_.full_range();
    }
    taxon_range non_target_taxa() const {
        return taxa_.subrange_from(1);
    }
    //-----------------------------------------------------
    taxon_range target_taxa() const {
        return taxa_.subrange_until(0);
    }


    //---------------------------------------------------------------
    void reset_parent(const taxon& tax, taxon_id parentId) {
        taxa_.reset_parent(tax.id(), parentId);
        mark_cached_lineages_outdated();
    }
    //-----------------------------------------------------
    void reset_parent(const taxon& tax, const taxon& parent) {
        taxa_.reset_parent(tax.id(), parent.id());
        mark_cached_lineages_outdated();
    }

    //---------------------------------------------------------------
    void mark_cached_lineages_outdated() {
        ranksCache_.mark_outdated();
        targetLineages_.mark_outdated();
    }

    //---------------------------------------------------------------
    void update_cached_lineages(taxon_rank rank) const {
        ranksCache_.update(rank);
        targetLineages_.update(target_count());
    }

    //-----------------------------------------------------
    full_lineage
    lineage(const taxon* tax) const noexcept {
        return tax ? lineage(*tax) : full_lineage();
    }
    full_lineage
    lineage(const taxon& tax) const noexcept {
        return taxa_.lineage(tax);
    }
    //-----------------------------------------------------
    const ranked_lineage&
    ranks(const taxon* tax) const noexcept {
        return ranksCache_[tax];
    }
    const ranked_lineage&
    ranks(const taxon& tax) const noexcept {
        return ranksCache_[tax];
    }
    const ranked_lineage&
    ranks(target_id tgt) const noexcept {
        return targetLineages_[tgt];
    }

    //-----------------------------------------------------
    const taxon*
    parent(const taxon* tax) const noexcept {
        return tax ? parent(*tax) : nullptr;
    }
    const taxon*
    parent(const taxon& tax) const noexcept {
        return taxa_.parent(tax);
    }
    //-----------------------------------------------------
    const taxon*
    ancestor(const taxon* tax, taxon_rank r) const noexcept {
        return tax ? ancestor(*tax,r) : nullptr;
    }
    const taxon*
    ancestor(const taxon& tax, taxon_rank r) const noexcept {
        return ranksCache_[tax][int(r)];
    }
    const taxon*
    ancestor(target_id tgt, taxon_rank r) const noexcept {
        return targetLineages_[tgt][int(r)];
    }

    //-----------------------------------------------------
    const taxon*
    next_ranked_ancestor(const taxon* tax) const noexcept {
        return tax ? next_ranked_ancestor(*tax) : nullptr;
    }
    const taxon*
    next_ranked_ancestor(const taxon& tax) const noexcept {
        if(tax.rank() != taxon_rank::none) return &tax;
        for(const taxon* a : ranksCache_[tax]) {
            if(a) return a;
        }
        return nullptr;
    }
    const taxon*
    lowest_ranked_ancestor(target_id tgt, taxon_rank lowest) const noexcept {
        const auto& lineage = targetLineages_[tgt];
        for(int rank = int(lowest); rank < int(taxon_rank::none); ++rank) {
            if(lineage[rank])
                return lineage[rank];
        }
        return nullptr;
    }

    //---------------------------------------------------------------
    const taxon*
    lca(const taxon* ta, const taxon* tb) const {
        return (ta && tb) ? lca(*ta,*tb) : nullptr;
    }
    const taxon*
    lca(const taxon& ta, const taxon& tb) const {
        return taxa_.lca(ta,tb);
    }

    //---------------------------------------------------------------
    const taxon*
    ranked_lca(const taxon* ta, const taxon* tb) const {
        return (ta && tb) ? ranked_lca(*ta,*tb) : nullptr;
    }
    const taxon*
    ranked_lca(const taxon& ta, const taxon& tb) const {
        return ranksCache_.ranked_lca(ta,tb);
    }
    const taxon*
    ranked_lca(target_id ta, target_id tb, taxon_rank lowest) const {
        return targetLineages_.ranked_lca(ta, tb, lowest);
    }


    //---------------------------------------------------------------
    /**
     * @return number of times a taxon is covered by any target in the DB
     */
    std::uint_least64_t
    coverage_of_taxon(const taxon* covered) const {
        return covered ? coverage_of_taxon(*covered) : 0;
    }
    //-----------------------------------------------------
    std::uint_least64_t
    coverage_of_taxon(const taxon& covered) const {
        auto coverage = std::uint_least64_t(0);

        for(const auto& tax : taxa_) {
            if(tax.rank() == taxon_rank::Sequence) {
                for(const taxon* t : taxa_.lineage(tax)) {
                    if(t == &covered) ++coverage;
                }
            }
        }
        return coverage;
    }

    //-----------------------------------------------------
    /**
     * @return true, if taxon is covered by any target in the DB
     */
    bool covers(const taxon* covered) const {
        return covered ? covers(*covered) : false;
    }
    //-----------------------------------------------------
    bool covers(const taxon& covered) const {
        for(const auto& tax : taxa_) {
            if(tax.rank() == taxon_rank::Sequence) {
                for(const taxon* t : taxa_.lineage(tax)) {
                    if(t == &covered) return true;
                }
            }
        }
        return false;
    }


    //---------------------------------------------------------------
    void
    query_gpu_async(query_batch<location>& queryBatch,
                    gpu_id hostId,
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
    void read_single(const std::string& filename, gpu_id gpuId, scope what);

public:
    /****************************************************************
     * @brief   read all database parts from binary files
     ****************************************************************/
    void read(const std::string& filename, gpu_id numGPUs, scope what = scope::everything);


private:
    /****************************************************************
     * @brief   write database to binary file
     ****************************************************************/
    void write_single(const std::string& filename, gpu_id gpuId) const;

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
    /*************************************************************************//**
    *
    * @brief concurrency-safe ranked lineage cache
    *
    *****************************************************************************/
    class ranked_lineages_of_targets
    {
    public:
        //---------------------------------------------------------------
        using ranked_lineage = taxonomy::ranked_lineage;
        using taxon_rank     = taxonomy::rank;

    public:
        //---------------------------------------------------------------
        explicit
        ranked_lineages_of_targets(const taxonomy& taxa)
        :
            taxa_(taxa),
            lins_{},
            outdated_(true)
        {}

        //---------------------------------------------------------------
        ranked_lineages_of_targets(const ranked_lineages_of_targets&) = delete;

        ranked_lineages_of_targets(ranked_lineages_of_targets&& src):
            taxa_(src.taxa_),
            lins_{std::move(src.lins_)},
            outdated_(src.outdated_)
        {}

        ranked_lineages_of_targets& operator = (const ranked_lineages_of_targets&) = delete;
        ranked_lineages_of_targets& operator = (ranked_lineages_of_targets&&) = delete;


        //---------------------------------------------------------------
        void mark_outdated() {
            outdated_ = true;
        }

        //---------------------------------------------------------------
        void update(target_id numTargets = 0) {
            if(!outdated_) return;

            if(numTargets == 0) numTargets = lins_.size();
            lins_.clear();
            for(size_t tgt = 0; tgt < numTargets; ++tgt) {
                lins_.emplace_back(taxa_.ranks(taxon_id_of_target(tgt)));
            }
            outdated_ = false;
        }

        //-----------------------------------------------------
        void clear() {
            lins_.clear();
            outdated_ = true;
        }

        //---------------------------------------------------------------
        const ranked_lineage&
        operator [] (target_id tgt) const {
            assert(outdated_ == false);
            return lins_[tgt];
        }

        //---------------------------------------------------------------
        const taxon*
        ranked_lca(target_id a, target_id b, taxon_rank lowest) const {
            assert(outdated_ == false);
            return taxa_.ranked_lca(lins_[a], lins_[b], lowest);
        }

        //---------------------------------------------------------------
        const std::vector<ranked_lineage>& lineages() const {
            return lins_;
        }

    private:
        //---------------------------------------------------------------
        const taxonomy& taxa_;
        std::vector<ranked_lineage> lins_;
        bool outdated_;
    };


private:
    //---------------------------------------------------------------
    sketcher targetSketcher_;
    sketcher querySketcher_;
    std::atomic<std::uint64_t> targetCount_;
    mutable feature_store featureStore_;
    std::vector<const taxon*> targets_;
    taxonomy taxa_;
    mutable ranked_lineages_cache ranksCache_;
    mutable ranked_lineages_of_targets targetLineages_;
    std::map<taxon_name,const taxon*> name2tax_;

    std::mutex name2taxMtx;
    std::mutex taxaMtx;
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
make_database(const std::string& filename,
              database::scope = database::scope::everything,
              info_level = info_level::moderate);



} // namespace mc

#endif

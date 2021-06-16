/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2021 Robin Kobus  (kobus@uni-mainz.de)
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

#ifndef MC_GPU_HASH_MAP_H_
#define MC_GPU_HASH_MAP_H_


#include "config.h"
#include "query_batch.cuh"
#include "sequence_batch.cuh"
#include "stat_combined.cuh"
#include "taxonomy.h"

#include <atomic>
#include <limits>
#include <memory>


namespace mc {


/*************************************************************************//**
 * TODO
 * @brief   (integer) key -> value hashed multimap
 *          optimized for many values per key (pay attention to max_bucket_size()!
 *          Each bucket contains only one key and all values mapped to that key.
 *          Buckets may also be occupied AND empty ('key only').
 *
 * @details Hash conflicts are resolved with open addressing.
 *
 * @tparam  Key:    key type
 * @tparam  ValueT: value type
 *
 *****************************************************************************/
template<
    class Key,
    class ValueT
>
class gpu_hashmap
{
    class build_hash_table;
    class query_hash_table;

    //-----------------------------------------------------
    /// @brief needed for batched, asynchonous insertion into build hash table
    struct insert_buffer
    {
        insert_buffer() :
            seqBatches_{{
                {MAX_TARGETS_PER_BATCH, MAX_LENGTH_PER_BATCH},
                {MAX_TARGETS_PER_BATCH, MAX_LENGTH_PER_BATCH}
            }},
            currentSeqBatch_(0)
        {};

        sequence_batch<policy::Host>& current_seq_batch() noexcept {
            return seqBatches_[currentSeqBatch_];
        }

        void switch_seq_batch() noexcept {
            currentSeqBatch_ ^= 1;
        }

        std::array<sequence_batch<policy::Host>,2> seqBatches_;
        unsigned currentSeqBatch_;
    };

    //---------------------------------------------------------------
    static constexpr std::size_t
    bucket_size_bits() noexcept {
        return 8;
    }

public:
    //---------------------------------------------------------------
    using key_type           = Key;
    using value_type         = ValueT;
    using bucket_size_type   = loclist_size_t;
    using size_type          = size_t;
    using feature_count_type = size_type;
    using target_id          = mc::target_id;
    using window_id          = mc::window_id;
    using sketcher           = mc::sketcher;

    using ranked_lineage = taxonomy::ranked_lineage;

    //---------------------------------------------------------------
    static constexpr std::size_t
    max_bucket_size() noexcept {
        return (std::size_t(1) << bucket_size_bits()) - 1;
    }

public:
    //---------------------------------------------------------------
    gpu_hashmap();
    ~gpu_hashmap();
    gpu_hashmap(gpu_hashmap&&);
    gpu_hashmap(const gpu_hashmap&) = delete;

private:
    //---------------------------------------------------------------
    /**
     * @brief set number of parts and number of used gpus
     */
    void config_num_gpus(part_id numParts, unsigned replication = 1)
    {
        const size_t numGPUs = size_t(numParts) * replication;

        if(numGPUs > num_gpus())
            throw std::runtime_error{"Number of GPUs must be greater than number of parts"};

        numParts_ = numParts;
        numGPUs_ = part_id(numGPUs);
    }

public:
    //---------------------------------------------------------------
    part_id num_parts() const noexcept { return numParts_; }
    //---------------------------------------------------------------
    part_id num_gpus() const noexcept { return numGPUs_; }

    //---------------------------------------------------------------
    bool add_target_failed(part_id gpuId) const noexcept;
    //---------------------------------------------------------------
    bool check_load_factor(part_id gpuId) const noexcept;

    //---------------------------------------------------------------
    void pop_status();
    void pop_status(part_id gpuId);

    //---------------------------------------------------------------
    size_type bucket_count() const noexcept;
    //-----------------------------------------------------
    size_type key_count() noexcept;
    //-----------------------------------------------------
    size_type dead_feature_count() noexcept { /*TODO count tombstones */ return 0;}
    //-----------------------------------------------------
    size_type value_count() noexcept;
    //-----------------------------------------------------
    bool empty() noexcept {
        return key_count() < 1;
    }

    //---------------------------------------------------------------
    void clear() {};
    void clear_without_deallocation() {
        // not neccssary for gpu version
        clear();
    };

    //---------------------------------------------------------------
    static bucket_size_type
    max_supported_locations_per_feature() noexcept {
        return (max_bucket_size() - 1);
    }
    //-----------------------------------------------------
    void max_locations_per_feature(bucket_size_type n)
    {
        if(n < 1) n = 1;
        if(n >= max_supported_locations_per_feature()) {
            n = max_supported_locations_per_feature();
        }
        maxLocationsPerFeature_ = n;
    }
    //-----------------------------------------------------
    bucket_size_type
    max_locations_per_feature() const noexcept {
        return maxLocationsPerFeature_;
    }

    //-----------------------------------------------------
    feature_count_type
    remove_features_with_more_locations_than(bucket_size_type) {
        std::cerr << "-remove-overpopulated-features is not supported in GPU version\n";
        return 0;
    };


    //---------------------------------------------------------------
    feature_count_type
    remove_ambiguous_features(taxon_rank, bucket_size_type, const taxonomy_cache&) {
        std::cerr << "-remove-ambig-features is not supported in GPU version\n";
        return 0;
    }


    //---------------------------------------------------------------
    static constexpr float default_max_load_factor() noexcept {
        return 0.95;
    }
    //-----------------------------------------------------
    void max_load_factor(float lf) {
        if(lf > 1.0f) lf = 1.0f;
        if(lf < 0.1f) lf = 0.1f;

        using std::abs;
        if(abs(maxLoadFactor_ - lf) > 0.00001f) {
            maxLoadFactor_ = lf;
        }
    }
    //-----------------------------------------------------
    float max_load_factor() const noexcept {
        return maxLoadFactor_;
    }

    //---------------------------------------------------------------
    statistics_accumulator_gpu<policy::Host>
    location_list_size_statistics();

    //---------------------------------------------------------------
    void print_feature_map(std::ostream&) const {
        std::cerr << "feature map is not available in GPU version\n";
    }

    //---------------------------------------------------------------
    void print_feature_counts(std::ostream&) const {
        std::cerr << "featurecounts are not available in GPU version\n";
    }

    /****************************************************************
     * @brief allocate gpu hash tables for database building
     */
    void initialize_build_hash_tables(part_id numGPUs);

    /****************************************************************
     * @brief set number of db parts and total number of gpus;
     *        resize vector of gpu hash tables before database loading
     */
    void prepare_for_query_hash_tables(part_id numParts, unsigned replication);

    /****************************************************************
     * @brief split sequence into batches and insert into build hash table
     */
    window_id add_target(
        part_id gpuId,
        const sequence& seq,
        target_id tgt,
        const sketcher& targetSketcher);
    //-----------------------------------------------------
    window_id add_target(
        part_id gpuId,
        sequence::const_iterator first,
        sequence::const_iterator last,
        target_id tgt,
        const sketcher& targetSketcher);

    //---------------------------------------------------------------
    void insert(
        part_id gpuId, sequence_batch<policy::Host>& seqBatchHost,
        const sketcher& targetSketcher);

    //-----------------------------------------------------
    void wait_until_add_target_complete(part_id gpuId, const sketcher& targetSketcher);

private:
    //---------------------------------------------------------------
    template<class Hashtable>
    void query_hashtables_async(
        const std::vector<Hashtable>& hashtables,
        query_batch<value_type>& batch,
        part_id hostId,
        const sketcher& querySketcher,
        bool copyAllHits,
        taxon_rank lowestRank) const;

public:
    //-----------------------------------------------------
    void query_async(
        query_batch<value_type>& batch,
        part_id hostId,
        const sketcher& querySketcher,
        bool copyAllHits,
        taxon_rank lowestRank) const;


    /****************************************************************
     * @brief deserialize hashmap from input stream
     */
     friend void read_binary(std::istream& is, gpu_hashmap& m, part_id gpuId,
                             concurrent_progress& readingProgress)
    {
        m.deserialize(is, gpuId, readingProgress);
    }

    /****************************************************************
     * @brief serialize hashmap to output stream
     */
    friend void write_binary(std::ostream& os, gpu_hashmap& m, part_id gpuId) {
        m.serialize(os, gpuId);
    }

    //---------------------------------------------------------------
    void copy_target_lineages_to_gpus(const std::vector<ranked_lineage>& lins);

    //---------------------------------------------------------------
    /**
     * @brief enable access between consecutive gpus of each replicated db
     */
    void enable_peer_access();


private:
    //---------------------------------------------------------------
    void deserialize(std::istream& is, part_id gpuId, concurrent_progress& readingProgress);

    //---------------------------------------------------------------
    /**
     * @brief binary serialization of all non-emtpy buckets
     */
    void serialize(std::ostream& os, part_id gpuId);


private:
    //---------------------------------------------------------------
    part_id numGPUs_;
    part_id numParts_;
    float maxLoadFactor_;
    std::uint64_t maxLocationsPerFeature_;

    std::vector<insert_buffer> insertBuffers_;

    /**
     * @brief multi value hashtables for building,
     *        which maps feature -> values
     */
    std::vector<std::unique_ptr<build_hash_table>> buildHashTables_;
    /**
     * @brief single value hashtable for querying,
     *        which maps feature -> values pointer
     */
    std::vector<std::unique_ptr<query_hash_table>> queryHashTables_;

    std::vector<ranked_lineage *> lineages_;
};


} // namespace mc


#endif

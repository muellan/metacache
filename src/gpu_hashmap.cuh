#ifndef MC_GPU_HASH_MAP_H_
#define MC_GPU_HASH_MAP_H_

#include <limits>
#include <memory>
#include <atomic>

#include "config.h"
#include "sequence_batch.cuh"
#include "query_batch.cuh"
#include "taxonomy.h"
#include "stat_combined.cuh"

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

public:
    //---------------------------------------------------------------
    unsigned num_parts() const noexcept {
        return numGPUs_;
    }

    //---------------------------------------------------------------
    bool valid() const noexcept { return valid_; };

    //---------------------------------------------------------------
    void pop_status();
    void pop_status(gpu_id gpuId);

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
    remove_ambiguous_features(taxon_rank, bucket_size_type, const ranked_lineages_of_targets&) {
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
     * @brief allocate gpu hash table for database building
     */
    gpu_id initialize_build_hash_tables(gpu_id numGPUs);

    /****************************************************************
     * @brief split sequence into batches and insert into build hash table
     */
    window_id add_target(
        gpu_id gpuId,
        const sequence& seq,
        target_id tgt,
        const sketcher& targetSketcher);
    //-----------------------------------------------------
    window_id add_target(
        gpu_id gpuId,
        sequence::const_iterator first,
        sequence::const_iterator last,
        target_id tgt,
        const sketcher& targetSketcher);

    //---------------------------------------------------------------
    void insert(
        gpu_id gpuId, sequence_batch<policy::Host>& seqBatchHost,
        const sketcher& targetSketcher);

    //-----------------------------------------------------
    void wait_until_add_target_complete(gpu_id gpuId, const sketcher& targetSketcher);

    //---------------------------------------------------------------
    void query_async(
        query_batch<value_type>& batch,
        gpu_id hostId,
        const sketcher& querySketcher,
        bool copyAllHits,
        taxon_rank lowestRank) const;


    /****************************************************************
     * @brief deserialize hashmap from input stream
     */
     friend void read_binary(std::istream& is, gpu_hashmap& m, gpu_id gpuId) {
        m.deserialize(is, gpuId);
    }

    /****************************************************************
     * @brief serialize hashmap to output stream
     */
    friend void write_binary(std::ostream& os, gpu_hashmap& m, gpu_id gpuId) {
        m.serialize(os, gpuId);
    }

    //---------------------------------------------------------------
    void copy_target_lineages_to_gpu(const std::vector<ranked_lineage>& lins, gpu_id gpuId);

    //---------------------------------------------------------------
    void enable_all_peer_access(gpu_id numGPUs);


private:
    //---------------------------------------------------------------
    void deserialize(std::istream& is, gpu_id gpuId);

    //---------------------------------------------------------------
    /**
     * @brief binary serialization of all non-emtpy buckets
     */
    void serialize(std::ostream& os, gpu_id gpuId);


private:
    //---------------------------------------------------------------
    gpu_id numGPUs_;
    float maxLoadFactor_;
    std::uint64_t maxLocationsPerFeature_;
    std::atomic_bool valid_;

    std::vector<insert_buffer> insertBuffers_;

    /**
     * @brief multi value hashtables for building,
     *        which maps feature -> values
     */
    std::vector<build_hash_table> buildHashTables_;
    /**
     * @brief single value hashtable for querying,
     *        which maps feature -> values pointer
     */
    std::vector<query_hash_table> queryHashTables_;
};


} // namespace mc


#endif

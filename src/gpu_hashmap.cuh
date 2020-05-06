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

public:
    //---------------------------------------------------------------
    using key_type         = Key;
    using value_type       = ValueT;
    using bucket_size_type = loclist_size_t;
    using size_type        = size_t;
    using target_id        = mc::target_id;
    using window_id        = mc::window_id;
    using sketcher         = mc::sketcher;

    using ranked_lineage = taxonomy::ranked_lineage;

public:
    //---------------------------------------------------------------
    gpu_hashmap();
    ~gpu_hashmap();
    gpu_hashmap(gpu_hashmap&&);
    gpu_hashmap(const gpu_hashmap&) = delete;

public:
    //---------------------------------------------------------------
    int num_gpus() const noexcept {
        return numGPUs_;
    }

    //---------------------------------------------------------------
    bool valid() const noexcept;

    //---------------------------------------------------------------
    void pop_status();
    void pop_status(int gpuId);

    //---------------------------------------------------------------
    static constexpr std::size_t
    max_bucket_size() noexcept {
        return std::numeric_limits<bucket_size_type>::max();
    }

    //---------------------------------------------------------------
    size_type bucket_count() const noexcept;
    //-----------------------------------------------------
    size_type key_count() noexcept;
    //-----------------------------------------------------
    size_type tombstone_count() noexcept { /*TODO*/ return 0;}
    //-----------------------------------------------------
    size_type value_count() noexcept;
    //-----------------------------------------------------
    bool empty() noexcept {
        return key_count() < 1;
    }

    //---------------------------------------------------------------
    static constexpr float default_max_load_factor() noexcept {
        return 0.95;
    }
    //-----------------------------------------------------
    float max_load_factor() const noexcept {
        return maxLoadFactor_;
    }

    //---------------------------------------------------------------
    statistics_accumulator_gpu<policy::Host>
    location_list_size_statistics();

    /****************************************************************
     * @brief allocate gpu hash table for database building
     */
    int initialize_build_hash_table(int numGPUs, std::uint64_t maxLocsPerFeature);

    /****************************************************************
     * @brief split sequence into batches and insert into build hash table
     */
    window_id add_target(
        int gpuId,
        const sequence& seq,
        target_id tgt,
        const sketcher& targetSketcher);
    //-----------------------------------------------------
    window_id add_target(
        int gpuId,
        sequence::const_iterator first,
        sequence::const_iterator last,
        target_id tgt,
        const sketcher& targetSketcher);

    //---------------------------------------------------------------
    void insert(
        int gpuId, sequence_batch<policy::Host>& seqBatchHost,
        const sketcher& targetSketcher);

    //-----------------------------------------------------
    void wait_until_add_target_complete(int gpuId, const sketcher& targetSketcher);
    void wait_until_add_target_complete(const sketcher& targetSketcher);

    //---------------------------------------------------------------
    void query_async(
        query_batch<value_type>& batch,
        const sketcher& querySketcher,
        bucket_size_type maxLocationPerFeature,
        bool copyAllHits,
        taxon_rank lowestRank) const;


    /****************************************************************
     * @brief deserialize hashmap from input stream
     */
     friend void read_binary(std::istream& is, gpu_hashmap& m) {
        m.deserialize(is);
    }

    /****************************************************************
     * @brief serialize hashmap to output stream
     */
    friend void write_binary(std::ostream& os, gpu_hashmap& m, int gpuId) {
        m.serialize(os, gpuId);
    }

    //---------------------------------------------------------------
    void copy_target_lineages_to_gpu(const std::vector<ranked_lineage>& lins);


private:
    //---------------------------------------------------------------
    void deserialize(std::istream& is);

    //---------------------------------------------------------------
    /**
     * @brief binary serialization of all non-emtpy buckets
     */
    void serialize(std::ostream& os, int gpuId);


private:
    //---------------------------------------------------------------
    int numGPUs_;
    float maxLoadFactor_;
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
    std::unique_ptr<query_hash_table> queryHashTable_;
};


} // namespace mc


#endif

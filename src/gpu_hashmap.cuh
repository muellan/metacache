#ifndef MC_GPU_HASH_MAP_H_
#define MC_GPU_HASH_MAP_H_

#include <limits>
#include <memory>

#include "config.h"
#include "sequence_batch.cuh"
#include "query_batch.cuh"
#include "taxonomy.h"

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
    class query_hash_table;

public:
    //---------------------------------------------------------------
    using key_type         = Key;
    using value_type       = ValueT;
    using bucket_size_type = loclist_size_t;

    using size_type        = size_t;

    using ranked_lineage = taxonomy::ranked_lineage;

public:
    //---------------------------------------------------------------
    gpu_hashmap();
    ~gpu_hashmap();
    gpu_hashmap(gpu_hashmap&&);
    gpu_hashmap(const gpu_hashmap&) = delete;

private:
    /*************************************************************************//**
    *
    * @brief batch contains features and locations of multiple targets
    *        allocated memory location depends on policy
    *
    *****************************************************************************/
    class feature_batch
    {
    public:
        using counter_type = uint32_t;

        //---------------------------------------------------------------
        feature_batch(counter_type maxFeatures = 0);
        //-----------------------------------------------------
        feature_batch(const feature_batch&) = delete;
        //-----------------------------------------------------
        feature_batch(feature_batch&& other) {
            maxFeatures_ = other.max_features();
            other.max_features(0);

            features_       = other.features();
            values_         = other.values();
            featureCounter_ = other.feature_counter();
        };

        //---------------------------------------------------------------
        ~feature_batch();

        //---------------------------------------------------------------
        counter_type max_features() const noexcept {
            return maxFeatures_;
        }
        //-----------------------------------------------------
        void max_features(counter_type n) noexcept {
            maxFeatures_ = n;
        }
        //---------------------------------------------------------------
        key_type * features() const noexcept {
            return features_;
        }
        //---------------------------------------------------------------
        value_type * values() const noexcept {
            return values_;
        }
        //---------------------------------------------------------------
        counter_type * feature_counter() const noexcept {
            return featureCounter_;
        }

    private:
        counter_type maxFeatures_;

        key_type     * features_;
        value_type   * values_;
        counter_type * featureCounter_;
    };

public:
    //---------------------------------------------------------------
    static constexpr std::size_t
    max_bucket_size() noexcept {
        return std::numeric_limits<bucket_size_type>::max();
    }

    //---------------------------------------------------------------
    size_type key_count() const noexcept;
    //-----------------------------------------------------
    size_type value_count() const noexcept;
    //-----------------------------------------------------
    bool empty() const noexcept {
        return key_count() < 1;
    }
    //---------------------------------------------------------------
    float load_factor() const noexcept;

    //---------------------------------------------------------------
    static constexpr float default_max_load_factor() noexcept {
        return 0.95;
    }
    //-----------------------------------------------------
    float max_load_factor() const noexcept {
        return maxLoadFactor_;
    }

    //---------------------------------------------------------------
    std::vector<key_type> insert(
        const sequence_batch<policy::Host>& seqBatchHost,
        const sketcher& targetSketcher);

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
    friend void write_binary(std::ostream& os, const gpu_hashmap& m) {
        m.serialize(os);
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
    void serialize(std::ostream& os) const;


private:
    //---------------------------------------------------------------
    /**
     * @brief single value hashtable for querying,
     *        which maps feature -> values pointer
     */
    std::unique_ptr<query_hash_table> queryHashTable_;

    float maxLoadFactor_;

    size_t maxBatches_;
    std::vector<sequence_batch<policy::Device>> seqBatches_;
    std::vector<feature_batch> featureBatches_;
};


} // namespace mc


#endif

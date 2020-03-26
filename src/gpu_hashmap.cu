#include <limits>

#include "gpu_hashmap.cuh"
#include "hash_dna.h"
#include "hash_int.h"
#include "sketch_database.h"
#include "gpu_hashmap_operations.cuh"

#include "../dep/warpcore/include/single_value_hash_table.cuh"
#include "../dep/warpcore/include/multi_value_hash_table.cuh"

namespace mc {


/*************************************************************************//**
 *
 * @brief   key -> values hashed multimap
 *          for building metacache db on GPU
 *
 * @details uses warpcore::MultiValueHashTable to map key -> locations
 *
 * @tparam  Key:    key type
 * @tparam  ValueT: value type
 *
 *****************************************************************************/
template<class Key, class ValueT>
class gpu_hashmap<Key,ValueT>::build_hash_table {

    using key_type   = Key;
    using value_type = ValueT;

    using ranked_lineage = taxonomy::ranked_lineage;

    using hash_table_t = warpcore::MultiValueHashTable<
        key_type, value_type,
        // warpcore::defaults::empty_key<key_type>(),       //=0
        key_type(-2),
        warpcore::defaults::tombstone_key<key_type>(),      //=-1
        warpcore::storage::multi_value::DynamicSlabListStore<value_type,34,12,10>>;

    using size_type  = typename hash_table_t::index_type;

public:
    build_hash_table(
        size_type key_capacity,
        size_type value_capacity,
        std::uint64_t maxLocsPerFeature
    ) :
        hashTable_{key_capacity, value_capacity},
        batchSize_{default_batch_size()},
        maxLocsPerFeature_{maxLocsPerFeature},
        seqBatches_{},
        currentSeqBatch_{0}
    {
        std::cerr << "hashtable status: " << hashTable_.pop_status() << "\n";

        std::cerr << "hashtable total: " << (hashTable_.bytes_total() >> 20) << " MB\n";

        seqBatches_.emplace_back(MAX_TARGETS_PER_BATCH, MAX_LENGTH_PER_BATCH);
        seqBatches_.emplace_back(MAX_TARGETS_PER_BATCH, MAX_LENGTH_PER_BATCH);

        cudaStreamCreate(&copyStream_); CUERR
        cudaStreamCreate(&insertStream_); CUERR

        // cudaDeviceSynchronize(); CUERR
    }

    //---------------------------------------------------------------
    static constexpr size_type default_batch_size() noexcept {
        return size_type(1) << 20;
    }
    //-----------------------------------------------------
    size_type batch_size() const noexcept {
        return batchSize_;
    }
    //---------------------------------------------------------------
    float load_factor() noexcept {
        return hashTable_.storage_density();
    }
    //---------------------------------------------------------------
    size_type bucket_count() noexcept {
        return hashTable_.key_capacity();
    }
    //-----------------------------------------------------
    size_type key_count() noexcept {
        return hashTable_.size_keys();
    }
    //-----------------------------------------------------
    size_type location_count() noexcept {
        return hashTable_.size_values();
    }

    //---------------------------------------------------------------
    void insert_async(
        sequence_batch<policy::Host>& seqBatchHost,
        const sketcher& targetSketcher
    ) {
        cudaEventSynchronize(seqBatches_[currentSeqBatch_].event());

        copy_host_to_device_async(
            seqBatchHost, seqBatches_[currentSeqBatch_], copyStream_);

        cudaEventRecord(seqBatchHost.event(), copyStream_); CUERR

        cudaStreamWaitEvent(insertStream_, seqBatchHost.event(), 0); CUERR

        // max 32*4 features => max window size is 128
        constexpr int warpsPerBlock = 1;
        constexpr int threadsPerBlock = 32*warpsPerBlock;
        constexpr int itemsPerThread = 4;

        //TODO increase grid in x and y dim
        const dim3 numBlocks{1024, seqBatches_[currentSeqBatch_].num_targets()};
        insert_features<threadsPerBlock,itemsPerThread>
            <<<numBlocks,threadsPerBlock,0,insertStream_>>>(
            hashTable_,
            seqBatches_[currentSeqBatch_].num_targets(),
            seqBatches_[currentSeqBatch_].target_ids(),
            seqBatches_[currentSeqBatch_].window_offsets(),
            seqBatches_[currentSeqBatch_].sequence(),
            seqBatches_[currentSeqBatch_].sequence_offsets(),
            targetSketcher.kmer_size(),
            targetSketcher.sketch_size(),
            targetSketcher.window_size(),
            targetSketcher.window_stride());

        cudaEventRecord(seqBatches_[currentSeqBatch_].event(), insertStream_);

        // cudaStreamSynchronize(insertStream_); CUERR

        currentSeqBatch_ ^= 1;
    }

    //-----------------------------------------------------
    void wait_until_insert_finished() const {
        cudaStreamSynchronize(insertStream_); CUERR
    }

    //---------------------------------------------------------------
    statistics_accumulator
    location_list_size_statistics() {
        auto priSize = statistics_accumulator{};

        cudaDeviceSynchronize(); CUERR

        key_type * keys = nullptr;
        size_type numKeys = 0;
        hashTable_.retrieve_all_keys(keys, numKeys); CUERR
        cudaMalloc(&keys, numKeys*sizeof(key_type)); CUERR
        hashTable_.retrieve_all_keys(keys, numKeys); CUERR

        size_type  * offsetBuffer_d = nullptr;
        size_type  * offsetBuffer_h = nullptr;
        void       * tmp = nullptr;
        size_type valuesCount = 0;
        size_type tmpSize = 0;

        cudaMalloc    (&offsetBuffer_d, batchSize_*sizeof(size_type)); CUERR
        cudaMallocHost(&offsetBuffer_h, batchSize_*sizeof(size_type)); CUERR

        hashTable_.retrieve(
            keys, batchSize_, offsetBuffer_d,
            nullptr, valuesCount,
            nullptr, tmpSize);
        CUERR

        cudaMalloc(&tmp, tmpSize); CUERR

        const size_type numCycles = numKeys / batchSize_;
        const size_type lastBatchSize = numKeys % batchSize_;

        for(size_type b = 0; b < numCycles; ++b) {
            hashTable_.retrieve(
                keys+b*batchSize_, batchSize_, offsetBuffer_d,
                nullptr, valuesCount,
                tmp, tmpSize);
            CUERR

            cudaMemcpy(offsetBuffer_h, offsetBuffer_d, batchSize_*sizeof(size_type),
                       cudaMemcpyDeviceToHost); CUERR

            priSize += offsetBuffer_h[0];
            for(size_type i = 1; i < batchSize_; ++i)
                priSize += offsetBuffer_h[i] - offsetBuffer_h[i-1];
        }

        hashTable_.retrieve(
            keys+numCycles*batchSize_, lastBatchSize, offsetBuffer_d,
            nullptr, valuesCount,
            tmp, tmpSize);
        CUERR

        cudaMemcpy(offsetBuffer_h, offsetBuffer_d, lastBatchSize*sizeof(size_type),
        cudaMemcpyDeviceToHost); CUERR

        priSize += offsetBuffer_h[0];
        for(size_type i = 1; i < lastBatchSize; ++i)
            priSize += offsetBuffer_h[i] - offsetBuffer_h[i-1];

        cudaFree(keys); CUERR
        cudaFree(offsetBuffer_d); CUERR
        cudaFreeHost(offsetBuffer_h); CUERR
        cudaFree(tmp); CUERR

        return priSize;
    }


private:
    //---------------------------------------------------------------
    void retrieve_and_write_binary(
        std::ostream& os,
        key_type * keys_d,
        size_type batchSize,
        key_type * keyBuffer_h,
        size_type * offsetBuffer_h,
        size_type * offsetBuffer_d,
        std::vector<bucket_size_type>& sizesBuffer,
        value_type *& valueBuffer_h,
        value_type *& valueBuffer_d,
        size_type& valuesCount,
        size_type& valuesAlloc,
        void * tmp,
        size_type tmpSize
    ) {
        hashTable_.retrieve(
            keys_d, batchSize, offsetBuffer_d,
            nullptr, valuesCount,
            tmp, tmpSize);
        CUERR

        if(valuesCount > valuesAlloc) {
            valuesAlloc = valuesCount * 1.1;
            cudaFreeHost(valueBuffer_h); CUERR
            cudaFree    (valueBuffer_d); CUERR
            cudaMallocHost(&valueBuffer_h, valuesAlloc*sizeof(value_type)); CUERR
            cudaMalloc    (&valueBuffer_d, valuesAlloc*sizeof(value_type)); CUERR
        }

        hashTable_.retrieve(
            keys_d, batchSize, offsetBuffer_d,
            valueBuffer_d, valuesCount,
            tmp, tmpSize);
        CUERR

        cudaMemcpy(keyBuffer_h, keys_d, batchSize*sizeof(key_type), cudaMemcpyDeviceToHost);
        write_binary(os, keyBuffer_h, batchSize);

        cudaMemcpy(offsetBuffer_h, offsetBuffer_d, batchSize*sizeof(size_type), cudaMemcpyDeviceToHost);
        sizesBuffer[0] = offsetBuffer_h[0];
        for(size_type i = 1; i < batchSize; ++i)
            sizesBuffer[i] = offsetBuffer_h[i] - offsetBuffer_h[i-1];

        write_binary(os, sizesBuffer.data(), batchSize);

        cudaMemcpy(valueBuffer_h, valueBuffer_d, valuesCount*sizeof(value_type), cudaMemcpyDeviceToHost);
        // write_binary(os, valueBuffer_h, valuesCount);
    }

public:
    //---------------------------------------------------------------
    void serialize(std::ostream& os)
    {
        cudaDeviceSynchronize(); CUERR

        using len_t = std::uint64_t;

        write_binary(os, len_t(key_count()));
        write_binary(os, len_t(location_count()));
        write_binary(os, len_t(batch_size()));

        key_type * keys_d = nullptr;
        size_type numKeys = 0;
        hashTable_.retrieve_all_keys(keys_d, numKeys); CUERR
        cudaMalloc(&keys_d, numKeys*sizeof(key_type)); CUERR
        hashTable_.retrieve_all_keys(keys_d, numKeys); CUERR

        std::vector<bucket_size_type> sizesBuffer(batchSize_);
        key_type   * keyBuffer_h = nullptr;
        size_type  * offsetBuffer_h = nullptr;
        size_type  * offsetBuffer_d = nullptr;
        value_type * valueBuffer_h = nullptr;
        value_type * valueBuffer_d = nullptr;
        void       * tmp = nullptr;
        size_type valuesCount = 0;
        size_type valuesAlloc = 0;
        size_type tmpSize = 0;
        std::vector<std::vector<value_type>> allValues{};

        cudaMallocHost(&keyBuffer_h, batchSize_*sizeof(key_type)); CUERR
        cudaMallocHost(&offsetBuffer_h, batchSize_*sizeof(size_type)); CUERR
        cudaMalloc    (&offsetBuffer_d, batchSize_*sizeof(size_type)); CUERR

        hashTable_.retrieve(
            keys_d, batchSize_, offsetBuffer_d,
            nullptr, valuesCount,
            nullptr, tmpSize);
        CUERR

        cudaMalloc(&tmp, tmpSize); CUERR

        const len_t numCycles = numKeys / batchSize_;
        const len_t lastBatchSize = numKeys % batchSize_;
        allValues.resize(numCycles+1);

        for(len_t b = 0; b < numCycles; ++b) {
            retrieve_and_write_binary(os,
                keys_d+b*batchSize_, batchSize_, keyBuffer_h,
                offsetBuffer_h, offsetBuffer_d, sizesBuffer,
                valueBuffer_h, valueBuffer_d, valuesCount, valuesAlloc,
                tmp, tmpSize);
            allValues[b] = std::vector<value_type>(valueBuffer_h, valueBuffer_h+valuesCount);
        }
        retrieve_and_write_binary(os,
            keys_d+numCycles*batchSize_, lastBatchSize, keyBuffer_h,
            offsetBuffer_h, offsetBuffer_d, sizesBuffer,
            valueBuffer_h, valueBuffer_d, valuesCount, valuesAlloc,
            tmp, tmpSize);
        allValues[numCycles] = std::vector<value_type>(valueBuffer_h, valueBuffer_h+valuesCount);

        for(const auto& someValues : allValues) {
            write_binary(os, someValues.data(), someValues.size());
        }

        cudaFree    (keys_d); CUERR
        cudaFreeHost(keyBuffer_h); CUERR
        cudaFreeHost(offsetBuffer_h); CUERR
        cudaFree    (offsetBuffer_d); CUERR
        cudaFreeHost(valueBuffer_h); CUERR
        cudaFree    (valueBuffer_d); CUERR
        cudaFree    (tmp); CUERR
    }

private:
    hash_table_t hashTable_;

    size_type batchSize_;
    std::uint64_t maxLocsPerFeature_;

    size_t maxBatches_;
    std::vector<sequence_batch<policy::Device>> seqBatches_;
    unsigned currentSeqBatch_;

    cudaStream_t copyStream_;
    cudaStream_t insertStream_;
};



/*************************************************************************//**
 *
 * @brief   key -> values hashed multimap
 *          loads metacache db to GPU to enable queries on GPU
 *
 * @details uses warpcore::SingleValueHashTable to map key -> locations pointer & size
 *          locations are stored in separate array
 *
 * @tparam  Key:    key type
 * @tparam  ValueT: value type
 *
 *****************************************************************************/
template<class Key, class ValueT>
class gpu_hashmap<Key,ValueT>::query_hash_table {

    using key_type   = Key;
    using value_type = std::uint64_t;
    using location_type = ValueT;
    using size_type  = size_t;

    using ranked_lineage = taxonomy::ranked_lineage;

    using hash_table_t = warpcore::SingleValueHashTable<
        key_type, value_type,
        // warpcore::defaults::empty_key<key_type>(),       //=0
        key_type(-2),
        warpcore::defaults::tombstone_key<key_type>(),      //=-1
        warpcore::defaults::probing_scheme_t<key_type, 8>,
        // warpcore::storage::key_value::SoAStore<key_type, value_type>>;
        warpcore::storage::key_value::AoSStore<key_type, value_type>>;

public:
    query_hash_table(size_t capacity) :
        hashTable_(capacity),
        numKeys_(0), numLocations_(0),
        locations_(nullptr)
    {
        size_t freeMemory = 0;
        size_t totalMemory = 0;
        cudaMemGetInfo(&freeMemory, &totalMemory); CUERR
        std::cerr << "freeMemory: " << (freeMemory >> 20) << " MB\n";

        // std::cerr << "capacity: " << hashTable_.capacity() << std::endl;
    }

    //---------------------------------------------------------------
    float load_factor() noexcept {
        return hashTable_.load_factor();
    }
    //---------------------------------------------------------------
    size_type bucket_count() noexcept {
        return hashTable_.capacity();
    }
    //-----------------------------------------------------
    size_type key_count() const noexcept {
        return numKeys_;
    }
    //-----------------------------------------------------
    size_type location_count() const noexcept {
        return numLocations_;
    }

    //---------------------------------------------------------------
    template<class result_type>
    void query_async(
        query_batch<result_type>& batch,
        const sketcher& querySketcher,
        bucket_size_type maxLocationPerFeature,
        bool copyAllHits,
        taxon_rank lowestRank) const
    {
        const cudaStream_t stream = batch.stream();

        // max 32*4 features => max window size is 128
        constexpr int warpsPerBlock = 2;
        constexpr int threadsPerBlock = 32*warpsPerBlock;
        constexpr int itemsPerThread = 4;

        const int numBlocks = (batch.num_queries_device()+warpsPerBlock-1) / warpsPerBlock;
        gpu_hahstable_query<threadsPerBlock,itemsPerThread><<<numBlocks,threadsPerBlock,0,stream>>>(
            hashTable_,
            batch.num_queries_device(),
            batch.sequence_offsets_device(),
            batch.sequences_device(),
            querySketcher.kmer_size(),
            querySketcher.sketch_size(),
            querySketcher.window_size(),
            querySketcher.window_stride(),
            locations_,
            maxLocationPerFeature,
            batch.query_results_device(),
            batch.result_counts_device()
        );
        // batch.sync_streams();
        // CUERR

        batch.compact_sort_and_copy_results_async(copyAllHits);

        batch.generate_and_copy_top_candidates_async(lineages_, lowestRank);

        // batch.sync_result_stream();
        // CUERR
    }

    //---------------------------------------------------------------
    template<class len_t>
    void deserialize(std::istream& is, len_t nkeys, len_t nlocations)
    {
        len_t keyBatchSize = 0;
        read_binary(is, keyBatchSize);

        //TODO tune sizes
        const size_t valBatchSize = 1UL << 20;

        cudaStream_t stream = 0;

        {//load hash table
            //allocate insert buffers
            key_type * h_keyBuffer;
            key_type * d_keyBuffer;
            cudaMallocHost(&h_keyBuffer, keyBatchSize*sizeof(key_type));
            cudaMalloc    (&d_keyBuffer, keyBatchSize*sizeof(key_type));
            uint64_t * h_offsetBuffer;
            uint64_t * d_offsetBuffer;
            cudaMallocHost(&h_offsetBuffer, keyBatchSize*sizeof(uint64_t));
            cudaMalloc    (&d_offsetBuffer, keyBatchSize*sizeof(uint64_t));

            std::vector<bucket_size_type> bsizeBuffer(keyBatchSize);

            using handler_type = warpcore::status_handlers::ReturnStatus;
            using handler_base_type = handler_type::base_type;

            handler_base_type * status;
            cudaMallocManaged(&status, keyBatchSize*sizeof(handler_base_type));
            cudaMemset(status, 0, keyBatchSize*sizeof(handler_base_type));

            const size_type probingLength = hashTable_.capacity();

            uint64_t locsOffset = 0;

            //load full batches
            const len_t numBatches = nkeys / keyBatchSize;
            for(len_t b = 0; b < numBatches; ++b) {
                //load batch
                read_binary(is, h_keyBuffer, keyBatchSize);
                read_binary(is, bsizeBuffer.data(), keyBatchSize);

                for(len_t i = 0; i < keyBatchSize; ++i) {
                    //store offset and size together in 64bit
                    //default is 56bit offset, 8bit size
                    h_offsetBuffer[i] = (locsOffset << sizeof(bucket_size_type)*CHAR_BIT)
                                        + bsizeBuffer[i];

                    locsOffset += bsizeBuffer[i];
                }

                //check status from previous batch
                //implicit sync
                const auto tableStatus = hashTable_.pop_status(stream);
                if(tableStatus.has_any()) {
                    std::cerr << tableStatus << '\n';
                    for(size_t j=0; j<keyBatchSize; ++j) {
                        if(status[j].has_any()) {
                            std::cerr << h_keyBuffer[j] << ' ' << status[j] << '\n';
                        }
                    }
                }

                //insert full batch
                cudaMemcpy(d_keyBuffer, h_keyBuffer, keyBatchSize*sizeof(key_type),
                            cudaMemcpyHostToDevice);
                cudaMemcpy(d_offsetBuffer, h_offsetBuffer, keyBatchSize*sizeof(uint64_t),
                            cudaMemcpyHostToDevice);
                // insert(d_keyBuffer, d_offsetBuffer, keyBatchSize);
                hashTable_.template insert<handler_type>(
                    d_keyBuffer, d_offsetBuffer, keyBatchSize, stream, probingLength, status);
            }

            //load last batch
            const size_t remainingSize = nkeys % keyBatchSize;
            if(remainingSize) {
                //load batch
                read_binary(is, h_keyBuffer, remainingSize);
                read_binary(is, bsizeBuffer.data(), remainingSize);

                for(len_t i = 0; i < remainingSize; ++i) {
                    //store offset and size together in 64bit
                    //default is 56bit offset, 8bit size
                    h_offsetBuffer[i] = (locsOffset << sizeof(bucket_size_type)*CHAR_BIT)
                                        + bsizeBuffer[i];

                    locsOffset += bsizeBuffer[i];
                }

                //check status from previous batch
                //implicit sync
                auto tableStatus = hashTable_.pop_status(stream);
                if(tableStatus.has_any()) {
                    std::cerr << tableStatus << '\n';
                    for(size_t j=0; j<keyBatchSize; ++j) {
                        if(status[j].has_any()) {
                            std::cerr << h_keyBuffer[j] << ' ' << status[j] << '\n';
                        }
                    }
                }

                //insert last batch
                cudaMemcpy(d_keyBuffer, h_keyBuffer, remainingSize*sizeof(key_type),
                            cudaMemcpyHostToDevice);
                cudaMemcpy(d_offsetBuffer, h_offsetBuffer, remainingSize*sizeof(uint64_t),
                            cudaMemcpyHostToDevice);
                // insert(d_keyBuffer, d_offsetBuffer, remainingSize);
                hashTable_.template insert<handler_type>(
                    d_keyBuffer, d_offsetBuffer, remainingSize, stream, probingLength, status);

                //check status from last batch
                //implicit sync
                tableStatus = hashTable_.pop_status(stream);
                if(tableStatus.has_any()) {
                    std::cerr << tableStatus << '\n';
                    for(size_t j=0; j<keyBatchSize; ++j) {
                        if(status[j].has_any()) {
                            std::cerr << h_keyBuffer[j] << ' ' << status[j] << '\n';
                        }
                    }
                }
            }

            cudaFreeHost(h_keyBuffer);
            cudaFree    (d_keyBuffer);
            cudaFreeHost(h_offsetBuffer);
            cudaFree    (d_offsetBuffer);
        }

        {//load locations
            //allocate large memory chunk for all locations,
            //individual buckets will then point into this array
            cudaMalloc(&locations_, nlocations*sizeof(location));

            //allocate buffer
            location * valueBuffers[2];
            cudaMallocHost(&valueBuffers[0], valBatchSize*sizeof(location));
            cudaMallocHost(&valueBuffers[1], valBatchSize*sizeof(location));

            cudaEvent_t events[2];
            cudaEventCreate(&events[0]);
            cudaEventCreate(&events[1]);

            //read batches of locations and copy to device
            auto locsOffset = locations_;
            const len_t numBatches = nlocations / valBatchSize;
            for(len_t i = 0; i < numBatches; ++i) {
                const len_t id = i % 2;
                cudaEventSynchronize(events[id]);
                read_binary(is, valueBuffers[id], valBatchSize);
                cudaMemcpyAsync(locsOffset, valueBuffers[id], valBatchSize*sizeof(location),
                                cudaMemcpyHostToDevice, stream);
                cudaEventRecord(events[id], stream);

                locsOffset += valBatchSize;
            }
            //read remaining locations and copy to device
            const size_t remainingSize = nlocations % valBatchSize;
            const len_t id = numBatches % 2;
            cudaEventSynchronize(events[id]);
            read_binary(is, valueBuffers[id], remainingSize);
            cudaMemcpyAsync(locsOffset, valueBuffers[id], remainingSize*sizeof(location),
                            cudaMemcpyHostToDevice, stream);

            cudaStreamSynchronize(stream);

            cudaFreeHost(valueBuffers[0]);
            cudaFreeHost(valueBuffers[1]);

            cudaEventDestroy(events[0]);
            cudaEventDestroy(events[1]);
        }
        CUERR

        size_t totalSize = hashTable_.capacity() * (sizeof(value_type) + sizeof(value_type))
                         + nlocations*sizeof(location);

        std::cerr << "hashtable total: " << (totalSize >> 20) << " MB\n";

        size_t freeMemory = 0;
        size_t totalMemory = 0;
        cudaMemGetInfo(&freeMemory, &totalMemory); CUERR
        std::cerr << "freeMemory: " << (freeMemory >> 20) << " MB\n";

        numKeys_ = nkeys;
        numLocations_ = nlocations;
    }

    //---------------------------------------------------------------
    void copy_target_lineages_to_gpu(const std::vector<ranked_lineage>& lins) {
        const size_t size = lins.size()*sizeof(ranked_lineage);
        cudaMalloc(&lineages_, size);
        cudaMemcpy(lineages_, lins.data(), size, cudaMemcpyHostToDevice);
    }

private:
    hash_table_t hashTable_;

    size_type numKeys_;
    size_type numLocations_;
    location * locations_;

    ranked_lineage * lineages_;
};



//---------------------------------------------------------------
template<class Key, class ValueT>
gpu_hashmap<Key,ValueT>::gpu_hashmap() :
    maxLoadFactor_(default_max_load_factor())
{}

//-----------------------------------------------------
template<class Key, class ValueT>
gpu_hashmap<Key,ValueT>::~gpu_hashmap() = default;

//-----------------------------------------------------
template<class Key, class ValueT>
gpu_hashmap<Key,ValueT>::gpu_hashmap(gpu_hashmap&&) = default;



//---------------------------------------------------------------
template<class Key, class ValueT>
size_t gpu_hashmap<Key,ValueT>::bucket_count() const noexcept {
    if(buildHashTable_) return buildHashTable_->bucket_count();
    if(queryHashTable_) return queryHashTable_->bucket_count();
    return 0;
}

//---------------------------------------------------------------
template<class Key, class ValueT>
size_t gpu_hashmap<Key,ValueT>::key_count() const noexcept {
    if(buildHashTable_) return buildHashTable_->key_count();
    if(queryHashTable_) return queryHashTable_->key_count();
    return 0;
}

//-----------------------------------------------------
template<class Key, class ValueT>
size_t gpu_hashmap<Key,ValueT>::value_count() const noexcept {
    if(buildHashTable_) return buildHashTable_->location_count();
    if(queryHashTable_) return queryHashTable_->location_count();
    return 0;
}

//---------------------------------------------------------------
template<class Key, class ValueT>
float gpu_hashmap<Key,ValueT>::load_factor() const noexcept {
    if(buildHashTable_) return buildHashTable_->load_factor();
    if(queryHashTable_) return queryHashTable_->load_factor();
    return -1;
}

//---------------------------------------------------------------
template<class Key, class ValueT>
statistics_accumulator gpu_hashmap<Key,ValueT>::location_list_size_statistics() const {
    if(buildHashTable_) return buildHashTable_->location_list_size_statistics();
    return statistics_accumulator{};
}



//---------------------------------------------------------------
template<class Key, class ValueT>
void gpu_hashmap<Key,ValueT>::initialize_hash_table(
    std::uint64_t maxLocsPerFeature)
{
    size_t freeMemory = 0;
    size_t totalMemory = 0;
    cudaMemGetInfo(&freeMemory, &totalMemory); CUERR
    std::cerr << "freeMemory: " << (freeMemory >> 20) << " MB\n";

    // keep 4 GB of memory free aside from hash table
    const size_t tableMemory = freeMemory - (1ULL << 32);

    constexpr size_t valueSize = sizeof(ValueT);

    const size_t keyCapacity   = tableMemory *  1/11 / (2*valueSize);
    const size_t valueCapacity = tableMemory * 10/11 / valueSize;

    std::cerr << "allocate hashtable for " << keyCapacity << " keys"
                                   " and " << valueCapacity << " values\n";
    buildHashTable_ = std::make_unique<build_hash_table>(
                          keyCapacity, valueCapacity, maxLocsPerFeature);

    cudaMemGetInfo(&freeMemory, &totalMemory); CUERR
    std::cerr << "freeMemory: " << (freeMemory >> 20) << " MB\n";
}


//---------------------------------------------------------------
template<class Key, class ValueT>
void gpu_hashmap<Key,ValueT>::insert(
    sequence_batch<policy::Host>& seqBatchHost,
    const sketcher& targetSketcher)
{
    buildHashTable_->insert_async(
        seqBatchHost,
        targetSketcher);
}
//-----------------------------------------------------
template<class Key, class ValueT>
void gpu_hashmap<Key,ValueT>::wait_until_insert_finished() const {
    buildHashTable_->wait_until_insert_finished();
}


//---------------------------------------------------------------
template<class Key, class ValueT>
void gpu_hashmap<Key,ValueT>::query_async(
    query_batch<value_type>& batch,
    const sketcher& querySketcher,
    bucket_size_type maxLocationPerFeature,
    bool copyAllHits,
    taxon_rank lowestRank) const
{
    queryHashTable_->query_async(
        batch,
        querySketcher,
        maxLocationPerFeature,
        copyAllHits,
        lowestRank);
}


//---------------------------------------------------------------
template<class Key, class ValueT>
void gpu_hashmap<Key,ValueT>::deserialize(std::istream& is)
{
    using len_t = std::uint64_t;

    len_t nkeys = 0;
    read_binary(is, nkeys);
    len_t nvalues = 0;
    read_binary(is, nvalues);

    std::cerr << "nkeys: " << nkeys << " nvalues: " << nvalues << std::endl;

    if(nkeys > 0) {
        //initialize hash table
        queryHashTable_ = std::make_unique<query_hash_table>(nkeys/maxLoadFactor_);

        queryHashTable_->deserialize(is, nkeys, nvalues);
    }
}


//---------------------------------------------------------------
/**
* @brief binary serialization of all non-emtpy buckets
*/
template<class Key, class ValueT>
void gpu_hashmap<Key,ValueT>::serialize(std::ostream& os) const
{
    buildHashTable_->serialize(os);
}


//---------------------------------------------------------------
template<class Key, class ValueT>
void gpu_hashmap<Key,ValueT>::copy_target_lineages_to_gpu(
    const std::vector<ranked_lineage>& lins)
{
    queryHashTable_->copy_target_lineages_to_gpu(lins);
}



//---------------------------------------------------------------
template class gpu_hashmap<kmer_type, location>;

} // namespace mc

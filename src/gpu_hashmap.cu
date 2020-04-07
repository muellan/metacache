#include <limits>

#include "gpu_hashmap.cuh"
#include "hash_dna.h"
#include "hash_int.h"
#include "sketch_database.h"
#include "gpu_hashmap_operations.cuh"

#include "../dep/warpcore/include/single_value_hash_table.cuh"
#include "../dep/warpcore/include/bucket_list_hash_table.cuh"

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

    using hash_table_t = warpcore::BucketListHashTable<
        key_type, value_type,
        // warpcore::defaults::empty_key<key_type>(),       //=0
        key_type(-2),
        warpcore::defaults::tombstone_key<key_type>(),      //=-1
        warpcore::storage::multi_value::BucketListStore<value_type,40,8,8>
    >;

    using size_type  = typename hash_table_t::index_type;
    using status_type  = typename warpcore::Status;

public:
    build_hash_table(
        size_type key_capacity,
        size_type value_capacity,
        std::uint64_t maxLocsPerFeature
    ) :
        hashTable_{key_capacity, value_capacity,
            warpcore::defaults::seed<key_type>(),   // seed
            1.051, 1, max_bucket_size(),            // grow factor, min & max bucket size
            maxLocsPerFeature},                     // max values per key
        batchSize_{default_batch_size()},
        seqBatches_{},
        currentSeqBatch_{0}
    {
        std::cerr << "hashtable status: " << hashTable_.pop_status() << "\n";

        seqBatches_.emplace_back(MAX_TARGETS_PER_BATCH, MAX_LENGTH_PER_BATCH);
        seqBatches_.emplace_back(MAX_TARGETS_PER_BATCH, MAX_LENGTH_PER_BATCH);

        cudaStreamCreate(&copyStream_); CUERR
        cudaStreamCreate(&insertStream_); CUERR

        // cudaDeviceSynchronize(); CUERR
    }

    //---------------------------------------------------------------
    bool validate() {
        if(hashTable_.peek_status() - status_type::max_values_for_key_reached())
            return false;
        return true;
    }

    //---------------------------------------------------------------
    status_type pop_status() {
        return hashTable_.pop_status();
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
        return hashTable_.num_keys();
    }
    //-----------------------------------------------------
    size_type location_count() noexcept {
        return hashTable_.num_values();
    }

    //---------------------------------------------------------------
    void insert_async(
        sequence_batch<policy::Host>& seqBatchHost,
        const sketcher& targetSketcher
    ) {
        // wait for previous insert of current batch to finish
        cudaStreamWaitEvent(copyStream_, seqBatches_[currentSeqBatch_].event(), 0); CUERR

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

        size_type  * numValuesBuffer_d = nullptr;
        size_type  * numValuesBuffer_h = nullptr;
        cudaMalloc    (&numValuesBuffer_d, batchSize_*sizeof(size_type)); CUERR
        cudaMallocHost(&numValuesBuffer_h, batchSize_*sizeof(size_type)); CUERR

        size_type * valuesCountPtr = nullptr;
        cudaMallocHost(&valuesCountPtr, sizeof(size_type)); CUERR
        *valuesCountPtr = 0;

        const size_type numCycles = numKeys / batchSize_;
        const size_type lastBatchSize = numKeys % batchSize_;

        for(size_type b = 0; b < numCycles; ++b) {
            hashTable_.num_values(
                keys+b*batchSize_, batchSize_, *valuesCountPtr,
                numValuesBuffer_d);
            CUERR

            cudaMemcpy(numValuesBuffer_h, numValuesBuffer_d, batchSize_*sizeof(size_type),
                       cudaMemcpyDeviceToHost); CUERR

            for(size_type i = 0; i < batchSize_; ++i)
                priSize += numValuesBuffer_h[i];
        }
        if(lastBatchSize) {
            hashTable_.num_values(
                keys+numCycles*batchSize_, lastBatchSize, *valuesCountPtr,
                numValuesBuffer_d);
            CUERR

            cudaMemcpy(numValuesBuffer_h, numValuesBuffer_d, lastBatchSize*sizeof(size_type),
                    cudaMemcpyDeviceToHost); CUERR

            for(size_type i = 0; i < lastBatchSize; ++i)
                priSize += numValuesBuffer_h[i];
        }

        cudaFree(keys); CUERR
        cudaFree(numValuesBuffer_d); CUERR
        cudaFreeHost(numValuesBuffer_h); CUERR
        cudaFreeHost(valuesCountPtr); CUERR

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
        size_type& valuesAlloc
    ) {
        // get valuesCount
        hashTable_.retrieve(
            keys_d, batchSize,
            offsetBuffer_d, offsetBuffer_d+1,
            nullptr, valuesCount);
        CUERR

        // reallocate if buffers to small
        if(valuesCount > valuesAlloc) {
            valuesAlloc = valuesCount * 1.1;
            cudaFreeHost(valueBuffer_h); CUERR
            cudaFree    (valueBuffer_d); CUERR
            cudaMallocHost(&valueBuffer_h, valuesAlloc*sizeof(value_type)); CUERR
            cudaMalloc    (&valueBuffer_d, valuesAlloc*sizeof(value_type)); CUERR
        }

        // get values
        hashTable_.retrieve(
            keys_d, batchSize,
            offsetBuffer_d, offsetBuffer_d+1,
            valueBuffer_d, valuesCount);
        CUERR

        cudaMemcpy(keyBuffer_h, keys_d, batchSize*sizeof(key_type), cudaMemcpyDeviceToHost);
        write_binary(os, keyBuffer_h, batchSize);

        cudaMemcpy(offsetBuffer_h, offsetBuffer_d, (batchSize+1)*sizeof(size_type), cudaMemcpyDeviceToHost);
        for(size_type i = 0; i < batchSize; ++i)
            sizesBuffer[i] = offsetBuffer_h[i+1] - offsetBuffer_h[i];

        write_binary(os, sizesBuffer.data(), batchSize);

        cudaMemcpy(valueBuffer_h, valueBuffer_d, valuesCount*sizeof(value_type), cudaMemcpyDeviceToHost);
        write_binary(os, valueBuffer_h, valuesCount);
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
        // get numKeys
        hashTable_.retrieve_all_keys(keys_d, numKeys); CUERR
        cudaMalloc(&keys_d, numKeys*sizeof(key_type)); CUERR
        // get keys
        hashTable_.retrieve_all_keys(keys_d, numKeys); CUERR

        std::vector<bucket_size_type> sizesBuffer(batchSize_);
        key_type   * keyBuffer_h = nullptr;
        size_type  * offsetBuffer_h = nullptr;
        size_type  * offsetBuffer_d = nullptr;
        value_type * valueBuffer_h = nullptr;
        value_type * valueBuffer_d = nullptr;

        cudaMallocHost(&keyBuffer_h, batchSize_*sizeof(key_type)); CUERR
        cudaMallocHost(&offsetBuffer_h, (batchSize_+1)*sizeof(size_type)); CUERR
        cudaMalloc    (&offsetBuffer_d, (batchSize_+1)*sizeof(size_type)); CUERR

        size_type valuesAlloc = 0;
        size_type * valuesCountPtr = nullptr;
        cudaMallocHost(&valuesCountPtr, sizeof(size_type)); CUERR
        *valuesCountPtr = 0;

        const len_t numCycles = numKeys / batchSize_;
        const len_t lastBatchSize = numKeys % batchSize_;

        for(len_t b = 0; b < numCycles; ++b) {
            retrieve_and_write_binary(os,
                keys_d+b*batchSize_, batchSize_, keyBuffer_h,
                offsetBuffer_h, offsetBuffer_d, sizesBuffer,
                valueBuffer_h, valueBuffer_d, *valuesCountPtr, valuesAlloc);
        }
        if(lastBatchSize) {
            retrieve_and_write_binary(os,
                keys_d+numCycles*batchSize_, lastBatchSize, keyBuffer_h,
                offsetBuffer_h, offsetBuffer_d, sizesBuffer,
                valueBuffer_h, valueBuffer_d, *valuesCountPtr, valuesAlloc);
        }

        cudaFree    (keys_d); CUERR
        cudaFreeHost(keyBuffer_h); CUERR
        cudaFreeHost(offsetBuffer_h); CUERR
        cudaFree    (offsetBuffer_d); CUERR
        cudaFreeHost(valuesCountPtr); CUERR
        if(valuesAlloc) {
            cudaFreeHost(valueBuffer_h); CUERR
            cudaFree    (valueBuffer_d); CUERR
        }
    }

private:
    hash_table_t hashTable_;

    size_type batchSize_;

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
    auto pop_status() {
        return hashTable_.pop_status();
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

private:
    //---------------------------------------------------------------
    template<class LenT, class Status>
    LenT deserialize_batch_of_buckets(
        std::istream& is,
        key_type * h_keyBuffer,
        key_type * d_keyBuffer,
        uint64_t * h_offsetBuffer,
        uint64_t * d_offsetBuffer,
        std::vector<bucket_size_type>& bsizeBuffer,
        LenT batchSize,
        location * valueBuffers[2],
        cudaEvent_t events[2],
        LenT valBatchSize,
        location * d_values,
        uint64_t valuesOffset,
        Status * status,
        cudaStream_t stream)
    {
        using len_t = LenT;
        using handler_type = warpcore::status_handlers::ReturnStatus;

        const size_type probingLength = hashTable_.capacity();

        auto batchValuesOffset = valuesOffset;

        //load batch
        read_binary(is, h_keyBuffer, batchSize);
        read_binary(is, bsizeBuffer.data(), batchSize);

        for(len_t i = 0; i < batchSize; ++i) {
            //store offset and size together in 64bit
            //default is 56bit offset, 8bit size
            h_offsetBuffer[i] = (valuesOffset << sizeof(bucket_size_type)*CHAR_BIT)
                                + bsizeBuffer[i];

            valuesOffset += bsizeBuffer[i];
        }

        //check status from previous batch
        //implicit sync
        const auto tableStatus = hashTable_.pop_status(stream);
        if(tableStatus.has_any()) {
            std::cerr << tableStatus << '\n';
            for(size_t j=0; j<batchSize; ++j) {
                if(status[j].has_any()) {
                    std::cerr << h_keyBuffer[j] << ' ' << status[j] << '\n';
                }
            }
        }

        //insert batch
        cudaMemcpy(d_keyBuffer, h_keyBuffer, batchSize*sizeof(key_type),
                    cudaMemcpyHostToDevice);
        cudaMemcpy(d_offsetBuffer, h_offsetBuffer, batchSize*sizeof(uint64_t),
                    cudaMemcpyHostToDevice);
        // insert(d_keyBuffer, d_offsetBuffer, batchSize);
        hashTable_.template insert<handler_type>(
            d_keyBuffer, d_offsetBuffer, batchSize, stream, probingLength, status);


        std::uint64_t batchValuesCount = valuesOffset - batchValuesOffset;
        //read batches of locations and copy to device
        const len_t numBatches = batchValuesCount / valBatchSize;
        const size_t remainingSize = batchValuesCount % valBatchSize;

        d_values += batchValuesOffset;

        for(len_t i = 0; i < numBatches; ++i) {
            const len_t id = i % 2;
            cudaEventSynchronize(events[id]);
            read_binary(is, valueBuffers[id], valBatchSize);
            cudaMemcpyAsync(d_values, valueBuffers[id], valBatchSize*sizeof(location),
                            cudaMemcpyHostToDevice, stream);
            cudaEventRecord(events[id], stream);

            d_values += valBatchSize;
        }
        //read remaining locations and copy to device
        const len_t id = numBatches % 2;
        cudaEventSynchronize(events[id]);
        read_binary(is, valueBuffers[id], remainingSize);
        cudaMemcpyAsync(d_values, valueBuffers[id], remainingSize*sizeof(location),
                        cudaMemcpyHostToDevice, stream);


        return batchValuesCount;
    }

public:
    //---------------------------------------------------------------
    template<class len_t>
    void deserialize(std::istream& is, len_t nkeys, len_t nlocations)
    {
        len_t batchSize = 0;
        read_binary(is, batchSize);

        //TODO tune sizes
        const len_t valBatchSize = 1UL << 20;

        cudaStream_t stream = 0;

        //allocate large memory chunk for all locations,
        //individual buckets will then point into this array
        cudaMalloc(&locations_, nlocations*sizeof(location)); CUERR
        uint64_t locsOffset = 0;

        {//load hash table
            //allocate insert buffers
            key_type * h_keyBuffer;
            key_type * d_keyBuffer;
            cudaMallocHost(&h_keyBuffer, batchSize*sizeof(key_type));
            cudaMalloc    (&d_keyBuffer, batchSize*sizeof(key_type));
            uint64_t * h_offsetBuffer;
            uint64_t * d_offsetBuffer;
            cudaMallocHost(&h_offsetBuffer, batchSize*sizeof(uint64_t));
            cudaMalloc    (&d_offsetBuffer, batchSize*sizeof(uint64_t));
            location * valueBuffers[2];
            cudaMallocHost(&valueBuffers[0], valBatchSize*sizeof(location));
            cudaMallocHost(&valueBuffers[1], valBatchSize*sizeof(location));
            cudaEvent_t events[2];
            cudaEventCreate(&events[0]);
            cudaEventCreate(&events[1]);
            CUERR

            std::vector<bucket_size_type> bsizeBuffer(batchSize);

            using handler_type = warpcore::status_handlers::ReturnStatus;
            using handler_base_type = handler_type::base_type;

            handler_base_type * status;
            cudaMallocManaged(&status, batchSize*sizeof(handler_base_type));
            cudaMemset(status, 0, batchSize*sizeof(handler_base_type));

            //load full batches
            const len_t numBatches = nkeys / batchSize;
            for(len_t b = 0; b < numBatches; ++b) {
                auto batchValuesCount = deserialize_batch_of_buckets(is,
                    h_keyBuffer, d_keyBuffer, h_offsetBuffer, d_offsetBuffer,
                    bsizeBuffer, batchSize,
                    valueBuffers, events, valBatchSize, locations_, locsOffset,
                    status, stream);

                locsOffset += batchValuesCount;
            }

            //load last batch
            const size_t remainingSize = nkeys % batchSize;
            if(remainingSize) {
                auto batchValuesCount = deserialize_batch_of_buckets(is,
                    h_keyBuffer, d_keyBuffer, h_offsetBuffer, d_offsetBuffer,
                    bsizeBuffer, remainingSize,
                    valueBuffers, events, valBatchSize, locations_, locsOffset,
                    status, stream);

                locsOffset += batchValuesCount;

                //check status from last batch
                //implicit sync
                const auto tableStatus = hashTable_.pop_status(stream);
                if(tableStatus.has_any()) {
                    std::cerr << tableStatus << '\n';
                    for(size_t j=0; j<batchSize; ++j) {
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

            cudaFreeHost(valueBuffers[0]);
            cudaFreeHost(valueBuffers[1]);
            cudaEventDestroy(events[0]);
            cudaEventDestroy(events[1]);
            CUERR
        }

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
    maxLoadFactor_(default_max_load_factor()),
    valid_(true)
{}

//-----------------------------------------------------
template<class Key, class ValueT>
gpu_hashmap<Key,ValueT>::~gpu_hashmap() = default;

//-----------------------------------------------------
template<class Key, class ValueT>
gpu_hashmap<Key,ValueT>::gpu_hashmap(gpu_hashmap&&) = default;



//---------------------------------------------------------------
template<class Key, class ValueT>
bool gpu_hashmap<Key,ValueT>::valid() const noexcept {
    return valid_;
}

//---------------------------------------------------------------
template<class Key, class ValueT>
void gpu_hashmap<Key,ValueT>::pop_status() {
    if(buildHashTable_)
        std::cerr << "hashtable status: " << buildHashTable_->pop_status() << "\n";
    else if(queryHashTable_)
        std::cerr << "hashtable status: " << queryHashTable_->pop_status() << "\n";
}

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

    const size_t keyCapacity   = tableMemory *  2/13 / (2*valueSize);
    const size_t valueCapacity = tableMemory * 11/13 / valueSize;

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
    if(valid_ && buildHashTable_->validate()) {
        buildHashTable_->insert_async(
            seqBatchHost,
            targetSketcher);
    }
    else {
        valid_ = false;
    }
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

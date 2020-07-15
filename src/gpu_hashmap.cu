#include <limits>

#include "gpu_hashmap.cuh"
#include "hash_dna.h"
#include "hash_int.h"
#include "database.h"
#include "gpu_hashmap_operations.cuh"
#include "stat_combined.cuh"

#include "../dep/warpcore/include/single_value_hash_table.cuh"
#include "../dep/warpcore/include/bucket_list_hash_table.cuh"

namespace mc {


template<class SizeT, class BucketSizeT>
__global__
void calculate_sizes_kernel(SizeT * d_offsets, BucketSizeT * d_sizes, SizeT batchSize)
{
    const auto tid = blockDim.x * blockIdx.x + threadIdx.x;

    if(tid < batchSize) {
        d_sizes[tid] = d_offsets[tid+1] - d_offsets[tid];
    }
}


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
        warpcore::storage::multi_value::BucketListStore<
            value_type,40,bucket_size_bits(),bucket_size_bits()
        >
    >;

    using size_type  = typename hash_table_t::index_type;
    using status_type  = typename warpcore::Status;

public:
    build_hash_table(
        size_type key_capacity,
        size_type value_capacity,
        std::uint64_t maxLocationsPerFeature
    ) :
        hashTable_{key_capacity, value_capacity,
            warpcore::defaults::seed<key_type>(),   // seed
            1.051, 1, max_bucket_size(),            // grow factor, min & max bucket size
            // 1.075, 3, 26,                           // grow factor, min & max bucket size
            maxLocationsPerFeature},                // max values per key
        batchSize_{default_batch_size()},
        seqBatches_{},
        currentSeqBatch_{0}
    {
        std::cerr << "hashtable status: " << hashTable_.pop_status() << "\n";

        seqBatches_.emplace_back(MAX_TARGETS_PER_BATCH, MAX_LENGTH_PER_BATCH);
        seqBatches_.emplace_back(MAX_TARGETS_PER_BATCH, MAX_LENGTH_PER_BATCH);

        cudaStreamCreate(&copyStream_); CUERR
        cudaStreamCreate(&insertStream_); CUERR
        cudaStreamCreate(&statusStream_); CUERR

        // cudaDeviceSynchronize(); CUERR
    }

    //---------------------------------------------------------------
    bool validate() {
        if(hashTable_.peek_status(statusStream_) - status_type::max_values_for_key_reached())
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
    size_type bucket_count() const noexcept {
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

        constexpr int maxSketchSize = 16;

        // max 32*4 characters per warp, so max window size is 128
        if(targetSketcher.window_size() <= 128 && targetSketcher.sketch_size() <= maxSketchSize) {
            constexpr int warpsPerBlock = 2;
            constexpr int threadsPerBlock = 32*warpsPerBlock;

            const dim3 numBlocks{1024, seqBatches_[currentSeqBatch_].num_targets()};
            insert_features<threadsPerBlock,maxSketchSize>
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
        }
        else {
            std::cerr << "Max window size is 128!\n";
            std::cerr << "Max sketch size is " << maxSketchSize << "\n";
        }

        cudaEventRecord(seqBatches_[currentSeqBatch_].event(), insertStream_); CUERR

        // cudaStreamSynchronize(insertStream_); CUERR

        currentSeqBatch_ ^= 1;
    }

    //-----------------------------------------------------
    void wait_until_insert_finished() const {
        cudaStreamSynchronize(insertStream_); CUERR
    }

    //---------------------------------------------------------------
    statistics_accumulator_gpu<policy::Host>
    location_list_size_statistics()
    {
        cudaDeviceSynchronize(); CUERR

        key_type * keys = nullptr;
        size_type numKeys = hashTable_.num_keys(); CUERR
        cudaMalloc(&keys, numKeys*sizeof(key_type)); CUERR
        hashTable_.retrieve_all_keys(keys, numKeys); CUERR

        size_type  * numValuesBuffer_d = nullptr;
        cudaMalloc(&numValuesBuffer_d, batchSize_*sizeof(size_type)); CUERR

        size_type * valuesCountPtr = nullptr;
        cudaMallocHost(&valuesCountPtr, sizeof(size_type)); CUERR
        *valuesCountPtr = 0;

        statistics_accumulator_gpu<policy::Device> accumulator_d{};

        const size_type numCycles = numKeys / batchSize_;
        const size_type lastBatchSize = numKeys % batchSize_;

        for(size_type b = 0; b < numCycles; ++b) {
            hashTable_.num_values(
                keys+b*batchSize_, batchSize_, *valuesCountPtr,
                numValuesBuffer_d);
            CUERR

            accumulator_d.accumulate(numValuesBuffer_d, batchSize_);
        }
        if(lastBatchSize) {
            hashTable_.num_values(
                keys+numCycles*batchSize_, lastBatchSize, *valuesCountPtr,
                numValuesBuffer_d);
            CUERR

            accumulator_d.accumulate(numValuesBuffer_d, lastBatchSize);
        }

        cudaFree(keys); CUERR
        cudaFree(numValuesBuffer_d); CUERR
        cudaFreeHost(valuesCountPtr); CUERR

        statistics_accumulator_gpu<policy::Host> accumulator_h{};
        accumulator_h = accumulator_d;

        return accumulator_h;
    }


private:
    class retrieval_buffer {
    public:
        retrieval_buffer(size_type batchSize) :
            valuesAlloc_{0},
            d_values_{nullptr},
            h_values_{nullptr}
        {
            cudaMallocHost(&valuesCountPtr_, sizeof(size_type)); CUERR
            valuesCountPtr_[0] = 0;

            cudaMallocHost(&h_keys_, batchSize*sizeof(key_type)); CUERR
            cudaMalloc    (&d_offsets_, (batchSize+1)*sizeof(size_type)); CUERR
            cudaMalloc    (&d_sizes_, batchSize*sizeof(bucket_size_type)); CUERR
            cudaMallocHost(&h_sizes_, batchSize*sizeof(bucket_size_type)); CUERR

            cudaStreamCreate(&stream_);
        }

        ~retrieval_buffer() {
            cudaFreeHost(valuesCountPtr_); CUERR
            cudaFreeHost(h_keys_); CUERR
            cudaFree    (d_offsets_); CUERR
            cudaFree    (d_sizes_); CUERR
            cudaFreeHost(h_sizes_); CUERR
            if(valuesAlloc_) {
                cudaFree    (d_values_); CUERR
                cudaFreeHost(h_values_); CUERR
            }

            cudaStreamDestroy(stream_);
        }

        size_type * values_count() const noexcept { return valuesCountPtr_; }
        size_type * d_offsets() const noexcept { return d_offsets_; }
        key_type * h_keys() const noexcept { return h_keys_; }
        bucket_size_type * d_sizes() const noexcept { return d_sizes_; }
        bucket_size_type * h_sizes() const noexcept { return h_sizes_; }
        value_type * d_values() const noexcept { return d_values_; }
        value_type * h_values() const noexcept { return h_values_; }
        cudaStream_t stream() const noexcept { return stream_; }

        void resize() {
            if(*values_count() > valuesAlloc_) {
                valuesAlloc_ = *values_count() * 1.1;
                cudaFreeHost(h_values()); CUERR
                cudaFree    (d_values()); CUERR
                cudaMallocHost(&h_values_, valuesAlloc_*sizeof(value_type)); CUERR
                cudaMalloc    (&d_values_, valuesAlloc_*sizeof(value_type)); CUERR
            }
        }

    private:
        size_type valuesAlloc_;
        size_type * valuesCountPtr_;

        key_type   * h_keys_;
        size_type  * d_offsets_;
        bucket_size_type * d_sizes_;
        bucket_size_type * h_sizes_;
        value_type * d_values_;
        value_type * h_values_;

        cudaStream_t stream_;
    };

    //---------------------------------------------------------------
    void retrieve_and_write_binary(
        std::ostream& os,
        key_type * d_keys,
        retrieval_buffer& buffer,
        size_type batchSize,
        std::mutex& mtx
    ) {
        // get valuesCount
        hashTable_.num_values(
            d_keys, batchSize,
            *(buffer.values_count()),
            buffer.d_offsets()+1,
            buffer.stream());
        cudaStreamSynchronize(buffer.stream()); CUERR

        // reallocate if buffers to small
        buffer.resize();

        // get values
        hashTable_.retrieve(
            d_keys, batchSize,
            buffer.d_offsets(), buffer.d_offsets()+1,
            buffer.d_values(), *(buffer.values_count()),
            buffer.stream());

        calculate_sizes_kernel<<<SDIV(batchSize, MAXBLOCKSIZE), MAXBLOCKSIZE, 0, buffer.stream()>>>(
            buffer.d_offsets(), buffer.d_sizes(), batchSize);

        cudaMemcpyAsync(buffer.h_keys(), d_keys, batchSize*sizeof(key_type),
            cudaMemcpyDeviceToHost, buffer.stream());
        cudaMemcpyAsync( buffer.h_sizes(),  buffer.d_sizes(), batchSize*sizeof(bucket_size_type),
            cudaMemcpyDeviceToHost, buffer.stream());
        cudaMemcpyAsync( buffer.h_values(), buffer. d_values(), *buffer.values_count()*sizeof(value_type),
            cudaMemcpyDeviceToHost, buffer.stream());

        cudaStreamSynchronize(buffer.stream()); CUERR

        const auto tableStatus = hashTable_.pop_status(buffer.stream());
        if(tableStatus.has_any())
            std::cerr << tableStatus << '\n';

        std::lock_guard<std::mutex> lock(mtx);
        write_binary(os, buffer.h_keys(), batchSize);
        write_binary(os, buffer.h_sizes(), batchSize);
        write_binary(os, buffer.h_values(), *buffer.values_count());
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

        size_type numKeys = hashTable_.num_keys(); CUERR
        // allocate buffers
        key_type * d_keys;
        cudaMalloc(&d_keys, numKeys*sizeof(key_type)); CUERR
        retrieval_buffer buffer0(batchSize_);
        retrieval_buffer buffer1(batchSize_);
        // get keys
        hashTable_.retrieve_all_keys(d_keys, numKeys); CUERR

        const len_t numCycles = numKeys / batchSize_;
        const len_t lastBatchSize = numKeys % batchSize_;
        std::mutex mtx;

        int gpuId = -1;
        cudaGetDevice(&gpuId);

        auto retriever0 = std::async(std::launch::async, [&] {
            cudaSetDevice(gpuId);

            for(len_t b = 0; b < numCycles; b+=2) {
                retrieve_and_write_binary(os,
                    d_keys + b * batchSize_,
                    buffer0, batchSize_, mtx);
            }
        });
        auto retriever1 = std::async(std::launch::async, [&] {
            cudaSetDevice(gpuId);

            for(len_t b = 1; b < numCycles; b+=2) {
                retrieve_and_write_binary(os,
                    d_keys + b * batchSize_,
                    buffer1, batchSize_, mtx);
            }
        });

        retriever0.get();
        retriever1.get();

        if(lastBatchSize) {
            retrieve_and_write_binary(os,
                d_keys + numCycles * batchSize_,
                buffer0, lastBatchSize, mtx);
        }

        cudaFree(d_keys); CUERR
    }

private:
    hash_table_t hashTable_;

    size_type batchSize_;

    size_t maxBatches_;
    std::vector<sequence_batch<policy::Device>> seqBatches_;
    unsigned currentSeqBatch_;

    cudaStream_t copyStream_;
    cudaStream_t insertStream_;
    cudaStream_t statusStream_;
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
        locations_(nullptr),
        lineages_(nullptr)
    {}

    //---------------------------------------------------------------
    auto pop_status() {
        return hashTable_.pop_status();
    }

    //---------------------------------------------------------------
    float load_factor() noexcept {
        return hashTable_.load_factor();
    }
    //---------------------------------------------------------------
    size_type bucket_count() const noexcept {
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
    ranked_lineage * lineages() const noexcept {
        return lineages_;
    }

    /*************************************************************************//**
    *
    * @brief   query all windows in batch using one warp per window
    *
    * @details saves sketches to gpu memory in case of multi-gpu query
    *
    *****************************************************************************/
    void query_sequences_async(
        uint32_t numWindows,
        const typename query_batch<location_type>::query_gpu_data& gpuData,
        const sketcher& querySketcher,
        bucket_size_type maxLocationsPerFeature) const
    {
        constexpr int maxSketchSize = 16;

        // max 32*4 characters per warp, so max window size is 128
        if(querySketcher.window_size() <= 128 && querySketcher.sketch_size() <= maxSketchSize) {
            constexpr int warpsPerBlock = 2;
            constexpr int threadsPerBlock = 32*warpsPerBlock;

            const int numBlocks = (numWindows+warpsPerBlock-1) / warpsPerBlock;
            gpu_hahstable_query<threadsPerBlock,maxSketchSize>
                <<<numBlocks,threadsPerBlock,0,gpuData.workStream_>>>(
                hashTable_,
                numWindows,
                gpuData.sequenceOffsets_,
                gpuData.sequences_,
                gpuData.sketches_,
                querySketcher.kmer_size(),
                querySketcher.sketch_size(),
                querySketcher.window_size(),
                querySketcher.window_stride(),
                locations_,
                maxLocationsPerFeature,
                gpuData.queryResults_,
                gpuData.resultCounts_
            );
        }
        else {
            std::cerr << "Max window size is 128!\n";
            std::cerr << "Max sketch size is " << maxSketchSize << "\n";
        }
    }

    /*************************************************************************//**
    *
    * @brief   query sketches of all windows in batch using one warp per window
    *
    *****************************************************************************/
    void query_sketches_async(
        uint32_t numWindows,
        const typename query_batch<location_type>::query_gpu_data& gpuData,
        const sketcher& querySketcher,
        bucket_size_type maxLocationsPerFeature) const
    {
        constexpr int maxSketchSize = 16;

        // max 32*4 characters per warp, so max window size is 128
        if(querySketcher.window_size() <= 128 && querySketcher.sketch_size() <= maxSketchSize) {
            constexpr int warpsPerBlock = 2;
            constexpr int threadsPerBlock = 32*warpsPerBlock;

            const int numBlocks = (numWindows+warpsPerBlock-1) / warpsPerBlock;
            gpu_hahstable_query<threadsPerBlock,maxSketchSize>
                <<<numBlocks,threadsPerBlock,0,gpuData.workStream_>>>(
                hashTable_,
                numWindows,
                gpuData.sketches_,
                querySketcher.kmer_size(),
                querySketcher.sketch_size(),
                querySketcher.window_size(),
                querySketcher.window_stride(),
                locations_,
                maxLocationsPerFeature,
                gpuData.queryResults_,
                gpuData.resultCounts_
            );
        }
        else {
            std::cerr << "Max window size is 128!\n";
            std::cerr << "Max sketch size is " << maxSketchSize << "\n";
        }
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
    template<class LenT>
    void deserialize(std::istream& is, LenT nkeys, LenT nlocations)
    {
        using len_t = LenT;

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
    numGPUs_(0),
    maxLoadFactor_(default_max_load_factor()),
    maxLocationsPerFeature_(max_supported_locations_per_feature()),
    valid_(true)
{
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount); CUERR
    numGPUs_ = deviceCount;

    std::cerr << "found " << numGPUs_ << " CUDA devices\n";
}

//-----------------------------------------------------
template<class Key, class ValueT>
gpu_hashmap<Key,ValueT>::~gpu_hashmap() = default;

//-----------------------------------------------------
template<class Key, class ValueT>
gpu_hashmap<Key,ValueT>::gpu_hashmap(gpu_hashmap&& other) :
    numGPUs_{other.numGPUs_},
    maxLoadFactor_{other.maxLoadFactor_},
    maxLocationsPerFeature_{other.maxLocationsPerFeature_},
    valid_{other.valid_.exchange(false)},
    buildHashTables_{std::move(other.buildHashTables_)},
    queryHashTables_{std::move(other.queryHashTables_)}
{};



//---------------------------------------------------------------
template<class Key, class ValueT>
void gpu_hashmap<Key,ValueT>::pop_status(part_id gpuId) {
    if(gpuId < buildHashTables_.size()) {
        cudaSetDevice(gpuId); CUERR
        std::cerr
            << "gpu " << gpuId
            << " hashtable status: " << buildHashTables_[gpuId].pop_status() << "\n";
    }
    else if(gpuId < queryHashTables_.size()) {
        cudaSetDevice(gpuId); CUERR
        std::cerr
            << "gpu " << gpuId
            << " hashtable status: " << queryHashTables_[gpuId].pop_status() << "\n";
    }
}
//---------------------------------------------------------------
template<class Key, class ValueT>
void gpu_hashmap<Key,ValueT>::pop_status() {
    for(part_id gpuId = 0; gpuId < buildHashTables_.size(); ++gpuId) {
        cudaSetDevice(gpuId); CUERR
        std::cerr
            << "gpu " << gpuId
            <<  " hashtable status: " << buildHashTables_[gpuId].pop_status() << "\n";
    }
    for(part_id gpuId = 0; gpuId < queryHashTables_.size(); ++gpuId) {
        cudaSetDevice(gpuId); CUERR
        std::cerr
            << "gpu " << gpuId
            <<  " hashtable status: " << queryHashTables_[gpuId].pop_status() << "\n";
    }
}

//---------------------------------------------------------------
template<class Key, class ValueT>
size_t gpu_hashmap<Key,ValueT>::bucket_count() const noexcept {
    size_t count = 0;
    for(part_id gpuId = 0; gpuId < buildHashTables_.size(); ++gpuId) {
        cudaSetDevice(gpuId); CUERR
        count += buildHashTables_[gpuId].bucket_count();
    }
    for(part_id gpuId = 0; gpuId < queryHashTables_.size(); ++gpuId) {
        cudaSetDevice(gpuId); CUERR
        count += queryHashTables_[gpuId].bucket_count();
    }
    return count;
}

//---------------------------------------------------------------
template<class Key, class ValueT>
size_t gpu_hashmap<Key,ValueT>::key_count() noexcept {
    size_t count = 0;
    for(part_id gpuId = 0; gpuId < buildHashTables_.size(); ++gpuId) {
        cudaSetDevice(gpuId); CUERR
        count += buildHashTables_[gpuId].key_count();
    }
    for(part_id gpuId = 0; gpuId < queryHashTables_.size(); ++gpuId) {
        cudaSetDevice(gpuId); CUERR
        count += queryHashTables_[gpuId].key_count();
    }
    return count;
}

//-----------------------------------------------------
template<class Key, class ValueT>
size_t gpu_hashmap<Key,ValueT>::value_count() noexcept {
    size_t count = 0;
    for(part_id gpuId = 0; gpuId < buildHashTables_.size(); ++gpuId) {
        cudaSetDevice(gpuId); CUERR
        count += buildHashTables_[gpuId].location_count();
    }
    for(part_id gpuId = 0; gpuId < queryHashTables_.size(); ++gpuId) {
        cudaSetDevice(gpuId); CUERR
        count += queryHashTables_[gpuId].location_count();
    }
    return count;
}

//---------------------------------------------------------------
template<class Key, class ValueT>
statistics_accumulator_gpu<policy::Host>
gpu_hashmap<Key,ValueT>::location_list_size_statistics() {
    statistics_accumulator_gpu<policy::Host> totalAccumulator = {};

    for(part_id gpuId = 0; gpuId < numGPUs_; ++gpuId) {
        cudaSetDevice(gpuId); CUERR
        auto accumulator = buildHashTables_[gpuId].location_list_size_statistics();

        std::cout
            << "------------------------------------------------\n"
            << "gpu " << gpuId << ":\n"
            << "hashtable status     " << buildHashTables_[gpuId].pop_status() << '\n'
            << "buckets              " << buildHashTables_[gpuId].bucket_count() << '\n'
            << "bucket size          " << "max: " << accumulator.max()
                                       << " mean: " << accumulator.mean()
                                       << " +/- " << accumulator.stddev()
                                       << " <> " << accumulator.skewness() << '\n'
            << "features             " << std::uint64_t(accumulator.size()) << '\n'
            << "dead features        " << dead_feature_count() << '\n'
            << "locations            " << std::uint64_t(accumulator.sum()) << '\n';
            // << "load                 " << buildHashTables_[gpuId].load_factor() << '\n';

        totalAccumulator += accumulator;
    }

    return totalAccumulator;
}



//---------------------------------------------------------------
template<class Key, class ValueT>
part_id gpu_hashmap<Key,ValueT>::initialize_build_hash_tables(part_id numGPUs)
{
    if(numGPUs < numGPUs_)
        numGPUs_ = numGPUs;

    std::cerr << "using " << numGPUs_ << " CUDA devices\n";

    insertBuffers_.reserve(numGPUs_);

    for(part_id gpuId = 0; gpuId < numGPUs_; ++gpuId) {
        cudaSetDevice(gpuId); CUERR

        size_t freeMemory = 0;
        size_t totalMemory = 0;
        cudaMemGetInfo(&freeMemory, &totalMemory); CUERR
        std::cerr << "gpu " << gpuId << " freeMemory: " << helpers::B2GB(freeMemory) << " GB\n";

        // keep 1 GB of memory free aside from hash table
        const size_t tableMemory = freeMemory - (1ULL << 30);

        constexpr size_t valueSize = sizeof(ValueT);

        const size_t keyCapacity   = tableMemory *  2/13 / (2*valueSize);
        const size_t valueCapacity = tableMemory * 11/13 / valueSize;

        std::cerr << "gpu " << gpuId
                  << " allocate hashtable for " << keyCapacity << " keys"
                                       " and " << valueCapacity << " values\n";
        buildHashTables_.emplace_back(keyCapacity, valueCapacity, maxLocationsPerFeature_);

        cudaMemGetInfo(&freeMemory, &totalMemory); CUERR
        std::cerr << "gpu " << gpuId << " freeMemory: " << helpers::B2GB(freeMemory) << " GB\n";

        // allocate host buffers
        insertBuffers_.emplace_back();
    }

    return numGPUs_;
}


//---------------------------------------------------------------
template<class Key, class ValueT>
window_id gpu_hashmap<Key,ValueT>::add_target(
    part_id gpuId, const sequence& seq, target_id tgt, const sketcher& targetSketcher)
{
    using std::begin;
    using std::end;

    return add_target(gpuId, begin(seq), end(seq), tgt, targetSketcher);
}
//-----------------------------------------------------
template<class Key, class ValueT>
window_id gpu_hashmap<Key,ValueT>::add_target(
    part_id gpuId,
    sequence::const_iterator first,
    sequence::const_iterator last,
    target_id tgt,
    const sketcher& targetSketcher)
{
    // std::cerr << "add target " << tgt << " to gpu " << gpuId << "\n";

    using std::distance;

    window_id totalWindows = 0;

    for(window_id processedWindows = 0;
        distance(first, last) >= targetSketcher.kmer_size();
        first += processedWindows*targetSketcher.window_stride())
    {
        //fill sequence batch
        processedWindows = insertBuffers_[gpuId].current_seq_batch().add_target(
            first, last, tgt, totalWindows, targetSketcher);

        // if no windows were processed batch must be full
        if(!processedWindows && insertBuffers_[gpuId].current_seq_batch().num_targets()) {
            // std::cerr << "gpu " << gpuId << " insert\n";
            insert(gpuId, insertBuffers_[gpuId].current_seq_batch(), targetSketcher);
            insertBuffers_[gpuId].switch_seq_batch();
            insertBuffers_[gpuId].current_seq_batch().clear();
        }

        totalWindows += processedWindows;
    }

    return totalWindows;
}


//---------------------------------------------------------------
template<class Key, class ValueT>
void gpu_hashmap<Key,ValueT>::insert(
    part_id gpuId,
    sequence_batch<policy::Host>& seqBatchHost,
    const sketcher& targetSketcher)
{
    cudaSetDevice(gpuId); CUERR

    if(valid_ && buildHashTables_[gpuId].validate()) {
        buildHashTables_[gpuId].insert_async(
            seqBatchHost, targetSketcher);
    }
    else {
        valid_ = false;
    }
}
//-----------------------------------------------------
template<class Key, class ValueT>
void gpu_hashmap<Key,ValueT>::wait_until_add_target_complete(
    part_id gpuId, const sketcher& targetSketcher)
{
    if(gpuId < numGPUs_) {
        cudaSetDevice(gpuId); CUERR

        if(insertBuffers_[gpuId].current_seq_batch().num_targets()) {
            insert(gpuId, insertBuffers_[gpuId].current_seq_batch(), targetSketcher);
        }

        buildHashTables_[gpuId].wait_until_insert_finished();
    }
}


//---------------------------------------------------------------
template<class Key, class ValueT>
void gpu_hashmap<Key,ValueT>::query_async(
    query_batch<value_type>& batch,
    part_id hostId,
    const sketcher& querySketcher,
    bool copyAllHits,
    taxon_rank lowestRank) const
{
    for(part_id gpuId = 0; gpuId < batch.num_gpus(); ++gpuId)
    {
        cudaSetDevice(gpuId); CUERR

        batch.wait_for_allhits_copied(gpuId);

        if(gpuId == 0) {
            batch.copy_queries_to_device_async(hostId);

            queryHashTables_[gpuId].query_sequences_async(
                batch.host_data(hostId).num_windows(),
                batch.gpu_data(gpuId),
                querySketcher,
                maxLocationsPerFeature_);
        }
        else {
            queryHashTables_[gpuId].query_sketches_async(
                batch.host_data(hostId).num_windows(),
                batch.gpu_data(gpuId),
                querySketcher,
                maxLocationsPerFeature_);
        }
        batch.mark_query_finished(gpuId);

        // batch.sync_work_stream(gpuId); CUERR

        if(gpuId < batch.num_gpus()-1)
            batch.copy_queries_to_next_device_async(hostId, gpuId);

        batch.compact_sort_and_copy_allhits_async(hostId, gpuId, copyAllHits);

        batch.generate_and_copy_top_candidates_async(
            hostId, gpuId, queryHashTables_[gpuId].lineages(), lowestRank);

        // batch.sync_copy_stream(gpuId); CUERR
    }
}


//---------------------------------------------------------------
template<class Key, class ValueT>
void gpu_hashmap<Key,ValueT>::deserialize(std::istream& is, part_id gpuId)
{
    using len_t = std::uint64_t;

    len_t nkeys = 0;
    read_binary(is, nkeys);
    len_t nvalues = 0;
    read_binary(is, nvalues);

    std::cerr << "\n\t#features: " << nkeys << " #locations: " << nvalues << "\n";

    if(nkeys > 0) {
        cudaSetDevice(gpuId); CUERR
        std::cerr << "\tloading database to gpu " << gpuId << "\n";

        size_t freeMemory = 0;
        size_t totalMemory = 0;
        cudaMemGetInfo(&freeMemory, &totalMemory); CUERR
        std::cerr << "\tfreeMemory: " << helpers::B2GB(freeMemory) << " GB\n";

        //initialize hash table
        queryHashTables_.emplace_back(nkeys/maxLoadFactor_);

        std::cerr << "\tfeatures capacity: " << queryHashTables_.back().bucket_count() << "\n";

        size_t indexSize = queryHashTables_.back().bucket_count() * (sizeof(value_type) + sizeof(value_type));
        std::cerr << "\tindex size: " << helpers::B2GB(indexSize) << " GB\n";

        // load hash table
        queryHashTables_.back().deserialize(is, nkeys, nvalues);

        size_t valuesSize = nvalues*sizeof(location);
        std::cerr << "\tlocations size: " << helpers::B2GB(valuesSize) << " GB\n";
        std::cerr << "\ttotal size: " << helpers::B2GB(indexSize + valuesSize) << " GB\n";

        cudaMemGetInfo(&freeMemory, &totalMemory); CUERR
        std::cerr << "\tfreeMemory: " << helpers::B2GB(freeMemory) << " GB\n";
    }
}


//---------------------------------------------------------------
/**
* @brief binary serialization of all non-empty buckets
*/
template<class Key, class ValueT>
void gpu_hashmap<Key,ValueT>::serialize(std::ostream& os, part_id gpuId)
{
    cudaSetDevice(gpuId); CUERR
    buildHashTables_[gpuId].serialize(os);
}


//---------------------------------------------------------------
template<class Key, class ValueT>
void gpu_hashmap<Key,ValueT>::copy_target_lineages_to_gpu(
    const std::vector<ranked_lineage>& lins,
    part_id gpuId)
{
    cudaSetDevice(gpuId); CUERR
    queryHashTables_[gpuId].copy_target_lineages_to_gpu(lins);
}


//---------------------------------------------------------------
template<class Key, class ValueT>
void gpu_hashmap<Key,ValueT>::enable_all_peer_access(part_id numGPUs)
{
    for (part_id srcId = 0; srcId < numGPUs; ++srcId) {
        cudaSetDevice(srcId);
        for (part_id dstId = 0; dstId < numGPUs; ++dstId) {
            if (srcId != dstId) {
                 cudaDeviceEnablePeerAccess(dstId, 0);
            }
        }
    }
}



//---------------------------------------------------------------
template class gpu_hashmap<kmer_type, location>;

} // namespace mc

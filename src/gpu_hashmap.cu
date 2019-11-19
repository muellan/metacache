#include <limits>

#include "gpu_hashmap.cuh"
#include "hash_dna.h"
#include "hash_int.h"
#include "sketch_database.h"
#include "gpu_engine.cuh"

#include "../dep/warpcore/include/warpcore.cuh"

namespace mc {


//---------------------------------------------------------------
template<
    class Key,
    class ValueT,
    class Hash,
    class KeyEqual,
    class BucketSizeT
>
gpu_hashmap<Key,ValueT,Hash,KeyEqual,BucketSizeT>::feature_batch::feature_batch(feature_batch::counter_type maxFeatures) :
    maxFeatures_{maxFeatures}
{
    if(maxFeatures_) {
        cudaMalloc(&features_, maxFeatures_*sizeof(Key));
        cudaMalloc(&values_, maxFeatures_*sizeof(ValueT));
        cudaMalloc(&featureCounter_, sizeof(feature_batch::counter_type));
    }
    CUERR
}
//---------------------------------------------------------------
template<
    class Key,
    class ValueT,
    class Hash,
    class KeyEqual,
    class BucketSizeT
>
gpu_hashmap<Key,ValueT,Hash,KeyEqual,BucketSizeT>::feature_batch::~feature_batch() {
    if(maxFeatures_) {
        cudaFree(features_);
        cudaFree(values_);
        cudaFree(featureCounter_);
    }
    CUERR
}


/*************************************************************************//**
 * TODO
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
template<
    class Key,
    class ValueT,
    class Hash,
    class KeyEqual,
    class BucketSizeT
>
class gpu_hashmap<Key,ValueT,Hash,KeyEqual,BucketSizeT>::hash_table {

    using key_type   = Key;
    using value_type = std::uint64_t;
    using location_type = ValueT;
    using size_type  = size_t;

    // using hasher_t = warpcore::hashers::MuellerHash;
    using hasher_t = warpcore::defaults::hasher_t<Key>;

    // using probing_t = warpcore::probing_schemes::DoubleHashing<hasher_t, hasher_t>;
    using probing_t = warpcore::probing_schemes::QuadraticProbing<hasher_t>;

    using hash_table_t = warpcore::SingleValueHashTable<
        key_type, value_type,
        // warpcore::defaults::probing_t<key_type>,
        probing_t,
        // warpcore::defaults::empty_key<key_type>(),          //=0
        key_type(-2),
        warpcore::defaults::tombstone_key<key_type>()>;     //=-1

public:
    hash_table(size_t capacity) :
        hashTable_(capacity),
        numKeys_(0), numLocations_(0),
        locations_(nullptr)
    {
        // std::cerr << "capacity: " << hashTable_.capacity() << std::endl;
    }

    //---------------------------------------------------------------
    float load_factor() const noexcept {
        return hashTable_.load_factor();
    }
    //---------------------------------------------------------------
    size_type key_count() const noexcept {
        return numKeys_;
    }
    //-----------------------------------------------------
    size_type location_count() const noexcept {
        return numLocations_;
    }

    //---------------------------------------------------------------
    template<class result_type>
    void query(query_batch<result_type>& batch,
               const sketcher& querySketcher,
               bucket_size_type maxLocationPerFeature) const
    {
        const cudaStream_t stream = batch.stream();
        batch.copy_queries_to_device_async();

        // max 32*4 features => max window size is 128
        #define BLOCK_THREADS 32
        #define ITEMS_PER_THREAD 4

        const size_type probingLength = 4*hashTable_.capacity();

        const size_t numBlocks = batch.num_queries();
        gpu_hahstable_query<BLOCK_THREADS,ITEMS_PER_THREAD><<<numBlocks,BLOCK_THREADS,0,stream>>>(
            hashTable_,
            probingLength,
            batch.num_queries(),
            batch.query_ids_device(),
            batch.encode_offsets_device(),
            batch.encoded_seq_device(),
            batch.encoded_ambig_device(),
            querySketcher.kmer_size(),
            querySketcher.sketch_size(),
            querySketcher.window_size(),
            querySketcher.window_stride(),
            locations_,
            maxLocationPerFeature,
            batch.query_results_device(),
            batch.result_counts_device(),
            batch.result_offsets_device()
        );

        batch.compact_sort_and_copy_results_async();

        //TODO async
        batch.sync_stream();
        CUERR
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
                    d_keyBuffer, d_offsetBuffer, keyBatchSize, probingLength, stream, status);
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
                    d_keyBuffer, d_offsetBuffer, remainingSize, probingLength, stream, status);

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

        numKeys_ = nkeys;
        numLocations_ = nlocations;
    }

private:
    hash_table_t hashTable_;

    size_type numKeys_;
    size_type numLocations_;
    location * locations_;
};


//---------------------------------------------------------------
template<
    class Key,
    class ValueT,
    class Hash,
    class KeyEqual,
    class BucketSizeT
>
gpu_hashmap<Key,ValueT,Hash,KeyEqual,BucketSizeT>::gpu_hashmap() :
    maxLoadFactor_(default_max_load_factor()),
    hash_{}, keyEqual_{},
    seqBatches_{}, featureBatches_{}
{
    cudaSetDevice(0);

    //FIXME unsafe
    //possible fix: limit number of windows in sequence batch
    //shortest possible window size to create full sketch can be encoded in 2 blocks
    // const size_t winLen = 2;
    //average should be bigger than this
    const size_t winLen = 5;

    const size_t maxWindows = MAX_ENCODE_LENGTH_PER_BATCH / winLen;

    //TODO get sketch size from sketcher
    sketch_size_type sketchSize = 16;

    const size_t maxFeatures = maxWindows * sketchSize;
    std::cout << "max features: " << maxFeatures << '\n';

    seqBatches_.emplace_back(MAX_TARGETS_PER_BATCH, MAX_ENCODE_LENGTH_PER_BATCH);

    featureBatches_.emplace_back(maxFeatures);
}

//-----------------------------------------------------
template<
    class Key,
    class ValueT,
    class Hash,
    class KeyEqual,
    class BucketSizeT
>
gpu_hashmap<Key,ValueT,Hash,KeyEqual,BucketSizeT>::~gpu_hashmap() = default;

//-----------------------------------------------------
template<
    class Key,
    class ValueT,
    class Hash,
    class KeyEqual,
    class BucketSizeT
>
gpu_hashmap<Key,ValueT,Hash,KeyEqual,BucketSizeT>::gpu_hashmap(gpu_hashmap&&) = default;



//---------------------------------------------------------------
template<
class Key,
class ValueT,
class Hash,
class KeyEqual,
class BucketSizeT
>
size_t gpu_hashmap<Key,ValueT,Hash,KeyEqual,BucketSizeT>::key_count() const noexcept {
    return hashTable_ ? hashTable_->key_count() : 0;
}

//-----------------------------------------------------
template<
class Key,
class ValueT,
class Hash,
class KeyEqual,
class BucketSizeT
>
size_t gpu_hashmap<Key,ValueT,Hash,KeyEqual,BucketSizeT>::value_count() const noexcept {
    return hashTable_ ? hashTable_->location_count() : 0;
}

//---------------------------------------------------------------
template<
class Key,
class ValueT,
class Hash,
class KeyEqual,
class BucketSizeT
>
float gpu_hashmap<Key,ValueT,Hash,KeyEqual,BucketSizeT>::load_factor() const noexcept {
    return hashTable_ ? hashTable_->load_factor() : -1;
}



//---------------------------------------------------------------
template<
    class Key,
    class ValueT,
    class Hash,
    class KeyEqual,
    class BucketSizeT
>
std::vector<Key> gpu_hashmap<Key,ValueT,Hash,KeyEqual,BucketSizeT>::insert(
    const sequence_batch<policy::Host>& seqBatchHost,
    const sketcher& targetSketcher
) {
    using counter_type = typename feature_batch::counter_type;

    cudaStream_t stream = 0;

    seqBatches_[0].num_targets(seqBatchHost.num_targets());

    //copy batch to gpu
    cudaMemcpyAsync(seqBatches_[0].target_ids(), seqBatchHost.target_ids(),
                    seqBatchHost.num_targets()*sizeof(target_id),
                    cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(seqBatches_[0].window_offsets(), seqBatchHost.window_offsets(),
                    seqBatchHost.num_targets()*sizeof(window_id),
                    cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(seqBatches_[0].encode_offsets(), seqBatchHost.encode_offsets(),
                    (seqBatchHost.num_targets()+1)*sizeof(encodinglen_t),
                    cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(seqBatches_[0].encoded_seq(), seqBatchHost.encoded_seq(),
                    seqBatchHost.encode_offsets()[seqBatchHost.num_targets()]*sizeof(encodedseq_t),
                    cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(seqBatches_[0].encoded_ambig(), seqBatchHost.encoded_ambig(),
                    seqBatchHost.encode_offsets()[seqBatchHost.num_targets()]*sizeof(encodedambig_t),
                    cudaMemcpyHostToDevice, stream);

    //initialize counter
    cudaMemsetAsync(featureBatches_[0].feature_counter(), 0, sizeof(counter_type), stream);
    // cudaStreamSynchronize(stream);
    // CUERR

    // max 32*4 features => max window size is 128
    #define BLOCK_THREADS 32
    #define ITEMS_PER_THREAD 4

    //TODO increase grid in x and y dim
    constexpr dim3 numBlocks{1,1};
    extract_features<BLOCK_THREADS,ITEMS_PER_THREAD><<<numBlocks,BLOCK_THREADS,0,stream>>>(
        seqBatches_[0].num_targets(),
        seqBatches_[0].target_ids(),
        seqBatches_[0].window_offsets(),
        seqBatches_[0].encode_offsets(),
        seqBatches_[0].encoded_seq(),
        seqBatches_[0].encoded_ambig(),
        targetSketcher.kmer_size(),
        targetSketcher.sketch_size(),
        targetSketcher.window_size(),
        targetSketcher.window_stride(),
        featureBatches_[0].features(),
        featureBatches_[0].values(),
        featureBatches_[0].feature_counter());
    //TODO remove later ------------------------------------------------------------
    //neccessary sync to get feature counter
    cudaStreamSynchronize(stream);
    CUERR

    counter_type h_featureCounter = 0;
    cudaMemcpyAsync(&h_featureCounter, featureBatches_[0].feature_counter(),
                sizeof(counter_type), cudaMemcpyDeviceToHost, stream);
    //neccessary sync to get feature counter
    cudaStreamSynchronize(stream);
    std::cout << "Counter: " << h_featureCounter << std::endl;

    // kmer_type * h_features;
    // cudaMallocHost(&h_features, numFeatures*sizeof(Key));
    std::vector<Key> h_features(h_featureCounter);
    std::vector<ValueT> h_values(h_featureCounter);
    cudaMemcpyAsync(h_features.data(), featureBatches_[0].features(),
                h_featureCounter*sizeof(Key), cudaMemcpyDeviceToHost, stream);
    cudaMemcpyAsync(h_values.data(), featureBatches_[0].values(),
                h_featureCounter*sizeof(ValueT), cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream);
    CUERR

    //print features
    // for(size_t i=0; i<h_featureCounter; ++i) {
    //     std:: cout << h_features[i] << ' ';
    //     if((i+1) % sketchSize_ == 0)
    //         std::cout << '\n';
    // }
    // std::cout << std::endl;

    // for(size_t i=0; i<h_featureCounter; ++i) {
    //     std:: cout << '(' << h_values[i].tgt << ',' << h_values[i].win << ") " ;
    //     if((i+1) % sketchSize_ == 0)
    //         std::cout << '\n';
    // }
    // std::cout << std::endl;

    //TODO remove later ------------------------------------------------------------

    return h_features;
}


//---------------------------------------------------------------
template<
    class Key,
    class ValueT,
    class Hash,
    class KeyEqual,
    class BucketSizeT
>
void gpu_hashmap<Key,ValueT,Hash,KeyEqual,BucketSizeT>::query(query_batch<value_type>& batch, const sketcher& querySketcher, bucket_size_type maxLocationPerFeature) const {
    hashTable_->query(batch, querySketcher, maxLocationPerFeature);
}


//---------------------------------------------------------------
template<
    class Key,
    class ValueT,
    class Hash,
    class KeyEqual,
    class BucketSizeT
>
void gpu_hashmap<Key,ValueT,Hash,KeyEqual,BucketSizeT>::deserialize(
    std::istream& is)
{
    using len_t = std::uint64_t;

    len_t nkeys = 0;
    read_binary(is, nkeys);
    len_t nvalues = 0;
    read_binary(is, nvalues);

    std::cout << "nkeys: " << nkeys << " nvalues: " << nvalues << std::endl;

    if(nkeys > 0) {
        //initialize hash table
        hashTable_ = std::make_unique<hash_table>(nkeys/maxLoadFactor_);

        hashTable_->deserialize(is, nkeys, nvalues);
    }
}


//---------------------------------------------------------------
/**
* @brief binary serialization of all non-emtpy buckets
*/
template<
    class Key,
    class ValueT,
    class Hash,
    class KeyEqual,
    class BucketSizeT
>
void gpu_hashmap<Key,ValueT,Hash,KeyEqual,BucketSizeT>::serialize(
    std::ostream& os) const
{
    //TODO
}



//---------------------------------------------------------------
template class gpu_hashmap<
        kmer_type,
        location,
        // uint64_t,
        same_size_hash<kmer_type>,
        std::equal_to<kmer_type>,
        unsigned char
        >;

} // namespace mc

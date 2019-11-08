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

    using hash_table_t = warpcore::SingleValueHashTable<
        key_type, value_type,
        warpcore::defaults::probing_t<key_type>,
        // warpcore::defaults::empty_key<key_type>(),          //=0
        key_type(-2),
        warpcore::defaults::tombstone_key<key_type>()>;     //=-1

public:
    hash_table(size_t capacity) :
        hashTable_(capacity),
        numKeys_(0), numLocations_(0),
        locations_(nullptr)
    {}

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
    void insert(key_type *keys_in, value_type *values_in, size_type size_in)
    {
        size_type probingLength = 4*hashTable_.capacity();
        //TODO
        cudaStream_t stream = 0;
        hashTable_.insert(keys_in, values_in, size_in, probingLength, stream);
    }

    //---------------------------------------------------------------
    template<class result_type>
    void query(query_batch<result_type>& batch,
               const sketcher& querySketcher,
               bucket_size_type maxLocationPerFeature) const
    {
        //TODO
        cudaStream_t stream = batch.stream();
        batch.copy_queries_to_device_async();

        // max 32*4 features => max window size is 128
        #define BLOCK_THREADS 32
        #define ITEMS_PER_THREAD 4

        size_type probingLength = 4*hashTable_.capacity();
        //TODO increase grid size
        constexpr size_t numBlocks = 32;
        gpu_hahstable_query<BLOCK_THREADS,ITEMS_PER_THREAD><<<numBlocks,BLOCK_THREADS,0,stream>>>(
            hashTable_,
            probingLength,
            batch.num_queries(),
            batch.encode_offsets_device(),
            batch.encoded_seq_device(),
            batch.encoded_ambig_device(),
            querySketcher.kmer_size(),
            querySketcher.sketch_size(),
            querySketcher.window_size(),
            querySketcher.window_stride(),
            locations_,
            maxLocationPerFeature,
            batch.query_results_device());

        batch.copy_results_to_host_tmp_async();

        batch.sort_results();

        batch.copy_results_to_host_async();
        //TODO async
        batch.sync_stream();
    }

    //---------------------------------------------------------------
    template<class len_t>
    void deserialize(std::istream& is, len_t nkeys, len_t nlocations)
    {
        //TODO tune sizes
        const size_t keyBatchSize = 1UL << 20;
        const size_t valBatchSize = 1UL << 20;

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

            using handler_type = warpcore::status_handlers::ReturnStatus;
            using handler_base_type = handler_type::base_type;

            handler_base_type * status;
            cudaMallocManaged(&status, keyBatchSize*sizeof(uint64_t));
            cudaMemset(status, 1, keyBatchSize*sizeof(uint64_t));

            uint64_t locsOffset = 0;
            for(len_t i = 0; i < nkeys; ++i) {
                key_type key;
                bucket_size_type nlocs = 0;
                read_binary(is, key);
                read_binary(is, nlocs);

                h_keyBuffer[i % keyBatchSize] = key;
                //store offset and size together in 64bit
                //default is 56bit offset, 8bit size
                h_offsetBuffer[i % keyBatchSize] =
                    (locsOffset << sizeof(bucket_size_type)*CHAR_BIT) + nlocs;

                locsOffset += nlocs;

                //insert buffer if full
                //TODO overlap async
                if((i+1) % keyBatchSize == 0) {
                    cudaMemcpy(d_keyBuffer, h_keyBuffer, keyBatchSize*sizeof(key_type),
                                cudaMemcpyHostToDevice);
                    cudaMemcpy(d_offsetBuffer, h_offsetBuffer, keyBatchSize*sizeof(uint64_t),
                                cudaMemcpyHostToDevice);
                    // insert(d_keyBuffer, d_offsetBuffer, keyBatchSize);
                    size_type probingLength = hashTable_.capacity();
                    hashTable_.template insert<handler_type>(d_keyBuffer, d_offsetBuffer, keyBatchSize, probingLength, 0, status);
                    cudaDeviceSynchronize();

                    for(size_t j=0; j<keyBatchSize; ++j) {
                        if(status[j % keyBatchSize].has_any()) {
                            std::cout << h_keyBuffer[j % keyBatchSize] << ' ' << status[j % keyBatchSize] << std::endl;
                        }
                    }
                }
            }
            const size_t remainingSize = nkeys % keyBatchSize;
            //insert remaining pairs in buffer
            cudaMemcpy(d_keyBuffer, h_keyBuffer, remainingSize*sizeof(key_type),
                        cudaMemcpyHostToDevice);
            cudaMemcpy(d_offsetBuffer, h_offsetBuffer, remainingSize*sizeof(uint64_t),
                        cudaMemcpyHostToDevice);
            // insert(d_keyBuffer, d_offsetBuffer, remainingSize);
            size_type probingLength = hashTable_.capacity();
            hashTable_.template insert<handler_type>(d_keyBuffer, d_offsetBuffer, remainingSize, probingLength, 0, status);
            cudaDeviceSynchronize();

            for(size_t j=0; j<remainingSize; ++j) {
                if(status[j % keyBatchSize].has_any()) {
                    std::cout << h_keyBuffer[j % keyBatchSize] << ' ' << status[j % keyBatchSize] << std::endl;
                }
            }
            // std::cout << hashTable_.peek_status() << std::endl;

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
            location * valueBuffer;
            cudaMallocHost(&valueBuffer, valBatchSize*sizeof(location));

            //read batch of locations and copy to device
            //TODO overlap async
            auto locsOffset = locations_;
            for(size_t i = 0; i < nlocations/valBatchSize; ++i) {
                read_binary(is, valueBuffer, valBatchSize);
                cudaMemcpy(locsOffset, valueBuffer, valBatchSize*sizeof(location),
                            cudaMemcpyHostToDevice);
                locsOffset += valBatchSize;
            }
            //read remaining locations and copy to device
            const size_t remainingSize = nlocations % valBatchSize;
            read_binary(is, valueBuffer, remainingSize);
            cudaMemcpy(locsOffset, valueBuffer, remainingSize*sizeof(location),
                        cudaMemcpyHostToDevice);

            cudaFreeHost(valueBuffer);
        }

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

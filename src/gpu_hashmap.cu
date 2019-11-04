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


//---------------------------------------------------------------
template<
    class Key,
    class ValueT,
    class Hash,
    class KeyEqual,
    class BucketSizeT
>
class gpu_hashmap<Key,ValueT,Hash,KeyEqual,BucketSizeT>::single_value_hash_table {

    using key_type   = Key;
    using value_type = std::uint64_t;
    using size_type  = size_t;

    using hash_table_t = warpcore::SingleValueHashTable<
        key_type, value_type,
        warpcore::defaults::probing_t<key_type>,
        warpcore::defaults::empty_key<key_type>(),          //=0
        warpcore::defaults::tombstone_key<key_type>()>;     //=-1

public:
    single_value_hash_table(size_t capacity) : hashTable_(capacity) {}

    //---------------------------------------------------------------
    float load_factor() const noexcept {
        return hashTable_.load_factor();
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
    void query(query_batch<result_type>& batch, const sketcher& querySketcher) const {
        //TODO
        cudaStream_t stream = 0;
        batch.copy_queries_to_device(stream);

        // max 32*4 features => max window size is 128
        #define BLOCK_THREADS 32
        #define ITEMS_PER_THREAD 4

        size_type probingLength = 4*hashTable_.capacity();
        //TODO increase grid in x and y dim
        constexpr dim3 numBlocks{1,1};
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
            batch.query_results_device());

        batch.copy_results_to_host(stream);
    }

private:
    hash_table_t hashTable_;
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
    numKeys_(0), numValues_(0), maxLoadFactor_(default_max_load_factor()),
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
void gpu_hashmap<Key,ValueT,Hash,KeyEqual,BucketSizeT>::query(query_batch<value_type>& batch, const sketcher& querySketcher) const {
    hashTable_->query(batch, querySketcher);
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

    //TODO tune sizes
    const size_t keyBatchSize = 1UL << 20;
    const size_t valBatchSize = 1UL << 20;

    len_t nkeys = 0;
    read_binary(is, nkeys);
    len_t nvalues = 0;
    read_binary(is, nvalues);

    if(nkeys > 0) {
        {//load hash table
            //initialize hash table
            hashTable_ = std::make_unique<single_value_hash_table>(nkeys/maxLoadFactor_);

            //allocate insert buffers
            key_type * h_keyBuffer;
            key_type * d_keyBuffer;
            cudaMallocHost(&h_keyBuffer, keyBatchSize*sizeof(key_type));
            cudaMalloc    (&d_keyBuffer, keyBatchSize*sizeof(key_type));
            uint64_t * h_offsetBuffer;
            uint64_t * d_offsetBuffer;
            cudaMallocHost(&h_offsetBuffer, keyBatchSize*sizeof(uint64_t));
            cudaMalloc    (&d_offsetBuffer, keyBatchSize*sizeof(uint64_t));

            uint64_t valuesOffset = 0;
            for(len_t i = 0; i < nkeys; ++i) {
                key_type key;
                bucket_size_type nvals = 0;
                read_binary(is, key);
                read_binary(is, nvals);

                h_keyBuffer[i % keyBatchSize] = key;
                //store offset and size together in 64bit
                //default is 56bit offset, 8bit size
                h_offsetBuffer[i % keyBatchSize] =
                    (valuesOffset << sizeof(bucket_size_type)*CHAR_BIT) + nvals;;

                valuesOffset += nvals;

                //insert buffer if full
                //TODO overlap async
                if(i % keyBatchSize == 0) {
                    cudaMemcpy(d_keyBuffer, h_keyBuffer, keyBatchSize*sizeof(key_type),
                               cudaMemcpyHostToDevice);
                    cudaMemcpy(d_offsetBuffer, h_offsetBuffer, keyBatchSize*sizeof(uint64_t),
                               cudaMemcpyHostToDevice);
                    hashTable_->insert(d_keyBuffer, d_offsetBuffer, keyBatchSize);
                }
            }
            //insert remaining pairs in buffer
            cudaMemcpy(d_keyBuffer, h_keyBuffer, keyBatchSize*sizeof(key_type),
                       cudaMemcpyHostToDevice);
            cudaMemcpy(d_offsetBuffer, h_offsetBuffer, keyBatchSize*sizeof(uint64_t),
                       cudaMemcpyHostToDevice);
            hashTable_->insert(d_keyBuffer, d_offsetBuffer, nkeys % keyBatchSize);

            cudaFreeHost(h_keyBuffer);
            cudaFree    (d_keyBuffer);
            cudaFreeHost(h_offsetBuffer);
            cudaFree    (d_offsetBuffer);
        }

        {//load values
            //allocate large memory chunk for all values,
            //individual buckets will then point into this array
            cudaMalloc(&values_, nvalues*sizeof(value_type));

            //allocate buffer
            value_type * valueBuffer;
            cudaMallocHost(&valueBuffer, valBatchSize*sizeof(value_type));

            //read batch of values and copy to device
            //TODO overlap async
            auto valuesOffset = values_;
            for(size_t i = 0; i < nvalues/valBatchSize; ++i) {
                read_binary(is, valueBuffer, valBatchSize);
                cudaMemcpy(valuesOffset, valueBuffer, valBatchSize*sizeof(value_type),
                           cudaMemcpyHostToDevice);
                valuesOffset += valBatchSize;
            }
            //read remaining values and copy to device
            const size_t remainingSize = nvalues % valBatchSize;
            read_binary(is, valueBuffer, remainingSize);
            cudaMemcpy(valuesOffset, valueBuffer, remainingSize*sizeof(value_type),
                       cudaMemcpyHostToDevice);

            cudaFreeHost(valueBuffer);
        }

        numKeys_ = nkeys;
        numValues_ = nvalues;
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

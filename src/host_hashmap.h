#ifndef MC_HOST_HASH_MAP_H_
#define MC_HOST_HASH_MAP_H_

#include "config.h"
#include "hash_multimap.h"
#include "stat_combined.h"
#include "batch_processing.h"
#include "taxonomy.h"


namespace mc {


template<
    class ValueT
>
class host_hashmap
{
public:
    //---------------------------------------------------------------
    using location         = ValueT;
    using match_locations  = std::vector<location>;
    using bucket_size_type = mc::loclist_size_t;

    using sketch  = typename sketcher::sketch_type;  //range of features
    using feature = typename sketch::value_type;

private:
    //-----------------------------------------------------
    /// @brief "heart of the database": maps features to target locations
    using hash_table = hash_multimap<feature,location, //key, value
                              feature_hash,               //key hasher
                              std::equal_to<feature>,     //key comparator
                              chunk_allocator<location>,  //value allocator
                              std::allocator<feature>,    //bucket+key allocator
                              bucket_size_type>;          //location list size

    //-----------------------------------------------------
    /// @brief needed for batched, asynchonous insertion into feature_store
    struct window_sketch
    {
        window_sketch() :
            tgt{}, win{}, sk{} {};

        window_sketch(target_id tgt, window_id win, sketch sk) :
            tgt{tgt}, win{win}, sk{std::move(sk)} {};

        target_id tgt;
        window_id win;
        sketch sk;
    };

    using sketch_batch = std::vector<window_sketch>;

public:
    //---------------------------------------------------------------
    using feature_count_type = typename hash_table::size_type;


    //---------------------------------------------------------------
    /** @brief used for query result storage/accumulation
     */
    class matches_sorter {
        friend class host_hashmap;

    public:
        void sort() {
            merge_sort(locs_, offsets_, temp_);
        }

        void clear() {
            locs_.clear();
            offsets_.clear();
            offsets_.resize(1, 0);
        }

        bool empty() const noexcept { return locs_.empty(); }
        auto size()  const noexcept { return locs_.size(); }

        auto begin() const noexcept { return locs_.begin(); }
        auto end()   const noexcept { return locs_.end(); }

        const match_locations&
        locations() const noexcept { return locs_; }

    private:
        static void
        merge_sort(match_locations& inout,
                   std::vector<size_t>& offsets,
                   match_locations& temp)
        {
            if(offsets.size() < 3) return;
            temp.resize(inout.size());

            int numChunks = offsets.size()-1;
            for(int s = 1; s < numChunks; s *= 2) {
                for(int i = 0; i < numChunks; i += 2*s) {
                    auto begin = offsets[i];
                    auto mid = i + s <= numChunks ? offsets[i + s] : offsets[numChunks];
                    auto end = i + 2*s <= numChunks ? offsets[i + 2*s] : offsets[numChunks];
                    std::merge(inout.begin()+begin, inout.begin()+mid,
                               inout.begin()+mid, inout.begin()+end,
                               temp.begin()+begin);
                }
                std::swap(inout, temp);
            }
        }

        match_locations locs_; // match locations from hashmap
        std::vector<std::size_t> offsets_;  // bucket sizes for merge sort
        match_locations temp_; // temp buffer for merge sort
    };

public:
    //---------------------------------------------------------------
    explicit
    host_hashmap() :
        hashTable_{},
        inserter_{}
    {
        hashTable_.max_load_factor(default_max_load_factor());
    }

    host_hashmap(const host_hashmap&) = delete;
    host_hashmap(host_hashmap&& other) :
        hashTable_{std::move(other.hashTable_)},
        inserter_{std::move(other.inserter_)}
    {}

    host_hashmap& operator = (const host_hashmap&) = delete;
    host_hashmap& operator = (host_hashmap&&)      = default;

    ~host_hashmap() {
        destroy_inserter();
    }


    //---------------------------------------------------------------
    unsigned num_gpus() const noexcept { return 0; }
    //---------------------------------------------------------------
    gpu_id initialize_build_hash_tables(gpu_id) { return 0; }


    //---------------------------------------------------------------
    static constexpr std::size_t
    max_bucket_size() noexcept {
        return hash_table::max_bucket_size();
    }


    //---------------------------------------------------------------
    void max_load_factor(float lf) {
        hashTable_.max_load_factor(lf);
    }
    //-----------------------------------------------------
    float max_load_factor() const noexcept {
        return hashTable_.max_load_factor();
    }
    //-----------------------------------------------------
    static constexpr float default_max_load_factor() noexcept {
        return 0.8f;
    }


    //---------------------------------------------------------------
    bool empty() const noexcept {
        return hashTable_.empty();
    }
    //---------------------------------------------------------------
    std::uint64_t key_count() const noexcept {
        return hashTable_.key_count();
    }
    //---------------------------------------------------------------
    std::uint64_t value_count() const noexcept {
        return hashTable_.value_count();
    }
    //---------------------------------------------------------------
    std::uint64_t bucket_count() const noexcept {
        return hashTable_.bucket_count();
    }
    //---------------------------------------------------------------
    std::uint64_t non_empty_bucket_count() const noexcept {
        return hashTable_.non_empty_bucket_count();
    }
    //---------------------------------------------------------------
    std::uint64_t dead_feature_count() const noexcept {
        return hashTable_.key_count() - hashTable_.non_empty_bucket_count();
    }

    //---------------------------------------------------------------
    void clear() {
        hashTable_.clear();
    }
    //---------------------------------------------------------------
    /**
     * @brief very dangerous! clears hash table without memory deallocation
     */
    void clear_without_deallocation() {
        hashTable_.clear_without_deallocation();
    }


    //---------------------------------------------------------------
    statistics_accumulator
    location_list_size_statistics() const {
        auto priSize = statistics_accumulator{};

        for(const auto& bucket : hashTable_) {
            if(!bucket.empty()) {
                priSize += bucket.size();
            }
        }

        return priSize;
    }


    //---------------------------------------------------------------
    void print_feature_map(std::ostream& os) const {
        for(const auto& bucket : hashTable_) {
            if(!bucket.empty()) {
                os << std::int_least64_t(bucket.key()) << " -> ";
                for(location p : bucket) {
                    os << '(' << std::int_least64_t(p.tgt)
                       << ',' << std::int_least64_t(p.win) << ')';
                }
                os << '\n';
            }
        }
    }


    //---------------------------------------------------------------
    void print_feature_counts(std::ostream& os) const {
        for(const auto& bucket : hashTable_) {
            if(!bucket.empty()) {
                os << std::int_least64_t(bucket.key()) << " -> "
                   << std::int_least64_t(bucket.size()) << '\n';
            }
        }
    }


    //---------------------------------------------------------------
    static bucket_size_type
    max_supported_locations_per_feature() noexcept {
        return (hash_table::max_bucket_size() - 1);
    }
    //-----------------------------------------------------
    void max_locations_per_feature(bucket_size_type n)
    {
        if(n < 1) n = 1;
        if(n >= max_supported_locations_per_feature()) {
            n = max_supported_locations_per_feature();
        }
        else if(n < maxLocsPerFeature_) {
            hashTable_.shrink_all(n);
        }
        maxLocsPerFeature_ = n;
    }
    //-----------------------------------------------------
    bucket_size_type
    max_locations_per_feature() const noexcept {
        return maxLocsPerFeature_;
    }


    // ----------------------------------------------------------------------------
    feature_count_type
    remove_features_with_more_locations_than(bucket_size_type n)
    {
        //note that features are not really removed, because the hashmap
        //does not support erasing keys; instead all values belonging to
        //the key are cleared and the key is kept without values
        feature_count_type rem = 0;
        for(auto i = hashTable_.begin(), e = hashTable_.end(); i != e; ++i) {
            if(i->size() > n) {
                hashTable_.clear(i);
                ++rem;
            }
        }
        return rem;
    }


    // ----------------------------------------------------------------------------
    feature_count_type
    remove_ambiguous_features(taxon_rank r, bucket_size_type maxambig,
                              ranked_lineages_of_targets& targetLineages)
    {
        feature_count_type rem = 0;

        if(maxambig == 0) maxambig = 1;

        if(r == taxon_rank::Sequence) {
            for(auto i = hashTable_.begin(), e = hashTable_.end(); i != e; ++i) {
                if(!i->empty()) {
                    std::set<target_id> targets;
                    for(auto loc : *i) {
                        targets.insert(loc.tgt);
                        if(targets.size() > maxambig) {
                            hashTable_.clear(i);
                            ++rem;
                            break;
                        }
                    }
                }
            }
        }
        else {
            for(auto i = hashTable_.begin(), e = hashTable_.end(); i != e; ++i) {
                if(!i->empty()) {
                    std::set<const taxon*> taxa;
                    for(auto loc : *i) {
                        taxa.insert(targetLineages[loc.tgt][int(r)]);
                        if(taxa.size() > maxambig) {
                            hashTable_.clear(i);
                            ++rem;
                            break;
                        }
                    }
                }
            }
        }
        return rem;
    }


    //---------------------------------------------------------------
    void wait_until_add_target_complete(gpu_id, const sketcher&) {
        destroy_inserter();
    }

private:
    //---------------------------------------------------------------
    void destroy_inserter() {
        // destroy inserter
        inserter_ = nullptr;
    }

public:
    //---------------------------------------------------------------
    bool valid() {
        return !(inserter_) || ((inserter_) && (inserter_->valid()));
    }


    //---------------------------------------------------------------
    /**
     * @brief adds sketches to database for all windows in sequence
     */
    window_id add_target(gpu_id,
                         const sequence& seq, target_id tgt,
                         const sketcher& targetSketcher)
    {
        if(!inserter_) make_sketch_inserter();

        window_id win = 0;
        targetSketcher.for_each_sketch(seq,
            [&, this] (auto&& sk) {
                if(inserter_->valid()) {
                    //insert sketch into batch
                    auto& sketch = inserter_->next_item();
                    sketch.tgt = tgt;
                    sketch.win = win;
                    sketch.sk = std::move(sk);
                }
                ++win;
            });
        return win;
    }


    //---------------------------------------------------------------
    void add_sketch_batch(const sketch_batch& batch) {
        for(const auto& windowSketch : batch) {
            //insert features from sketch into database
            for(const auto& f : windowSketch.sk) {
                auto it = hashTable_.insert(
                    f, location{windowSketch.win, windowSketch.tgt});
                if(it->size() > maxLocsPerFeature_) {
                    hashTable_.shrink(it, maxLocsPerFeature_);
                }
            }
        }
    }


    //---------------------------------------------------------------
    void make_sketch_inserter() {
        batch_processing_options<window_sketch> execOpt;
        execOpt.batch_size(1000);
        execOpt.queue_size(100);
        execOpt.concurrency(1,1);

        inserter_ = std::make_unique<batch_executor<window_sketch>>( execOpt,
            [&,this](int, const auto& batch) {
                this->add_sketch_batch(batch);
            });
    }


    //---------------------------------------------------------------
    template<class InputIterator>
    void
    accumulate_matches(const sketcher& querySketcher,
                       InputIterator queryBegin, InputIterator queryEnd,
                       matches_sorter& res) const
    {
        querySketcher.for_each_sketch(queryBegin, queryEnd,
            [this, &res] (const auto& sk) {
                for(auto f : sk) {
                    auto locs = hashTable_.find(f);
                    if(locs != hashTable_.end() && locs->size() > 0) {
                        res.locs_.insert(res.locs_.end(), locs->begin(), locs->end());
                        res.offsets_.emplace_back(res.locs_.size());
                    }
                }
            });
    }
    //---------------------------------------------------------------
    void
    accumulate_matches(const sketcher& querySketcher,
                       const sequence& query,
                       matches_sorter& res) const
    {
        using std::begin;
        using std::end;
        accumulate_matches(querySketcher, begin(query), end(query), res);
    }
    //---------------------------------------------------------------
    void
    query_host(const sketcher& querySketcher,
               const sequence& query1, const sequence& query2,
               matches_sorter& res) const
    {
        res.clear();

        accumulate_matches(querySketcher, query1, res);
        accumulate_matches(querySketcher, query2, res);

        res.sort();
    }



    //---------------------------------------------------------------
    friend void read_binary(std::istream& is, host_hashmap& m, gpu_id) {
        read_binary(is, m.hashTable_);
    }

    //---------------------------------------------------------------
    friend void write_binary(std::ostream& os, const host_hashmap& m, gpu_id) {
        write_binary(os, m.hashTable_);
    }


private:
    hash_table hashTable_;
    std::uint64_t maxLocsPerFeature_;

    std::unique_ptr<batch_executor<window_sketch>> inserter_;
};


} // namespace mc


#endif
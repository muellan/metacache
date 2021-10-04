/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2021 Robin Kobus  (kobus@uni-mainz.de)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *****************************************************************************/

#ifndef MC_HOST_HASH_MAP_H_
#define MC_HOST_HASH_MAP_H_

#include "batch_processing.h"
#include "config.h"
#include "hash_multimap.h"
#include "query_handler.h"
#include "stat_combined.h"
#include "taxonomy.h"

#include <set>
#include <vector>


namespace mc {


template<
    class ValueT
>
class host_hashmap
{
public:
    //---------------------------------------------------------------
    using location         = ValueT;
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


public:
    //---------------------------------------------------------------
    explicit
    host_hashmap() :
        maxLoadFactor_(default_max_load_factor()),
        maxLocationsPerFeature_{max_supported_locations_per_feature()},
        hashTables_{},
        sketchers_{},
        inserters_{}
    {}

    host_hashmap(const host_hashmap&) = delete;
    host_hashmap(host_hashmap&& other) :
        maxLoadFactor_{other.maxLoadFactor_},
        maxLocationsPerFeature_{other.maxLocationsPerFeature_},
        hashTables_{std::move(other.hashTables_)},
        sketchers_{std::move(other.sketchers_)},
        inserters_{std::move(other.inserters_)}
    {}

    host_hashmap& operator = (const host_hashmap&) = delete;
    host_hashmap& operator = (host_hashmap&&)      = default;

    ~host_hashmap() {
        for(size_t i = 0; i < inserters_.size(); ++i)
            destroy_inserter(i);
    }


    //---------------------------------------------------------------
    unsigned num_parts() const noexcept { return hashTables_.size(); }


    //---------------------------------------------------------------
    void initialize_build_hash_tables(part_id numParts, const sketcher& targetSketcher) {
        hashTables_.resize(numParts);
        inserters_.resize(numParts);
        sketchers_.resize(numParts, targetSketcher);

        for(auto& hashTable : hashTables_)
            hashTable.max_load_factor(maxLoadFactor_);
    }


    //---------------------------------------------------------------
    static constexpr std::size_t
    max_bucket_size() noexcept {
        return hash_table::max_bucket_size();
    }


    //---------------------------------------------------------------
    void max_load_factor(float lf) {
        if(lf > 1.0f) lf = 1.0f;
        if(lf < 0.1f) lf = 0.1f;

        using std::abs;
        if(abs(maxLoadFactor_ - lf) > 0.00001f) {
            maxLoadFactor_ = lf;
        }

        for(auto& hashTable : hashTables_)
            hashTable.max_load_factor(maxLoadFactor_);
    }
    //-----------------------------------------------------
    float max_load_factor() const noexcept {
        return maxLoadFactor_;
    }
    //-----------------------------------------------------
    static constexpr float default_max_load_factor() noexcept {
        return 0.8f;
    }


    //---------------------------------------------------------------
    bool empty() const noexcept {
        return key_count() < 1;
    }
    //---------------------------------------------------------------
    std::uint64_t key_count() const noexcept {
        std::uint64_t count = 0;
        for(const auto& hashTable : hashTables_)
            count += hashTable.key_count();
        return count;
    }
    //---------------------------------------------------------------
    std::uint64_t value_count() const noexcept {
        std::uint64_t count = 0;
        for(const auto& hashTable : hashTables_)
            count += hashTable.value_count();
        return count;
    }
    //---------------------------------------------------------------
    std::uint64_t bucket_count() const noexcept {
        std::uint64_t count = 0;
        for(const auto& hashTable : hashTables_)
            count += hashTable.bucket_count();
        return count;
    }
    //---------------------------------------------------------------
    std::uint64_t non_empty_bucket_count() const noexcept {
        std::uint64_t count = 0;
        for(const auto& hashTable : hashTables_)
            count += hashTable.non_empty_bucket_count();
        return count;
    }
    //---------------------------------------------------------------
    std::uint64_t dead_feature_count() const noexcept {
        return key_count() - non_empty_bucket_count();
    }

    //---------------------------------------------------------------
    void clear() {
        for(auto& hashTable : hashTables_)
            hashTable.clear();
    }
    //---------------------------------------------------------------
    /**
     * @brief very dangerous! clears hash table without memory deallocation
     */
    void clear_without_deallocation() {
        for(auto& hashTable : hashTables_)
            hashTable.clear_without_deallocation();
    }


    //---------------------------------------------------------------
    statistics_accumulator
    location_list_size_statistics() const {
        auto totalAccumulator = statistics_accumulator{};

        for(part_id part = 0; part < num_parts(); ++part) {
            auto accumulator = statistics_accumulator{};

            for(const auto& bucket : hashTables_[part]) {
                if(!bucket.empty()) {
                    accumulator += bucket.size();
                }
            }

            if(num_parts() > 1) {
                std::cout
                    << "------------------------------------------------\n"
                    << "database part " << part << ":\n"
                    << "buckets              " << hashTables_[part].bucket_count() << '\n'
                    << "bucket size          " << "max: " << accumulator.max()
                                               << " mean: " << accumulator.mean()
                                               << " +/- " << accumulator.stddev()
                                               << " <> " << accumulator.skewness() << '\n'
                    << "features             " << std::uint64_t(accumulator.size()) << '\n'
                    << "dead features        " << hashTables_[part].key_count() -
                                                  hashTables_[part].non_empty_bucket_count() << '\n'
                    << "locations            " << std::uint64_t(accumulator.sum()) << '\n';
                    // << "load                 " << hashTables_[part].load_factor() << '\n';
            }

            totalAccumulator.merge(accumulator);
        }

        return totalAccumulator;
    }


    //---------------------------------------------------------------
    void print_feature_map(std::ostream& os) const {
        for(part_id part = 0; part < num_parts(); ++part) {
            if(num_parts() > 1)
                os << "database part " << part << ":\n";

            for(const auto& bucket : hashTables_[part]) {
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
    }


    //---------------------------------------------------------------
    void print_feature_counts(std::ostream& os) const {
        for(part_id part = 0; part < num_parts(); ++part) {
            if(num_parts() > 1)
                os << "database part " << part << ":\n";

            for(const auto& bucket : hashTables_[part]) {
                if(!bucket.empty()) {
                    os << std::int_least64_t(bucket.key()) << " -> "
                    << std::int_least64_t(bucket.size()) << '\n';
                }
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
        else if(n < maxLocationsPerFeature_) {
            for(auto& hashTable : hashTables_)
                hashTable.shrink_all(n);
        }
        maxLocationsPerFeature_ = n;
    }
    //-----------------------------------------------------
    bucket_size_type
    max_locations_per_feature() const noexcept {
        return maxLocationsPerFeature_;
    }


    /**************************************************************************
     * @details note that features are not really removed, because the hashmap
     *          does not support erasing keys; instead all values belonging to
     *          the key are cleared and the key is kept without values
     */
    feature_count_type
    remove_features_with_more_locations_than(bucket_size_type n)
    {
        feature_count_type rem = 0;

        for(auto& hashTable : hashTables_) {
            for(auto i = hashTable.begin(), e = hashTable.end(); i != e; ++i) {
                if(i->size() > n) {
                    hashTable.clear(i);
                    ++rem;
                }
            }
        }

        return rem;
    }


    //-------------------------------------------------------------------------
    feature_count_type
    remove_ambiguous_features(taxon_rank r, bucket_size_type maxambig,
                              const taxonomy_cache& taxonomy)
    {
        feature_count_type rem = 0;

        if(maxambig == 0) maxambig = 1;

        for(auto& hashTable : hashTables_) {
            if(r == taxon_rank::Sequence) {
                for(auto i = hashTable.begin(), e = hashTable.end(); i != e; ++i) {
                    if(!i->empty()) {
                        std::set<target_id> targets;
                        for(auto loc : *i) {
                            targets.insert(loc.tgt);
                            if(targets.size() > maxambig) {
                                hashTable.clear(i);
                                ++rem;
                                break;
                            }
                        }
                    }
                }
            }
            else {
                for(auto i = hashTable.begin(), e = hashTable.end(); i != e; ++i) {
                    if(!i->empty()) {
                        std::unordered_set<const taxon*> taxa;
                        for(auto loc : *i) {
                            taxa.insert(taxonomy.cached_ancestor(loc.tgt, r));
                            if(taxa.size() > maxambig) {
                                hashTable.clear(i);
                                ++rem;
                                break;
                            }
                        }
                    }
                }
            }
        }

        return rem;
    }


    //---------------------------------------------------------------
    void wait_until_add_target_complete(part_id part, const sketcher&) {
        destroy_inserter(part);
    }

private:
    //---------------------------------------------------------------
    void destroy_inserter(part_id part) {
        inserters_[part] = nullptr;
    }

public:
    //---------------------------------------------------------------
    bool add_target_failed(part_id partId) const noexcept {
        if(inserters_[partId])
            return !(inserters_[partId]->valid());
        else
            return false;
    }

    //---------------------------------------------------------------
    bool check_load_factor(part_id) const {
        // guaranteed by hash table
        return true;
    }

    //---------------------------------------------------------------
    /**
     * @brief adds sketches to database for all windows in sequence
     */
    window_id add_target(part_id part,
                         const sequence& seq, target_id tgt,
                         const sketcher&)
    {
        if(!inserters_[part]) make_sketch_inserter(part);

        window_id win = 0;
        sketchers_[part].for_each_sketch(seq,
            [&, this] (const auto& sk) {
                if(inserters_[part]->valid()) {
                    //insert sketch into batch
                    auto& sketch = inserters_[part]->next_item();
                    sketch.tgt = tgt;
                    sketch.win = win;
                    sketch.sk = sk;
                }
                ++win;
            });

        return win;
    }

private:
    //---------------------------------------------------------------
    void add_sketch_batch(part_id part, const sketch_batch& batch) {
        for(const auto& windowSketch : batch) {
            //insert features from sketch into database
            for(const auto& f : windowSketch.sk) {
                auto it = hashTables_[part].insert(
                    f, location{windowSketch.win, windowSketch.tgt});
                if(it->size() > maxLocationsPerFeature_) {
                    hashTables_[part].shrink(it, maxLocationsPerFeature_);
                }
            }
        }
    }


    //---------------------------------------------------------------
    void make_sketch_inserter(part_id part) {
        batch_processing_options<window_sketch> execOpt;
        execOpt.batch_size(1000);
        execOpt.queue_size(100);
        execOpt.concurrency(1,1);

        inserters_[part] = std::make_unique<batch_executor<window_sketch>>( execOpt,
            [&,this,part](int, const auto& batch) {
                this->add_sketch_batch(part, batch);
                return true;
            });
    }


    //---------------------------------------------------------------
    /**
     * @brief accumulate matches from first db part,
     *        keep sketches in case of multiple db parts
     */
    sketch
    accumulate_matches(const sequence& query1, const sequence& query2,
                       query_handler<location>& queryHandler,
                       bool keepSketches) const
    {
        using std::begin;
        using std::end;

        sketch allWindowSketch;

        const auto& sketcher = queryHandler.querySketcher;
        auto& sorter = queryHandler.matchesSorter;

        sketcher.for_each_sketch(begin(query1), end(query1),
            [&, this] (const auto& sk) {
                for(auto f : sk) {
                    auto locs = hashTables_[0].find(f);
                    if(locs != hashTables_[0].end() && locs->size() > 0) {
                        sorter.append(locs->begin(), locs->end());
                    }
                }

                if(keepSketches)
                    allWindowSketch.insert(allWindowSketch.end(), sk.begin(), sk.end());
            });

        sketcher.for_each_sketch(begin(query2), end(query2),
            [&, this] (const auto& sk) {
                for(auto f : sk) {
                    auto locs = hashTables_[0].find(f);
                    if(locs != hashTables_[0].end() && locs->size() > 0) {
                        sorter.append(locs->begin(), locs->end());
                    }
                }

                if(keepSketches)
                    allWindowSketch.insert(allWindowSketch.end(), sk.begin(), sk.end());
            });

        return allWindowSketch;
    }


    //---------------------------------------------------------------
    /**
     * @brief accumulate matches from other db parts
     */
    void
    accumulate_matches(part_id part,
                       query_handler<location>& queryHandler,
                       const sketch& allWindowSketch) const
    {
        auto& sorter = queryHandler.matchesSorter;

        for(auto f : allWindowSketch) {
            auto locs = hashTables_[part].find(f);
            if(locs != hashTables_[part].end() && locs->size() > 0) {
                sorter.append(locs->begin(), locs->end());
            }
        }
    }

public:
    //---------------------------------------------------------------
    classification_candidates
    query_host_hashmap(const sequence& query1, const sequence& query2,
                       query_handler<location>& queryHandler,
                       const taxonomy_cache& taxonomy,
                       const candidate_generation_rules& rules) const
    {
        auto& sorter = queryHandler.matchesSorter;

        sorter.clear();

        // accumulate and sort matches from first db part
        sorter.next();

        sketch allWindowSketch = accumulate_matches(
            query1, query2, queryHandler, num_parts() > 1);

        sorter.sort();

        // accumulate and sort matches from other db parts
        for(part_id part = 1; part < num_parts(); ++part) {
            sorter.next();

            accumulate_matches(part, queryHandler, allWindowSketch);

            sorter.sort();
        }

        return classification_candidates{taxonomy, queryHandler.matchesSorter.locations(), rules};
    }

    //---------------------------------------------------------------
    void prepare_for_query_hash_tables(part_id numParts, unsigned) {
        hashTables_.resize(numParts);

        for(auto& hashTable : hashTables_)
            hashTable.max_load_factor(maxLoadFactor_);
    }

    //---------------------------------------------------------------
    friend void read_binary(std::istream& is, host_hashmap& m, part_id part,
                            concurrent_progress& readingProgress)
    {
        read_binary(is, m.hashTables_[part], readingProgress);
    }

    //---------------------------------------------------------------
    friend void write_binary(std::ostream& os, const host_hashmap& m, part_id part) {
        write_binary(os, m.hashTables_[part]);
    }


private:
    float maxLoadFactor_;
    std::uint64_t maxLocationsPerFeature_;

    std::vector<hash_table> hashTables_;

    std::vector<sketcher> sketchers_;
    std::vector<std::unique_ptr<batch_executor<window_sketch>>> inserters_;
};


} // namespace mc


#endif

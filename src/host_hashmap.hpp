/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2022 Robin Kobus  (github.com/funatiq)
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

#ifndef MC_HOST_HASH_MAP_HPP_
#define MC_HOST_HASH_MAP_HPP_

#include "batch_processing.hpp"
#include "hash_multimap.hpp"
#include "query_handler.hpp"
#include "stat_combined.hpp"
#include "taxonomy.hpp"
#include "config.hpp"
#include "cmdline_utility.hpp"

#include <set>
#include <vector>


namespace mc {


template <class ValueT>
class host_hashmap
{
public:
    //---------------------------------------------------------------
    using location         = ValueT;
    using bucket_size_type = mc::loclist_size_t;

    using sketch  = typename sketcher::sketch_type;  // range of features
    using feature = typename sketcher::feature_type;

    using table_index = part_id;


private:
    //-----------------------------------------------------
    // / @brief "heart of the database": maps features to target locations
    using hash_table = hash_multimap<feature,location>; // key, value

    //-----------------------------------------------------
    // / @brief needed for batched, asynchonous insertion into feature_store
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
    host_hashmap () :
        maxLoadFactor_(default_max_load_factor()),
        maxLocationsPerFeature_{max_supported_locations_per_feature()},
        tables_{},
        sketchers_{},
        inserters_{}
    {}

    host_hashmap (const host_hashmap&) = delete;

    host_hashmap (host_hashmap&& other) :
        maxLoadFactor_{other.maxLoadFactor_},
        maxLocationsPerFeature_{other.maxLocationsPerFeature_},
        tables_{std::move(other.tables_)},
        sketchers_{std::move(other.sketchers_)},
        inserters_{std::move(other.inserters_)}
    {}

    host_hashmap& operator = (const host_hashmap&) = delete;
    host_hashmap& operator = (host_hashmap&&)      = default;


    //---------------------------------------------------------------
    unsigned table_count () const noexcept { return tables_.size(); }


    //---------------------------------------------------------------
    void initialize_tables (table_index tableCount)
    {
        tables_.resize(tableCount);
        inserters_.resize(tableCount);
        sketchers_.resize(tableCount);

        for (auto& table : tables_) {
            table.max_load_factor(maxLoadFactor_);
        }
    }

    //---------------------------------------------------------------
    void add_empty_table ()
    {
        tables_.emplace_back();
        tables_.back().max_load_factor(maxLoadFactor_);

        inserters_.emplace_back();
        sketchers_.emplace_back();
    }


    //---------------------------------------------------------------
    static constexpr std::size_t
    max_bucket_size () noexcept {
        return hash_table::max_bucket_size();
    }


    //---------------------------------------------------------------
    void max_load_factor (float lf) {
        if (lf > 1.0f) lf = 1.0f;
        if (lf < 0.1f) lf = 0.1f;

        using std::abs;
        if (abs(maxLoadFactor_ - lf) > 0.00001f) {
            maxLoadFactor_ = lf;
        }

        for (auto& table : tables_) {
            table.max_load_factor(maxLoadFactor_);
        }
    }
    //-----------------------------------------------------
    float max_load_factor () const noexcept {
        return maxLoadFactor_;
    }
    //-----------------------------------------------------
    static constexpr float default_max_load_factor () noexcept {
        return 0.8f;
    }


    //---------------------------------------------------------------
    bool empty () const noexcept {
        return key_count() < 1;
    }
    //---------------------------------------------------------------
    std::uint64_t key_count () const noexcept {
        std::uint64_t count = 0;
        for (const auto& table : tables_) {
            count += table.key_count();
        }
        return count;
    }
    //---------------------------------------------------------------
    std::uint64_t value_count () const noexcept {
        std::uint64_t count = 0;
        for (const auto& table : tables_) {
            count += table.value_count();
        }
        return count;
    }
    //---------------------------------------------------------------
    std::uint64_t bucket_count () const noexcept {
        std::uint64_t count = 0;
        for (const auto& table : tables_) {
            count += table.bucket_count();
        }
        return count;
    }
    //---------------------------------------------------------------
    std::uint64_t non_empty_bucket_count () const noexcept {
        std::uint64_t count = 0;
        for (const auto& table : tables_) {
            count += table.non_empty_bucket_count();
        }
        return count;
    }
    //---------------------------------------------------------------
    std::uint64_t dead_feature_count () const noexcept {
        const auto nebs = non_empty_bucket_count();
        const auto keys = key_count();
        return keys > nebs ? keys - nebs : 0;
    }

    //---------------------------------------------------------------
    std::size_t occupied_bytes () const noexcept
    {
        std::size_t bytes = 0;
        for (const auto& table : tables_) { bytes += table.occupied_bytes(); }
        return bytes;
    }
    //---------------------------------------------------------------
    std::size_t occupied_bytes (unsigned tableIdx) const noexcept
    {
        return (tableIdx < tables_.size()) ? tables_[tableIdx].occupied_bytes() : 0;
    }


    //---------------------------------------------------------------
    void clear () {
        for (auto& table : tables_) {
            table.clear();
        }
        inserters_.clear();
        sketchers_.clear();
    }

    //---------------------------------------------------------------
    /**
     * @brief very dangerous! clears hash table without memory deallocation
     */
    void clear_without_deallocation () {
        for (auto& table : tables_) {
            table.clear_without_deallocation();
        }
    }


    //---------------------------------------------------------------
    unsigned remove_empty_tables ()
    {
        unsigned removed = 0;

        auto i = tables_.begin();
        while (i < tables_.end()) {
            if (i->value_count() < 1 || i->key_count() < 1) {
                i = tables_.erase(i);
            } else {
                ++i;
            }
        } 

        return removed;
    }


    //---------------------------------------------------------------
    bool merge_reduce_max_tables_max_bytes (
        unsigned desiredCount, std::size_t maxBytes,
        int maxConcurrency = 1, info_level info = info_level::moderate)
    {
        if (desiredCount < 1) return false;
        if (maxBytes < 1) return false;

        if (desiredCount >  tables_.size()) return false;
        if (desiredCount == tables_.size()) return true;

        if (remove_empty_tables() && desiredCount <= tables_.size()) return true;

        const std::size_t expectedMerges = tables_.size() - desiredCount;
        const bool showInfo = info != info_level::silent;

        if (showInfo) {
            if (info == info_level::verbose) {   
                std::cout << "\n================================================\n"
                          << "merging " << tables_.size() << " tables into min. " << desiredCount 
                          << " parts with max. size of " << (maxBytes/1048576) << " MB):\n";
                std::cout << "------------------------------------------------\n";
                // std::cout << "bucket size: " << sizeof(typename hash_table::bucket_type)
                //           << " value size:  " << sizeof(typename hash_table::value_type) << "\n";
                int pid = 0;
                for (auto const& table : tables_) {
                    std::cout << (pid < 10 ? "[ " : "[") << pid << "]" 
                              << "  keys: " << table.key_count()
                              << "  vals: " << table.value_count()
                              << "  size: " << (table.occupied_bytes()/1048576) << " MB\n";
                    ++pid;
                }
            } else {
                show_progress_indicator(std::cout, 0.0f);
            }
        }

        std::atomic<int> merges {0};

        // merge tables pairwise 
        while (tables_.size() > desiredCount) {
            // sort tables by size: largest -> smallest
            std::sort(tables_.begin(), tables_.end(),
                [](const auto& a, const auto& b){
                    return a.occupied_bytes() > b.occupied_bytes();
                });

            const int numTbl = tables_.size();
            const int numRes = desiredCount;

            std::vector<std::thread> workers;
            const int maxMergeThreads = std::min(maxConcurrency,std::max(numRes,numTbl/2));
            // const int insertThreads = 1 + ((maxConcurrency-maxMergeThreads)/maxMergeThreads);
            workers.reserve(maxMergeThreads);
            
            for (int si = numTbl-1, ti = 0, mergesLeft = numTbl - numRes; 
                si > ti && mergesLeft > 0; --si, ++ti)
            {
                // respect size limit with 10% safety margin
                const auto newSize = static_cast<std::size_t>(
                                        1.1 * (tables_[si].occupied_bytes() +
                                               tables_[ti].occupied_bytes()));
                if (info == info_level::verbose) {   
                    std::cout << "  merge "
                              << "(" << ti << ") " << (tables_[ti].occupied_bytes()/1048576) << " MB + "
                              << "(" << si << ") " << (tables_[si].occupied_bytes()/1048576) << " MB -> "
                              << (newSize/1048576) << " MB\n";
                }
                if (newSize <= maxBytes) {
                    workers.emplace_back(std::thread{[&](int s, int t){
                        tables_[t].insert(std::move(tables_[s]));
                        tables_[s].clear();
                        merges++;
                    }, si, ti});
                    --mergesLeft;
                }
            }
            if (workers.empty()) break;
            for (auto& w : workers) { w.join(); }

            float progress = float(merges) / float(expectedMerges);
            if (showInfo) show_progress_indicator(std::cout, progress);

            remove_empty_tables();
        }
         
        if (showInfo) show_progress_indicator(std::cout, 1.0f);

        // sort all values for each key (important invariant for querying!)
        for (auto& table : tables_) {
            for (auto& b : table) {
                std::sort(b.begin(), b.end());
            }
        }

        if (showInfo) {
            if (info == info_level::verbose) {   
                std::cout << "------------------------------------------------\n";
                int pid = 0;
                for (auto const& table : tables_) {
                    std::cout << (pid < 10 ? "[ " : "[") << pid << "]" 
                              << "  keys: " << table.key_count()
                              << "  vals: " << table.value_count()
                              << "  size: " << (table.occupied_bytes()/1048576) << " MB\n";
                    ++pid;
                }
                std::cout << "================================================\n";
            } else {
                clear_current_line (std::cout);
            }
        }

        return merges > 0;
    }


    //---------------------------------------------------------------
    statistics_accumulator
    location_list_size_statistics () const {
        auto totalAccumulator = statistics_accumulator{};

        for (table_index t = 0; t < table_count(); ++t) {
            auto accumulator = statistics_accumulator{};

            for (const auto& bucket : tables_[t]) {
                if (not bucket.empty()) {
                    accumulator += bucket.size();
                }
            }

            if (table_count() > 1) {
                std::cout
                    << "------------------------------------------------\n"
                    << "database part " << (t+1) << " / " << (table_count()) << ":\n"
                    << "buckets            " << tables_[t].bucket_count() << '\n'
                    << "bucket size        " << "max: " << accumulator.max()
                                             << " mean: " << accumulator.mean()
                                             << " +/- " << accumulator.stddev()
                                             << " <> " << accumulator.skewness() << '\n'
                    << "features           " << std::uint64_t(accumulator.size()) << '\n'
                    << "dead features      " << (tables_[t].key_count() -
                                                 tables_[t].non_empty_bucket_count()) << '\n'
                    << "locations          " << std::uint64_t(accumulator.sum()) << '\n';
                    // << "load               " << tables_[t].load_factor() << '\n';
            }

            totalAccumulator.merge(accumulator);
        }

        return totalAccumulator;
    }


    //---------------------------------------------------------------
    void print_feature_map (std::ostream& os) const {
        for (table_index t = 0; t < table_count(); ++t) {
            if (table_count() > 1)
                os << "database part " << (t+1) << ":\n";

            for (const auto& bucket : tables_[t]) {
                if (not bucket.empty()) {
                    os << std::int_least64_t(bucket.key()) << " -> ";
                    for (location p : bucket) {
                        os << '(' << std::int_least64_t(p.tgt)
                        << ',' << std::int_least64_t(p.win) << ')';
                    }
                    os << '\n';
                }
            }
        }
    }


    //---------------------------------------------------------------
    void print_feature_counts (std::ostream& os) const {
        for (table_index t = 0; t < table_count(); ++t) {
            if (table_count() > 1)
                os << "database part " << (t+1) << ":\n";

            for (const auto& bucket : tables_[t]) {
                if (not bucket.empty()) {
                    os << std::int_least64_t(bucket.key()) << " -> "
                    << std::int_least64_t(bucket.size()) << '\n';
                }
            }
        }
    }


    //---------------------------------------------------------------
    static bucket_size_type
    max_supported_locations_per_feature () noexcept {
        return (hash_table::max_bucket_size() - 1);
    }
    //-----------------------------------------------------
    void max_locations_per_feature (bucket_size_type n)
    {
        if (n < 1) n = 1;
        if (n >= max_supported_locations_per_feature()) {
            n = max_supported_locations_per_feature();
        }
        else if (n < maxLocationsPerFeature_) {
            for (auto& table : tables_) {
                table.shrink_all(n);
            }
        }
        maxLocationsPerFeature_ = n;
    }
    //-----------------------------------------------------
    bucket_size_type
    max_locations_per_feature () const noexcept {
        return maxLocationsPerFeature_;
    }


    //-------------------------------------------------------------------------
    /**
     * @details note that features are not really removed, because the hashmap
     *          does not support erasing keys; instead all values belonging to
     *          the key are cleared and the key is kept without values
     */
    feature_count_type
    remove_features_with_more_locations_than (bucket_size_type n)
    {
        feature_count_type rem = 0;

        for (auto& table : tables_) {
            for (auto i = table.begin(), e = table.end(); i != e; ++i) {
                if (i->size() > n) {
                    table.clear(i);
                    ++rem;
                }
            }
        }

        return rem;
    }


    //-------------------------------------------------------------------------
    feature_count_type
    remove_ambiguous_features (
        taxon_rank r, bucket_size_type maxambig, const taxonomy_cache& taxonomy)
    {
        feature_count_type rem = 0;

        if (maxambig == 0) maxambig = 1;

        for (auto& table : tables_) {
            if (r == taxon_rank::Sequence) {
                for (auto i = table.begin(), e = table.end(); i != e; ++i) {
                    if (not i->empty()) {
                        std::set<target_id> targets;
                        for (auto loc : *i) {
                            targets.insert(loc.tgt);
                            if (targets.size() > maxambig) {
                                table.clear(i);
                                ++rem;
                                break;
                            }
                        }
                    }
                }
            }
            else {
                for (auto i = table.begin(), e = table.end(); i != e; ++i) {
                    if (not i->empty()) {
                        std::unordered_set<const taxon*> taxa;
                        for (auto loc : *i) {
                            taxa.insert(taxonomy.cached_ancestor(loc.tgt, r));
                            if (taxa.size() > maxambig) {
                                table.clear(i);
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
    void wait_until_add_target_complete (table_index idx, const sketching_opt&) {
        inserters_[idx] = nullptr;
    }


public:
    //---------------------------------------------------------------
    bool add_target_failed (table_index idx) const noexcept
    {
        if (inserters_[idx]) {
            return not (inserters_[idx]->valid());
        }
        return false;
    }

    //---------------------------------------------------------------
    bool check_load_factor (table_index) const {
        // guaranteed by hash table
        return true;
    }

    //---------------------------------------------------------------
    /**
     * @brief adds sketches to database for all windows in sequence
     */
    window_id add_target (
        table_index idx, const sequence& seq, target_id tgt, const sketching_opt& opt)
    {
        if (not inserters_[idx]) make_sketch_inserter(idx);

        window_id win = 0;
        sketchers_[idx].for_each_sketch(seq, opt,
            [&, this] (const auto& sk) {
                if (inserters_[idx]->valid()) {
                    // insert sketch into batch
                    auto& sketch = inserters_[idx]->next_item();
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
    void add_sketch_batch (table_index idx, const sketch_batch& batch)
    {
        for (const auto& windowSketch : batch) {
            // insert features from sketch into database
            for (const auto& f : windowSketch.sk) {
                auto it = tables_[idx].insert(
                    f, location{windowSketch.win, windowSketch.tgt});
                if (it->size() > maxLocationsPerFeature_) {
                    tables_[idx].shrink(it, maxLocationsPerFeature_);
                }
            }
        }
    }


    //---------------------------------------------------------------
    void make_sketch_inserter (table_index idx)
    {
        batch_processing_options<window_sketch> execOpt;
        execOpt.batch_size(1000);
        execOpt.queue_size(100);
        execOpt.concurrency(1,1);

        inserters_[idx] = std::make_unique<batch_executor<window_sketch>>( execOpt,
            [&,this,idx](int, const auto& batch) {
                this->add_sketch_batch(idx, batch);
                return true;
            });
    }


    //---------------------------------------------------------------
    /**
     * @brief accumulate matches from first db part,
     *        keep sketches in case of multiple db parts
     */
    sketch
    accumulate_matches (
        const sequence& query1, const sequence& query2,
        query_handler<location>& queryHandler,
        const sketching_opt& opt,
        bool keepSketches) const
    {
        using std::begin;
        using std::end;

        sketch allWindowSketch;

        auto& sketcher = queryHandler.querySketcher;
        auto& sorter = queryHandler.matchesSorter;

        sketcher.for_each_sketch(begin(query1), end(query1), opt,
            [&, this] (const auto& sk) {
                for (auto f : sk) {
                    auto locs = tables_[0].find(f);
                    if (locs != tables_[0].end() && locs->size() > 0) {
                        sorter.append(locs->begin(), locs->end());
                    }
                }

                if (keepSketches)
                    allWindowSketch.insert(allWindowSketch.end(), sk.begin(), sk.end());
            });

        sketcher.for_each_sketch(begin(query2), end(query2), opt,
            [&, this] (const auto& sk) {
                for (auto f : sk) {
                    auto locs = tables_[0].find(f);
                    if (locs != tables_[0].end() && locs->size() > 0) {
                        sorter.append(locs->begin(), locs->end());
                    }
                }

                if (keepSketches)
                    allWindowSketch.insert(allWindowSketch.end(), sk.begin(), sk.end());
            });

        return allWindowSketch;
    }


    //---------------------------------------------------------------
    /**
     * @brief accumulate matches from other table
     */
    void accumulate_matches (
        table_index idx,
        query_handler<location>& queryHandler,
        const sketch& allWindowSketch) const
    {
        auto& sorter = queryHandler.matchesSorter;

        for (auto f : allWindowSketch) {
            auto locs = tables_[idx].find(f);
            if (locs != tables_[idx].end() && locs->size() > 0) {
                sorter.append(locs->begin(), locs->end());
            }
        }
    }

public:
    //---------------------------------------------------------------
    void query_host_hashmap (
        const sequence& query1, const sequence& query2,
        query_handler<location>& queryHandler,
        const taxonomy_cache& taxonomy,
        const sketching_opt& opt,
        const candidate_generation_rules& rules) const
    {
        queryHandler.clear();

        auto& sorter = queryHandler.matchesSorter;

        // accumulate and sort matches from first table
        sorter.next();

        sketch allWindowSketch = accumulate_matches(
            query1, query2, queryHandler, opt, table_count() > 1);

        sorter.sort();

        // accumulate and sort matches from other tables
        for (table_index idx = 1; idx < table_count(); ++idx) {
            sorter.next();
            accumulate_matches(idx, queryHandler, allWindowSketch);
            sorter.sort();
        }

        queryHandler.classificationCandidates.insert(
            taxonomy, queryHandler.matchesSorter.locations(), rules);
    }

    //---------------------------------------------------------------
    void prepare_query_tables (table_index tableCount, unsigned) 
    {
        tables_.resize(tableCount);

        for (auto& table : tables_) {
            table.max_load_factor(maxLoadFactor_);
        }
    }

    //---------------------------------------------------------------
    friend void read_binary (
        std::istream& is, host_hashmap& m, table_index idx,
        concurrent_progress& readingProgress)
    {
        read_binary(is, m.tables_[idx], readingProgress);
    }

    //---------------------------------------------------------------
    friend void write_binary (
        std::ostream& os, const host_hashmap& m, table_index idx)
    {
        write_binary(os, m.tables_[idx]);
    }


private:
    float maxLoadFactor_;
    std::uint64_t maxLocationsPerFeature_;

    std::vector<hash_table> tables_;

    std::vector<sketcher> sketchers_;
    std::vector<std::unique_ptr<batch_executor<window_sketch>>> inserters_;
};


} // namespace mc


#endif

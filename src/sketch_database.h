/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2020 André Müller (muellan@uni-mainz.de)
 *                       & Robin Kobus  (kobus@uni-mainz.de)
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

#ifndef MC_SKETCH_DATABASE_H_
#define MC_SKETCH_DATABASE_H_

#include <climits>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <algorithm>
#include <type_traits>
#include <cstdint>
#include <string>
#include <limits>
#include <memory>
#include <future>
#include <chrono>

#include "version.h"
#include "config.h"
#include "io_error.h"
#include "io_options.h"
#include "stat_combined.h"
#include "taxonomy.h"
#include "hash_multimap.h"
#include "dna_encoding.h"
#include "typename.h"

#include "../dep/queue/readerwriterqueue.h"


namespace mc {

/*************************************************************************//**
 *
 * @brief  maps 'features' (e.g. hash values obtained by min-hashing)
 *         to 'locations' = positions in targets/reference sequences
 *
 *         not copyable, but movable
 *
 * @details terminology
 *  target          reference sequence whose sketches are stored in the DB
 *
 *  taxon           one node in a taxonomic hierarchy (sequence metadata is
 *                  also stored in taxa)
 *
 *  taxon_id        numeric taxon identifier
 *  taxon_name      alphanumeric sequence identifier (e.g. an NCBI accession)
 *
 *  query           sequence (usually short reads) that shall be matched against
 *                  the reference targets
 *
 *  window_id       window index (starting with 0) within a target
 *
 *  location        (window_id, target id) = "window within a target sequence"
 *
 *  full_lineage    path from root -> lowest taxon (vector of const taxon*);
 *                  may have arbitrary length
 *
 *  ranked_lineage  main rank taxon ids (array of const taxon*);
 *                  has fixed length; each index always refers to the same rank
 *
 *  sketcher        function object type, that maps reference sequence
 *                      windows (= sequence interval) to
 *                      sketches (= collection of features of the same type)
 *
 *  feature         single part of a sketch
 *
 *  feature_hash    hash function for (feature -> location) map
 *
 *  target id        internal target (reference sequence) identification

 *  window id        target window identification
 *
 *  loclist_size_t   bucket (location list) size tracking type
 *
 *****************************************************************************/
class database
{
public:
    //---------------------------------------------------------------
    // from global config
    using sequence         = mc::sequence;
    using sketcher         = mc::sketcher;
    using feature_hash     = mc::feature_hash;
    using target_id        = mc::target_id;
    using window_id        = mc::window_id;
    using bucket_size_type = mc::loclist_size_t;
    //-----------------------------------------------------
    // from taxonomy
    using taxon_id       = taxonomy::taxon_id;
    using taxon          = taxonomy::taxon;
    using taxon_rank     = taxonomy::rank;
    using taxon_name     = taxonomy::taxon_name;
    using ranked_lineage = taxonomy::ranked_lineage;
    using full_lineage   = taxonomy::full_lineage;
    using file_source    = taxonomy::file_source;
    using taxon_iterator = taxonomy::const_iterator;
    using taxon_range    = taxonomy::const_range;

    //-----------------------------------------------------
    using match_count_type = std::uint16_t;

    //---------------------------------------------------------------
    enum class scope { everything, metadata_only };


    //-----------------------------------------------------
    class target_limit_exceeded_error : public std::runtime_error {
    public:
        target_limit_exceeded_error():
            std::runtime_error{"target count limit exceeded"}
        {}
    };

    const taxon* taxon_of_target(target_id id) const {return targets_[id]; }

    //-----------------------------------------------------
    /** @brief internal location representation = (window index, target index)
     *         these are stored in the in-memory database and on disk
     */
    #pragma pack(push, 1)
    //avoid padding bits
    struct location
    {
        window_id win;
        target_id tgt;

        friend bool
        operator == (const location& a, const location& b) noexcept {
            return (a.tgt == b.tgt) && (a.win == b.win);
        }

        friend bool
        operator < (const location& a, const location& b) noexcept {
            if(a.tgt < b.tgt) return true;
            if(a.tgt > b.tgt) return false;
            return (a.win < b.win);
        }
    };
    //avoid padding bits
    #pragma pack(pop)

    using match_locations = std::vector<location>;


    //-----------------------------------------------------
    using sketch  = typename sketcher::sketch_type;  //range of features
    using feature = typename sketch::value_type;


private:
    static constexpr std::size_t insert_batch_size() noexcept { return 10000; };

    //use negative numbers for sequence level taxon ids
    static constexpr taxon_id
    taxon_id_of_target(target_id id) noexcept { return -taxon_id(id)-1; }


    //-----------------------------------------------------
    /// @brief "heart of the database": maps features to target locations
    using feature_store = hash_multimap<feature,location, //key, value
                              feature_hash,               //key hasher
                              std::equal_to<feature>,     //key comparator
                              chunk_allocator<location>,  //value allocator
                              std::allocator<feature>,    //bucket+key allocator
                              bucket_size_type>;          //location list size


    //-----------------------------------------------------
    /// @brief needed for batched, asynchonous insertion into feature_store
    struct window_sketch
    {
        window_sketch() = default;

        window_sketch(target_id tgt, window_id win, sketch sk) :
            tgt{tgt}, win{win}, sk{std::move(sk)} {};

        target_id tgt;
        window_id win;
        sketch sk;
    };

    using sketch_batch = std::vector<window_sketch>;
    using sketch_queue = moodycamel::BlockingReaderWriterQueue<sketch_batch>;


public:
    //---------------------------------------------------------------
    using feature_count_type = typename feature_store::size_type;


    //---------------------------------------------------------------
    /** @brief used for query result storage/accumulation
     */
    class matches_sorter {
        friend class database;

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


    //---------------------------------------------------------------
    explicit
    database(sketcher targetSketcher = sketcher{}) :
        database{targetSketcher, targetSketcher}
    {}
    //-----------------------------------------------------
    explicit
    database(sketcher targetSketcher, sketcher querySketcher) :
        targetSketcher_{std::move(targetSketcher)},
        querySketcher_{std::move(querySketcher)},
        maxLocsPerFeature_(max_supported_locations_per_feature()),
        features_{},
        targets_{},
        taxa_{},
        ranksCache_{taxa_, taxon_rank::Sequence},
        targetLineages_{taxa_},
        name2tax_{},
        batch_{},
        queue_(10),
        sketchingDone_{true},
        inserterThread_{}
    {
        features_.max_load_factor(default_max_load_factor());
    }

    database(const database&) = delete;
    database(database&& other) :
        targetSketcher_{std::move(other.targetSketcher_)},
        querySketcher_{std::move(other.querySketcher_)},
        maxLocsPerFeature_(other.maxLocsPerFeature_),
        features_{std::move(other.features_)},
        targets_{std::move(other.targets_)},
        taxa_{std::move(other.taxa_)},
        ranksCache_{std::move(other.ranksCache_)},
        targetLineages_{std::move(other.targetLineages_)},
        name2tax_{std::move(other.name2tax_)},
        batch_{std::move(other.batch_)},
        queue_{std::move(other.queue_)},
        sketchingDone_(other.sketchingDone_.load()),
        inserterThread_{std::move(other.inserterThread_)}
    {}

    database& operator = (const database&) = delete;
    database& operator = (database&&)      = default;

    ~database() {
        if(inserterThread_ && inserterThread_->valid()) {
            inserterThread_->get();
        }
    }


    //---------------------------------------------------------------
    /**
     * @return const ref to the object that transforms target sequence
     *         snippets into a collection of features
     */
    const sketcher&
    target_sketcher() const noexcept {
        return targetSketcher_;
    }
    //-----------------------------------------------------
    /**
     * @return const ref to the object that transforms query sequence
     *         snippets into a collection of features
     */
    const sketcher&
    query_sketcher() const noexcept {
        return querySketcher_;
    }
    //-----------------------------------------------------
    /**
     * @brief sets the object that transforms query sequence
     *        snippets into a collection of features
     */
    void query_sketcher(const sketcher& s) {
        querySketcher_ = s;
    }
    void query_sketcher(sketcher&& s) {
        querySketcher_ = std::move(s);
    }


    //---------------------------------------------------------------
    void max_locations_per_feature(bucket_size_type n)
    {
        if(n < 1) n = 1;
        if(n >= max_supported_locations_per_feature()) {
            n = max_supported_locations_per_feature();
        }
        else if(n < maxLocsPerFeature_) {
            for(auto i = features_.begin(), e = features_.end(); i != e; ++i) {
                if(i->size() > n) features_.shrink(i, n);
            }
        }
        maxLocsPerFeature_ = n;
    }
    //-----------------------------------------------------
    bucket_size_type
    max_locations_per_feature() const noexcept {
        return maxLocsPerFeature_;
    }
    //-----------------------------------------------------
    static bucket_size_type
    max_supported_locations_per_feature() noexcept {
        return (feature_store::max_bucket_size() - 1);
    }

    //-----------------------------------------------------
    feature_count_type
    remove_features_with_more_locations_than(bucket_size_type n)
    {
        //note that features are not really removed, because the hashmap
        //does not support erasing keys; instead all values belonging to
        //the key are cleared and the key is kept without values
        feature_count_type rem = 0;
        for(auto i = features_.begin(), e = features_.end(); i != e; ++i) {
            if(i->size() > n) {
                features_.clear(i);
                ++rem;
            }
        }
        return rem;
    }


    //---------------------------------------------------------------
    /**
     * @brief  removes features that have more than 'maxambig' different
     *         taxa on a certain taxonomic rank
     *         e.g. remove features that are present in more than 4 phyla
     *
     * @return number of features (hash table buckets) that were removed
     */
    feature_count_type
    remove_ambiguous_features(taxon_rank r, bucket_size_type maxambig)
    {
        feature_count_type rem = 0;

        if(taxa_.empty()) {
            throw std::runtime_error{"no taxonomy available!"};
        }

        if(maxambig == 0) maxambig = 1;

        if(r == taxon_rank::Sequence) {
            for(auto i = features_.begin(), e = features_.end(); i != e; ++i) {
                if(!i->empty()) {
                    std::set<target_id> targets;
                    for(auto loc : *i) {
                        targets.insert(loc.tgt);
                        if(targets.size() > maxambig) {
                            features_.clear(i);
                            ++rem;
                            break;
                        }
                    }
                }
            }
        }
        else {
            for(auto i = features_.begin(), e = features_.end(); i != e; ++i) {
                if(!i->empty()) {
                    std::set<const taxon*> taxa;
                    for(auto loc : *i) {
                        taxa.insert(targetLineages_[loc.tgt][int(r)]);
                        if(taxa.size() > maxambig) {
                            features_.clear(i);
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
    bool add_target(const sequence& seq, taxon_name sid,
                    taxon_id parentTaxid = 0,
                    file_source source = file_source{})
    {
        using std::begin;
        using std::end;

        //reached hard limit for number of targets
        if(targets_.size() >= max_target_count()) {
            throw target_limit_exceeded_error{};
        }

        if(seq.empty()) return false;

        //don't allow non-unique sequence ids
        if(name2tax_.find(sid) != name2tax_.end()) return false;

        const auto targetCount = target_id(targets_.size());
        const auto taxid = taxon_id_of_target(targetCount);

        //sketch sequence -> insert features
        source.windows = add_all_window_sketches(seq, targetCount);

        //insert sequence metadata as a new taxon
        if(parentTaxid < 1) parentTaxid = 0;
        auto nit = taxa_.emplace(
            taxid, parentTaxid, sid,
            taxon_rank::Sequence, std::move(source));

        //should never happen
        if(nit == taxa_.end()) {
            throw std::runtime_error{"target taxon could not be created"};
        }

        //allows lookup via sequence id (e.g. NCBI accession number)
        const taxon* newtax = &(*nit);
        name2tax_.insert({std::move(sid), newtax});

        //target id -> taxon lookup table
        targets_.push_back(newtax);

        targetLineages_.mark_outdated();

        return true;
    }


    //---------------------------------------------------------------
    std::uint64_t
    target_count() const noexcept {
        return targets_.size();
    }
    static constexpr std::uint64_t
    max_target_count() noexcept {
        return std::numeric_limits<target_id>::max();
    }
    static constexpr std::uint64_t
    max_windows_per_target() noexcept {
        return std::numeric_limits<window_id>::max();
    }

    //-----------------------------------------------------
    bool empty() const noexcept {
        return features_.empty();
    }


    //---------------------------------------------------------------
    void clear() {
        ranksCache_.clear();
        targetLineages_.clear();
        name2tax_.clear();
        features_.clear();
    }


    //-----------------------------------------------------
    /**
     * @brief very dangerous! clears feature map without memory deallocation
     */
    void clear_without_deallocation() {
        ranksCache_.clear();
        targetLineages_.clear();
        name2tax_.clear();
        features_.clear_without_deallocation();
    }


    //---------------------------------------------------------------
    static constexpr taxon_id no_taxon_id() noexcept {
        return taxonomy::none_id();
    }


    //---------------------------------------------------------------
    const taxon*
    taxon_with_id(taxon_id id) const noexcept {
        return taxa_[id];
    }
    //-----------------------------------------------------
    /**
     * @brief will only find sequence-level taxon names == sequence id
     */
    const taxon*
    taxon_with_name(const taxon_name& name) const noexcept {
        if(name.empty()) return nullptr;
        auto i = name2tax_.find(name);
        if(i == name2tax_.end()) return nullptr;
        return i->second;
    }
    //-----------------------------------------------------
    /**
     * @brief will find sequence-level taxon names with different versions
     */
    const taxon*
    taxon_with_similar_name(const taxon_name& name) const noexcept {
        if(name.empty()) return nullptr;
        auto i = name2tax_.upper_bound(name);
        if(i == name2tax_.end()) return nullptr;
        const auto s = name.size();
        if(0 != i->first.compare(0,s,name)) return nullptr;
        return i->second;
    }


    //---------------------------------------------------------------
    void reset_taxa_above_sequence_level(taxonomy&& tax) {
        taxa_.erase_above(taxon_rank::Sequence);
        for(auto& t : tax) {
            taxa_.insert_or_replace(std::move(t));
        }
        //re-initialize ranks cache
        mark_cached_lineages_outdated();
        update_cached_lineages(taxon_rank::Sequence);
    }


    //---------------------------------------------------------------
    std::uint64_t
    non_target_taxon_count() const noexcept {
        return taxa_.size() - targets_.size();
    }
    //-----------------------------------------------------
    taxon_range taxa() const {
        return taxa_.full_range();
    }
    taxon_range non_target_taxa() const {
        return taxa_.subrange_from(1);
    }
    //-----------------------------------------------------
    taxon_range target_taxa() const {
        return taxa_.subrange_until(0);
    }


    //---------------------------------------------------------------
    void reset_parent(const taxon& tax, taxon_id parentId) {
        taxa_.reset_parent(tax.id(), parentId);
        mark_cached_lineages_outdated();
    }
    //-----------------------------------------------------
    void reset_parent(const taxon& tax, const taxon& parent) {
        taxa_.reset_parent(tax.id(), parent.id());
        mark_cached_lineages_outdated();
    }

    //---------------------------------------------------------------
    void mark_cached_lineages_outdated() {
        ranksCache_.mark_outdated();
        targetLineages_.mark_outdated();
    }

    //---------------------------------------------------------------
    void update_cached_lineages(taxon_rank rank) const {
        ranksCache_.update(rank);
        targetLineages_.update(target_count());
    }

    //-----------------------------------------------------
    full_lineage
    lineage(const taxon* tax) const noexcept {
        return tax ? lineage(*tax) : full_lineage();
    }
    full_lineage
    lineage(const taxon& tax) const noexcept {
        return taxa_.lineage(tax);
    }
    //-----------------------------------------------------
    const ranked_lineage&
    ranks(const taxon* tax) const noexcept {
        return ranksCache_[tax];
    }
    const ranked_lineage&
    ranks(const taxon& tax) const noexcept {
        return ranksCache_[tax];
    }
    const ranked_lineage&
    ranks(target_id tgt) const noexcept {
        return targetLineages_[tgt];
    }

    //-----------------------------------------------------
    const taxon*
    parent(const taxon* tax) const noexcept {
        return tax ? parent(*tax) : nullptr;
    }
    const taxon*
    parent(const taxon& tax) const noexcept {
        return taxa_.parent(tax);
    }
    //-----------------------------------------------------
    const taxon*
    ancestor(const taxon* tax, taxon_rank r) const noexcept {
        return tax ? ancestor(*tax,r) : nullptr;
    }
    const taxon*
    ancestor(const taxon& tax, taxon_rank r) const noexcept {
        return ranksCache_[tax][int(r)];
    }
    const taxon*
    ancestor(target_id tgt, taxon_rank r) const noexcept {
        return targetLineages_[tgt][int(r)];
    }

    //-----------------------------------------------------
    const taxon*
    next_ranked_ancestor(const taxon* tax) const noexcept {
        return tax ? next_ranked_ancestor(*tax) : nullptr;
    }
    const taxon*
    next_ranked_ancestor(const taxon& tax) const noexcept {
        if(tax.rank() != taxon_rank::none) return &tax;
        for(const taxon* a : ranksCache_[tax]) {
            if(a) return a;
        }
        return nullptr;
    }
    const taxon*
    lowest_ranked_ancestor(target_id tgt, taxon_rank lowest) const noexcept {
        const auto& lineage = targetLineages_[tgt];
        for(int rank = int(lowest); rank < int(taxon_rank::none); ++rank) {
            if(lineage[rank])
                return lineage[rank];
        }
        return nullptr;
    }

    //---------------------------------------------------------------
    const taxon*
    lca(const taxon* ta, const taxon* tb) const {
        return (ta && tb) ? lca(*ta,*tb) : nullptr;
    }
    const taxon*
    lca(const taxon& ta, const taxon& tb) const {
        return taxa_.lca(ta,tb);
    }

    //---------------------------------------------------------------
    const taxon*
    ranked_lca(const taxon* ta, const taxon* tb) const {
        return (ta && tb) ? ranked_lca(*ta,*tb) : nullptr;
    }
    const taxon*
    ranked_lca(const taxon& ta, const taxon& tb) const {
        return ranksCache_.ranked_lca(ta,tb);
    }
    const taxon*
    ranked_lca(target_id ta, target_id tb, taxon_rank lowest) const {
        return targetLineages_.ranked_lca(ta, tb, lowest);
    }


    //---------------------------------------------------------------
    /**
     * @return number of times a taxon is covered by any target in the DB
     */
    std::uint_least64_t
    coverage_of_taxon(const taxon* covered) const {
        return covered ? coverage_of_taxon(*covered) : 0;
    }
    //-----------------------------------------------------
    std::uint_least64_t
    coverage_of_taxon(const taxon& covered) const {
        auto coverage = std::uint_least64_t(0);

        for(const auto& tax : taxa_) {
            if(tax.rank() == taxon_rank::Sequence) {
                for(const taxon* t : taxa_.lineage(tax)) {
                    if(t == &covered) ++coverage;
                }
            }
        }
        return coverage;
    }

    //-----------------------------------------------------
    /**
     * @return true, if taxon is covered by any target in the DB
     */
    bool covers(const taxon* covered) const {
        return covered ? covers(*covered) : false;
    }
    //-----------------------------------------------------
    bool covers(const taxon& covered) const {
        for(const auto& tax : taxa_) {
            if(tax.rank() == taxon_rank::Sequence) {
                for(const taxon* t : taxa_.lineage(tax)) {
                    if(t == &covered) return true;
                }
            }
        }
        return false;
    }


    //---------------------------------------------------------------
    template<class InputIterator>
    void
    accumulate_matches(InputIterator queryBegin, InputIterator queryEnd,
                       matches_sorter& res) const
    {
        querySketcher_.for_each_sketch(queryBegin, queryEnd,
            [this, &res] (auto sk) {
                 res.offsets_.reserve(res.offsets_.size() + sk.size());

                for(auto f : sk) {
                    auto locs = features_.find(f);
                    if(locs != features_.end() && locs->size() > 0) {
                        res.locs_.insert(res.locs_.end(), locs->begin(), locs->end());
                        res.offsets_.emplace_back(res.locs_.size());
                    }
                }
            });
    }

    //---------------------------------------------------------------
    void
    accumulate_matches(const sequence& query,
                       matches_sorter& res) const
    {
        using std::begin;
        using std::end;
        accumulate_matches(begin(query), end(query), res);
    }


    //---------------------------------------------------------------
    void max_load_factor(float lf) {
        features_.max_load_factor(lf);
    }
    //-----------------------------------------------------
    float max_load_factor() const noexcept {
        return features_.max_load_factor();
    }
    //-----------------------------------------------------
    static constexpr float default_max_load_factor() noexcept {
        return 0.8f;
    }


    //---------------------------------------------------------------
    /**
     * @brief   read database from binary file
     * @details Note that the map is not just de-serialized but
     *          rebuilt by inserting individual keys and values
     *          This should make DB files more robust against changes in the
     *          internal mapping structure.
     */
    void read(const std::string& filename, scope what = scope::everything)
    {
        std::ifstream is{filename, std::ios::in | std::ios::binary};

        if(!is.good()) {
            throw file_access_error{"can't open file " + filename};
        }

        //database version info
        using std::uint64_t;
        using std::uint8_t;
        uint64_t dbVer = 0;
        read_binary(is, dbVer);

        if(uint64_t( MC_DB_VERSION ) != dbVer) {
            throw file_read_error{
                "Database " + filename + " (version " + std::to_string(dbVer) + ")"
                + " is incompatible\nwith this version of MetaCache"
                + " (uses version " + std::to_string(MC_DB_VERSION) + ")" };
        }

        //data type info
        {
            //data type widths
            uint8_t featureSize = 0; read_binary(is, featureSize);
            uint8_t targetSize = 0;  read_binary(is, targetSize);
            uint8_t windowSize = 0;  read_binary(is, windowSize);
            uint8_t bucketSize = 0;  read_binary(is, bucketSize);
            uint8_t taxidSize = 0;   read_binary(is, taxidSize);
            uint8_t numTaxRanks = 0; read_binary(is, numTaxRanks);

            if( (sizeof(feature) != featureSize) ||
                (sizeof(target_id) != targetSize) ||
                (sizeof(bucket_size_type) != bucketSize) ||
                (sizeof(window_id) != windowSize) )
            {
                throw file_read_error{
                    "Database " + filename +
                    " is incompatible with this variant of MetaCache" +
                    " due to different data type sizes"};
            }

            if( (sizeof(taxon_id) != taxidSize) ||
                (taxonomy::num_ranks != numTaxRanks) )
            {
                throw file_read_error{
                    "Database " + filename +
                    " is incompatible with this variant of MetaCache" +
                    " due to different taxonomy data types"};
            }
        }

        clear();

        //sketching parameters
        read_binary(is, targetSketcher_);
        read_binary(is, querySketcher_);

        //target insertion parameters
        read_binary(is, maxLocsPerFeature_);

        //taxon metadata
        read_binary(is, taxa_);

        target_id targetCount = 0;
        read_binary(is, targetCount);
        if(targetCount < 1) return;

        //update target id -> target taxon lookup table
        targets_.reserve(targetCount);
        for(decltype(targetCount) i = 0 ; i < targetCount; ++i) {
            targets_.push_back(taxa_[taxon_id_of_target(i)]);
        }

        //sequence id lookup
        name2tax_.clear();
        for(const auto& t : taxa_) {
            if(t.rank() == taxon_rank::Sequence) {
                name2tax_.insert({t.name(), &t});
            }
        }

        mark_cached_lineages_outdated();
        update_cached_lineages(taxon_rank::Sequence);

        if(what == scope::metadata_only) return;

        //hash table
        read_binary(is, features_);
    }


    //-------------------------------------------------------------------
    /**
     * @brief   write database to binary file
     */
    void write(const std::string& filename) const
    {
        using std::uint64_t;
        using std::uint8_t;

        std::ofstream os{filename, std::ios::out | std::ios::binary};

        if(!os.good()) {
            throw file_access_error{"can't open file " + filename};
        }

        //database version info
        write_binary(os, uint64_t( MC_DB_VERSION ));

        //data type widths
        write_binary(os, uint8_t(sizeof(feature)));
        write_binary(os, uint8_t(sizeof(target_id)));
        write_binary(os, uint8_t(sizeof(window_id)));
        write_binary(os, uint8_t(sizeof(bucket_size_type)));
        write_binary(os, uint8_t(sizeof(taxon_id)));
        write_binary(os, uint8_t(taxonomy::num_ranks));

        //sketching parameters
        write_binary(os, targetSketcher_);
        write_binary(os, querySketcher_);

        //target insertion parameters
        write_binary(os, maxLocsPerFeature_);

        //taxon & target metadata
        write_binary(os, taxa_);
        write_binary(os, target_id(targets_.size()));

        //hash table
        write_binary(os, features_);
    }


    //---------------------------------------------------------------
    std::uint64_t bucket_count() const noexcept {
        return features_.bucket_count();
    }
    //---------------------------------------------------------------
    std::uint64_t feature_count() const noexcept {
        return features_.key_count();
    }
    //---------------------------------------------------------------
    std::uint64_t dead_feature_count() const noexcept {
        return features_.key_count() - features_.non_empty_bucket_count();
    }
    //---------------------------------------------------------------
    std::uint64_t location_count() const noexcept {
        return features_.value_count();
    }


    //---------------------------------------------------------------
    statistics_accumulator
    location_list_size_statistics() const {
        auto priSize = statistics_accumulator{};

        for(const auto& bucket : features_) {
            if(!bucket.empty()) {
                priSize += bucket.size();
            }
        }

        return priSize;
    }


    //---------------------------------------------------------------
    void print_feature_map(std::ostream& os) const {
        for(const auto& bucket : features_) {
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
        for(const auto& bucket : features_) {
            if(!bucket.empty()) {
                os << std::int_least64_t(bucket.key()) << " -> "
                   << std::int_least64_t(bucket.size()) << '\n';
            }
        }
    }


    //---------------------------------------------------------------
    void wait_until_add_target_complete() {
        enqueue_insertion_of_sketch_batch();
        sketchingDone_.store(true);
        if(inserterThread_ && inserterThread_->valid()) {
            inserterThread_->get();
        }
    }


    //---------------------------------------------------------------
    bool add_target_failed() {
        if(!inserterThread_) return false;

        if(inserterThread_->valid()) {
            auto status = inserterThread_->wait_for(std::chrono::minutes::zero());
            if(status == std::future_status::timeout) return false;
        }
        return true;
    }


private:
    //---------------------------------------------------------------
    window_id add_all_window_sketches(const sequence& seq, target_id tgt)
    {
        window_id win = 0;
        targetSketcher_.for_each_sketch(seq,
            [&, this] (auto sk) {
                //insert sketch into batch
                batch_.emplace_back(tgt, win, std::move(sk));

                //enqueue full batch and reset
                if(batch_.size() >= insert_batch_size()) {
                    this->enqueue_insertion_of_sketch_batch();
                    batch_.clear();
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
                auto it = features_.insert(
                    f, location{windowSketch.win, windowSketch.tgt});
                if(it->size() > maxLocsPerFeature_) {
                    features_.shrink(it, maxLocsPerFeature_);
                }
            }
        }
    }


    //---------------------------------------------------------------
    void spawn_sketch_inserter() {
        if(inserterThread_) return;

        inserterThread_ = std::make_unique<std::future<void>>(
            std::async(std::launch::async, [&]() {
                sketch_batch batch;

                while(!sketchingDone_.load() || queue_.peek()) {
                    if(queue_.wait_dequeue_timed(batch, std::chrono::milliseconds(10))) {
                        add_sketch_batch(batch);
                    }
                }

            })
        );
    }


    //---------------------------------------------------------------
    void enqueue_insertion_of_sketch_batch() {
        if(batch_.empty()) return;

        if(sketchingDone_.load()) {
            sketchingDone_.store(false);
            spawn_sketch_inserter();
        }

        //allow queue to grow
        queue_.enqueue(batch_);
    }


    /*************************************************************************//**
    *
    * @brief concurrency-safe ranked lineage cache
    *
    *****************************************************************************/
    class ranked_lineages_of_targets
    {
    public:
        //---------------------------------------------------------------
        using ranked_lineage = taxonomy::ranked_lineage;
        using taxon_rank     = taxonomy::rank;

    public:
        //---------------------------------------------------------------
        explicit
        ranked_lineages_of_targets(const taxonomy& taxa)
        :
            taxa_(taxa),
            lins_{},
            outdated_(true)
        {}

        //---------------------------------------------------------------
        ranked_lineages_of_targets(const ranked_lineages_of_targets&) = delete;

        ranked_lineages_of_targets(ranked_lineages_of_targets&& src):
            taxa_(src.taxa_),
            lins_{std::move(src.lins_)},
            outdated_(src.outdated_)
        {}

        ranked_lineages_of_targets& operator = (const ranked_lineages_of_targets&) = delete;
        ranked_lineages_of_targets& operator = (ranked_lineages_of_targets&&) = delete;


        //---------------------------------------------------------------
        void mark_outdated() {
            outdated_ = true;
        }

        //---------------------------------------------------------------
        void update(target_id numTargets = 0) {
            if(!outdated_) return;

            if(numTargets == 0) numTargets = lins_.size();
            lins_.clear();
            for(size_t tgt = 0; tgt < numTargets; ++tgt) {
                lins_.emplace_back(taxa_.ranks(taxon_id_of_target(tgt)));
            }
            outdated_ = false;
        }

        //-----------------------------------------------------
        void clear() {
            lins_.clear();
            outdated_ = true;
        }

        //---------------------------------------------------------------
        const ranked_lineage&
        operator [] (target_id tgt) const {
            assert(outdated_ == false);
            return lins_[tgt];
        }

        //---------------------------------------------------------------
        const taxon*
        ranked_lca(target_id a, target_id b, taxon_rank lowest) const {
            assert(outdated_ == false);
            return taxa_.ranked_lca(lins_[a], lins_[b], lowest);
        }

    private:
        //---------------------------------------------------------------
        const taxonomy& taxa_;
        std::vector<ranked_lineage> lins_;
        bool outdated_;
    };


private:
    //---------------------------------------------------------------
    sketcher targetSketcher_;
    sketcher querySketcher_;
    std::uint64_t maxLocsPerFeature_;
    feature_store features_;
    std::vector<const taxon*> targets_;
    taxonomy taxa_;
    mutable ranked_lineages_cache ranksCache_;
    mutable ranked_lineages_of_targets targetLineages_;
    std::map<taxon_name,const taxon*> name2tax_;

    sketch_batch batch_;
    sketch_queue queue_;
    std::atomic_bool sketchingDone_;
    std::unique_ptr<std::future<void>> inserterThread_;
};




/*************************************************************************//**
 *
 * @brief pull some types into global namespace
 *
 *****************************************************************************/
using match_locations = database::match_locations;




/*************************************************************************//**
 *
 * @brief reads database from file
 *
 *****************************************************************************/
database
make_database(const std::string& filename,
              database::scope = database::scope::everything,
              info_level = info_level::moderate);



} // namespace mc

#endif

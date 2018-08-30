/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2018 André Müller (muellan@uni-mainz.de)
 *                       & Robin Kobus  (rkobus@uni-mainz.de)
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

#include "version.h"
#include "io_error.h"
#include "io_options.h"
#include "stat_combined.h"
#include "taxonomy.h"
#include "hash_multimap.h"
#include "dna_encoding.h"
#include "typename.h"


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
 *  location        (const taxon*, window_id) = "window within a target sequence"
 *
 *  full_lineage    path from root -> lowest taxon (vector of const taxon*);
 *                  may have arbitrary length
 *
 *  ranked_lineage  main rank taxon ids (array of const taxon*);
 *                  has fixed length; each index always refers to the same rank
 *
 * @tparam SequenceType  reference and query sequence type, usually std::string
 *
 * @tparam Sketcher      function object type, that maps reference sequence
 *                       windows (= sequence interval) to
 *                       sketches (= collection of features of the same type)
 *
 * @tparam FeatureHash   hash function for feature map
 *                       (default: identity for integer features,
 *                       so h(x) = x mod tablesize)
 *
 * @tparam TargetId      internal target (reference sequence) identification

 * @tparam WindowId      target window identification
 *
 * @tparam BucketSizeT   bucket (location list) size tracking
 *
 *****************************************************************************/
template<
    class SequenceType,
    class Sketcher,
    class FeatureHash = std::hash<typename Sketcher::feature_type>,
    class TargetId = std::uint16_t,
    class WindowId = std::uint16_t,
    class BucketSizeT = std::uint8_t
>
class sketch_database
{
public:
    //---------------------------------------------------------------
    using sequence = SequenceType;
    using sketcher = Sketcher;
    using feature_hash = FeatureHash;
    //-----------------------------------------------------
    using target_id = TargetId;
    using window_id = WindowId;
    using bucket_size_type = BucketSizeT;
    using match_count_type = std::uint16_t;
    //-----------------------------------------------------
    using taxon_id    = taxonomy::taxon_id;
    using taxon       = taxonomy::taxon;
    using taxon_rank  = taxonomy::rank;
    using taxon_name  = taxonomy::taxon_name;
    //-----------------------------------------------------
    using ranked_lineage = taxonomy::ranked_lineage;
    using full_lineage   = taxonomy::full_lineage;
    using file_source    = taxonomy::file_source;
    using taxon_iterator = taxonomy::const_iterator;
    using taxon_range    = taxonomy::const_range;


    //---------------------------------------------------------------
    enum class scope {
        everything, metadata_only
    };


    //-----------------------------------------------------
    class target_limit_exceeded_error : public std::runtime_error {
    public:
        target_limit_exceeded_error():
            std::runtime_error{"target count limit exceeded"}
        {}
    };


private:
    //use negative numbers for sequence level taxon ids
    static constexpr taxon_id
    taxon_id_of_target(target_id id) noexcept { return -taxon_id(id)-1; }


    //-----------------------------------------------------
    /** @brief internal location representation = (target index, window index)
     *         these are stored in the in-memory database and on disk
     */
    #pragma pack(push, 1)
    //avoid padding bits
    struct target_location
    {
        using target_id_t = target_id;
        using window_id_t = window_id;

        constexpr
        target_location(target_id g = 0, window_id w = 0) noexcept :
            tgt{g}, win{w}
        {}

        target_id tgt;
        window_id win;

        friend bool
        operator < (const target_location& a, const target_location& b) noexcept {
            if(a.tgt < b.tgt) return true;
            if(a.tgt > b.tgt) return false;
            return (a.win < b.win);
        }

        friend void read_binary(std::istream& is, target_location& p) {
            read_binary(is, p.tgt);
            read_binary(is, p.win);
        }
        friend void write_binary(std::ostream& os, const target_location& p) {
            write_binary(os, p.tgt);
            write_binary(os, p.win);
        }
    };
    //avoid padding bits
    #pragma pack(pop)


public:
    //-----------------------------------------------------
    using sketch  = typename sketcher::sketch_type;  //range of features
    using feature = typename sketch::value_type;


private:
    //-----------------------------------------------------
    //"heart of the database": maps features to target locations
    using feature_store = hash_multimap<feature,target_location, //key, value
                              feature_hash,               //key hasher
                              std::equal_to<feature>,     //key comparator
                              chunk_allocator<target_location>,  //value allocator
                              std::allocator<feature>,    //bucket+key allocator
                              bucket_size_type>;          //location list size


public:
    //---------------------------------------------------------------
    using feature_count_type = typename feature_store::size_type;


    //-----------------------------------------------------
    /** @brief external location represenation
     *         = (target sequence taxon pointer, window index)
     *         these are used to return match results
     */
    struct location
    {
        constexpr
        location(const taxon* t = nullptr, window_id w = 0) noexcept :
            tax{t}, win{w}
        {}

        const taxon* tax;
        window_id win;

        friend bool
        operator == (const location& a, const location& b) noexcept {
            return (a.tax == b.tax) && (a.win == b.win);
        }

        friend bool
        operator < (const location& a, const location& b) noexcept {
            if(a.tax < b.tax) return true;
            if(a.tax > b.tax) return false;
            return (a.win < b.win);
        }
    };

    using match_locations = std::vector<location>;


    //---------------------------------------------------------------
    explicit
    sketch_database(sketcher targetSketcher = sketcher{}) :
        sketch_database{targetSketcher, targetSketcher}
    {}
    //-----------------------------------------------------
    explicit
    sketch_database(sketcher targetSketcher, sketcher querySketcher) :
        targetSketcher_{std::move(targetSketcher)},
        querySketcher_{std::move(querySketcher)},
        targetWindowSize_(128),
        targetWindowStride_(128-15),
        queryWindowSize_(targetWindowSize_),
        queryWindowStride_(targetWindowStride_),
        maxLocsPerFeature_(max_supported_locations_per_feature()),
        features_{},
        targets_{},
        taxa_{},
        ranksCache_{taxa_, taxon_rank::Sequence},
        name2tax_{}
    {
        features_.max_load_factor(default_max_load_factor());
    }

    sketch_database(const sketch_database&) = delete;
    sketch_database(sketch_database&&)      = default;

    sketch_database& operator = (const sketch_database&) = delete;
    sketch_database& operator = (sketch_database&&)      = default;


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
    /** @brief set size of target windows that are fed to the target sketcher */
    void target_window_size(std::size_t s) {
        targetWindowSize_ = s;
    }
    //-----------------------------------------------------
    /** @return size of target windows that are fed to the target sketcher */
    size_t target_window_size() const noexcept {
        return targetWindowSize_;
    }

    //---------------------------------------------------------------
    /** @brief set target window stride for target sketching */
    void target_window_stride(std::size_t s) {
        if(s < 1) s = 1;
        targetWindowStride_ = s;
    }
    //-----------------------------------------------------
    /** @return target window stride for target sketching */
    size_t target_window_stride() const noexcept {
        return targetWindowStride_;
    }


    //---------------------------------------------------------------
    /** @brief set size of query windows that are fed to the query sketcher */
    void query_window_size(std::size_t s) {
        queryWindowSize_ = s;
    }
    //-----------------------------------------------------
    /** @return size of query windows that are fed to the query sketcher */
    size_t query_window_size() const noexcept {
        return queryWindowSize_;
    }
    //---------------------------------------------------------------
    /** @brief set query window stride for query sketching */
    void query_window_stride(std::size_t s) {
        if(s < 1) s = 1;
        queryWindowStride_ = s;
    }
    //-----------------------------------------------------
    /** @return query window stride for query sketching */
    size_t query_window_stride() const noexcept {
        return queryWindowStride_;
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
                        taxa.insert(ranksCache_[*targets_[loc.tgt]][int(r)]);
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

        return true;
    }


    //---------------------------------------------------------------
    taxon_id
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
        name2tax_.clear();
        features_.clear();
    }


    //-----------------------------------------------------
    /**
     * @brief very dangerous! clears feature map without memory deallocation
     */
    void clear_without_deallocation() {
        ranksCache_.clear();
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
        ranksCache_.update();
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
        if(taxa_.reset_parent(tax.id(), parentId)) {
            ranksCache_.update(tax);
        }
    }
    //-----------------------------------------------------
    void reset_parent(const taxon& tax, const taxon& parent) {
        if(taxa_.reset_parent(tax.id(), parent.id())) {
            ranksCache_.update(tax);
        }
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

    //-----------------------------------------------------
    const taxon*
    next_ranked_ancestor(const taxon* tax) const noexcept {
        return tax ? next_ranked_ancestor(*tax) : nullptr;
    }
    const taxon*
    next_ranked_ancestor(const taxon& tax) const noexcept {
        if(tax.rank() == taxon_rank::Sequence) {
            for(const taxon* a : ranksCache_[tax]) {
                if(a && a->rank() != taxon_rank::none && a->rank() > tax.rank())
                    return a;
            }
            return nullptr;
        }
        return taxa_.next_ranked_ancestor(tax);
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
    template<class Consumer>
    void
    for_each_match(const sequence& query, Consumer&& consume) const
    {
        if(query.empty()) return;

        using std::begin;
        using std::end;
        for_each_match(begin(query), end(query),
                             std::forward<Consumer>(consume));
    }

    //-----------------------------------------------------
    template<class InputIterator, class Consumer>
    void
    for_each_match(InputIterator queryBegin, InputIterator queryEnd,
                   Consumer&& consume) const
    {
        for_each_window(queryBegin, queryEnd, queryWindowSize_, queryWindowStride_,
            [this, &consume] (InputIterator b, InputIterator e) {
                auto sk = querySketcher_(b,e);
                for(auto f : sk) {
                    auto locs = features_.find(f);
                    if(locs != features_.end()) {
                        for(const auto& loc : *locs) {
                            consume(loc);
                        }
                    }
                }
            });
    }


    //---------------------------------------------------------------
    template<class InputIterator>
    void
    accumulate_matches(InputIterator queryBegin, InputIterator queryEnd,
                       match_locations& res) const
    {
        for_each_match(queryBegin, queryEnd,
            [this, &res] (const target_location& loc) {
                res.emplace_back(targets_[loc.tgt], loc.win);
            });
    }

    //---------------------------------------------------------------
    void
    accumulate_matches(const sequence& query,
                       match_locations& res) const
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
        return 0.85f;
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
        read_binary(is, targetWindowSize_);
        read_binary(is, targetWindowStride_);

        read_binary(is, querySketcher_);
        read_binary(is, queryWindowSize_);
        read_binary(is, queryWindowStride_);

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

        //ranked linage cache
        ranksCache_.update();

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
        write_binary(os, targetWindowSize_);
        write_binary(os, targetWindowStride_);

        write_binary(os, querySketcher_);
        write_binary(os, queryWindowSize_);
        write_binary(os, queryWindowStride_);

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
                for(target_location p : bucket) {
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


private:
    //---------------------------------------------------------------
    window_id add_all_window_sketches(const sequence& seq, target_id tgt)
    {
        using iter_t = typename sequence::const_iterator;

        window_id win = 0;
        for_each_window(seq, targetWindowSize_, targetWindowStride_,
            [&, this] (iter_t b, iter_t e) {
                auto sk = targetSketcher_(b,e);
                //insert features from sketch into database
                for(const auto& f : sk) {
                    auto it = features_.insert(f, target_location{tgt, win});
                    if(it->size() > maxLocsPerFeature_) {
                        features_.shrink(it, maxLocsPerFeature_);
                    }
                }
                ++win;
            });
        return win;
    }


    //---------------------------------------------------------------
    sketcher targetSketcher_;
    sketcher querySketcher_;
    std::uint64_t targetWindowSize_;
    std::uint64_t targetWindowStride_;
    std::uint64_t queryWindowSize_;
    std::uint64_t queryWindowStride_;
    std::uint64_t maxLocsPerFeature_;
    target_id targetCount_;
    feature_store features_;
    std::vector<const taxon*> targets_;
    taxonomy taxa_;
    ranked_lineages_cache ranksCache_;
    std::map<taxon_name,const taxon*> name2tax_;
};





/*************************************************************************//**
 *
 * @brief reads database from file
 *
 *****************************************************************************/
template<class Database>
Database
make_database(const std::string& filename,
              typename Database::scope what = Database::scope::everything,
              info_level info = info_level::moderate)
{
    if(filename.empty()) throw file_access_error{"No database name given"};

    Database db;

    const bool showInfo = info != info_level::silent;

    if(showInfo) {
        std::cerr << "Reading database from file '"
                  << filename << "' ... " << std::flush;
    }
    try {
        db.read(filename, what);
        if(showInfo) std::cerr << "done." << std::endl;
    }
    catch(const file_access_error& e) {
        std::cerr << "FAIL" << std::endl;
        throw file_access_error{"Could not read database file '" + filename + "'"};
    }

    return db;
}



/*************************************************************************//**
 *
 * @brief prints database properties to stdout
 *
 *****************************************************************************/
template<class S, class K, class H, class G, class W, class L>
void print_static_properties(const sketch_database<S,K,H,G,W,L>& db)
{
    using db_t = sketch_database<S,K,H,G,W,L>;
    using target_id = typename db_t::target_id;
    using window_id = typename db_t::window_id;
    using feature_t = typename db_t::feature;
    using bkt_sz_t  = typename db_t::bucket_size_type;

    std::cout
        << "------------------------------------------------\n"
        << "MetaCache version    " << MC_VERSION_STRING << " (" << MC_VERSION << ")\n"
        << "database verion      " << MC_DB_VERSION << '\n'
        << "------------------------------------------------\n"
        << "sequence type        " << type_name<typename db_t::sequence>() << '\n'
        << "target id type       " << type_name<target_id>() << " " << (sizeof(target_id)*CHAR_BIT) << " bits\n"
        << "target limit         " << std::uint64_t(db.max_target_count()) << '\n'
        << "------------------------------------------------\n"
        << "window id type       " << type_name<window_id>() << " " << (sizeof(window_id)*CHAR_BIT) << " bits\n"
        << "window limit         " << std::uint64_t(db.max_windows_per_target()) << '\n'
        << "window length        " << db.target_window_size() << '\n'
        << "window stride        " << db.target_window_stride() << '\n'
        << "------------------------------------------------\n"
        << "sketcher type        " << type_name<typename db_t::sketcher>() << '\n'
        << "feature type         " << type_name<feature_t>() << " " << (sizeof(feature_t)*CHAR_BIT) << " bits\n"
        << "feature hash         " << type_name<typename db_t::feature_hash>() << '\n'
        << "kmer size            " << std::uint64_t(db.target_sketcher().kmer_size()) << '\n'
        << "kmer limit           " << std::uint64_t(db.target_sketcher().max_kmer_size()) << '\n'
        << "sketch size          " << db.target_sketcher().sketch_size() << '\n'
        << "------------------------------------------------\n"
        << "bucket size type     " << type_name<bkt_sz_t>() << " " << (sizeof(bkt_sz_t)*CHAR_BIT) << " bits\n"
        << "max. locations       " << std::uint64_t(db.max_locations_per_feature()) << '\n'
        << "location limit       " << std::uint64_t(db.max_supported_locations_per_feature()) << '\n'
        << "------------------------------------------------"
        << std::endl;
}




/*************************************************************************//**
 *
 * @brief prints database properties to stdout
 *
 *****************************************************************************/
template<class S, class K, class H, class G, class W, class L>
void print_content_properties(const sketch_database<S,K,H,G,W,L>& db)
{
    if(db.target_count() > 0) {

        std::uint64_t numRankedTargets = 0;
        for(const auto& t : db.target_taxa()) {
            if(t.has_parent()) ++numRankedTargets;
        }

        std::cout
        << "targets              " << db.target_count() << '\n'
        << "ranked targets       " << numRankedTargets << '\n'
        << "taxa in tree         " << db.non_target_taxon_count() << '\n';
    }

    if(db.feature_count() > 0) {
        auto lss = db.location_list_size_statistics();

        std::cout
        << "------------------------------------------------\n"
        << "buckets              " << db.bucket_count() << '\n'
        << "bucket size          " << "max: " << lss.max()
                                   << " mean: " << lss.mean()
                                   << " +/- " << lss.stddev()
                                   << " <> " << lss.skewness() << '\n'
        << "features             " << db.feature_count() << '\n'
        << "dead features        " << db.dead_feature_count() << '\n'
        << "locations            " << db.location_count() << '\n';
    }
    std::cout
        << "------------------------------------------------\n";
}


} // namespace mc

#endif

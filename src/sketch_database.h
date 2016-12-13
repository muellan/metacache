/*****************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016 André Müller (muellan@uni-mainz.de)
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
#include <typeinfo>

#include "version.h"
#include "io_error.h"
#include "stat_moments.h"
#include "taxonomy.h"
#include "hash_multimap.h"


namespace mc {

/*****************************************************************************
 *
 * @brief  maps 'features' (e.g. hash values obtained by min-hashing)
 *         to 'targets' = positions in reference targets
 *
 * @details
 *   terminology
 *   target:      reference sequence whose sketches are stored in the DB
 *
 *   query:       sequence (usually short reads) that shall be matched against
 *                the reference targets
 *
 *   target_id:   numeric database identifier of reference targets
 *                in the range [0,n-1] where n is the number of targets in the DB
 *
 *   window_id:   window index (starting with 0) within a reference target
 *
 *   sequence_id: alphanumeric sequence identifier (e.g. an NCBI accession)
 *
 *   taxon_id:    numeric taxon identifier
 *
 * @tparam
 *   SequenceType:  type of reference and query sequence, usually std::string
 *
 *   Sketcher:      function object type, that maps reference sequence
 *                  windows (= sequence interval) to
 *                  sketches (= collection of features of the same C++ type)
 *   TargetId:      type for reference sequence identification;
 *
 *                  may have heavy impact on memory footprint of database
 *
 *****************************************************************************/
template<
    class SequenceType,
    class Sketcher,
    class TargetId = std::uint16_t,
    class WindowId = std::uint16_t
>
class sketch_database
{
public:
    //---------------------------------------------------------------
    using sequence = SequenceType;
    using sketcher = Sketcher;
    //-----------------------------------------------------
    using target_id   = TargetId;
    using window_id   = WindowId;
    using sequence_id = std::string;
    //-----------------------------------------------------
    using taxon_id   = taxonomy::taxon_id;
    using taxon      = taxonomy::taxon;
    using taxon_rank = taxonomy::rank;
    //-----------------------------------------------------
    using ranked_lineage = taxonomy::ranked_lineage;
    using full_lineage   = taxonomy::full_lineage;


    //---------------------------------------------------------------
    struct sequence_origin {
        std::string filename = "";
        std::uint32_t index = 0;  //if file contains more than one sequence
    };

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
    //-----------------------------------------------------
    /**
     * @brief target meta-information
     */
    class target_property
    {
    public:
        using taxon_id = sketch_database::taxon_id;

        explicit
        target_property(std::string identifier = "",
                        taxon_id tax = 0,
                        sequence_origin origin = sequence_origin{})
        :
            id(std::move(identifier)),
            taxonId(tax),
            ranks{}, lineage{},
            origin{std::move(origin)}
        {}

        friend void
        read_binary(std::istream& is, target_property& p) {
            read_binary(is, p.id);
            read_binary(is, p.taxonId);
            read_binary(is, p.ranks);
            read_binary(is, p.lineage);
            read_binary(is, p.origin.filename);
            read_binary(is, p.origin.index);
        }

        friend void
        write_binary(std::ostream& os, const target_property& p) {
            write_binary(os, p.id);
            write_binary(os, p.taxonId);
            write_binary(os, p.ranks);
            write_binary(os, p.lineage);
            write_binary(os, p.origin.filename);
            write_binary(os, p.origin.index);
        }

        std::string id;
        taxon_id taxonId;
        ranked_lineage ranks;
        full_lineage lineage;
        sequence_origin origin;
    };


public:
    //-----------------------------------------------------
    using sketch  = typename sketcher::result_type;
    using feature = typename sketch::value_type;

    //-----------------------------------------------------
    struct location
    {
        constexpr
        location(target_id g = 0, window_id w = 0) noexcept :
            tgt{g}, win{w}
        {}

        target_id tgt;
        window_id win;

        friend bool
        operator < (const location& a, const location& b) noexcept {
            if(a.tgt < b.tgt) return true;
            if(a.tgt > b.tgt) return false;
            return (a.win < b.win);
        }

        friend void read_binary(std::istream& is, location& p) {
            read_binary(is, p.tgt);
            read_binary(is, p.win);
        }
        friend void write_binary(std::ostream& os, const location& p) {
            write_binary(os, p.tgt);
            write_binary(os, p.win);
        }
    };


private:
    //-----------------------------------------------------
    using feature_store = hash_multimap<feature,location>;


public:
    //-------------------------------------------------------------------
    using bucket_size_type  = typename feature_store::bucket_size_type;
    using match_result = std::map<location,std::uint16_t>; //features per target


    //---------------------------------------------------------------
    explicit
    sketch_database(sketcher targetSketcher = sketcher{}) :
        targetSketcher_{std::move(targetSketcher)},
        querySketcher_{targetSketcher_},
        targetWindowSize_(128),
        targetWindowStride_(128-15),
        queryWindowSize_(targetWindowSize_),
        queryWindowStride_(targetWindowStride_),
        maxLocsPerFeature_(max_supported_locations_per_feature()),
        nextTargetId_(0),
        targets_{},
        features_{},
        sid2gid_{},
        taxa_{}
    {
        features_.max_load_factor(0.8);
    }
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
        nextTargetId_(0),
        targets_{},
        features_{},
        sid2gid_{},
        taxa_{}
    {
        features_.max_load_factor(0.8);
    }

    sketch_database(const sketch_database&) = delete;
    sketch_database(sketch_database&&)      = default;

    sketch_database& operator = (const sketch_database&) = delete;
    sketch_database& operator = (sketch_database&&)      = default;


    //---------------------------------------------------------------
    const sketcher&
    target_sketcher() const noexcept {
        return targetSketcher_;
    }
    //-----------------------------------------------------
    const sketcher&
    query_sketcher() const noexcept {
        return querySketcher_;
    }
    //-----------------------------------------------------
    void
    query_sketcher(const sketcher& s) {
        querySketcher_ = s;
    }
    void
    query_sketcher(sketcher&& s) {
        querySketcher_ = std::move(s);
    }


    //---------------------------------------------------------------
    void target_window_size(std::size_t s) {
        targetWindowSize_ = s;
    }
    //-----------------------------------------------------
    size_t target_window_size() const noexcept {
        return targetWindowSize_;
    }
    //---------------------------------------------------------------
    void target_window_stride(std::size_t s) {
        if(s < 1) s = 1;
        targetWindowStride_ = s;
    }
    //-----------------------------------------------------
    size_t target_window_stride() const noexcept {
        return targetWindowStride_;
    }


    //---------------------------------------------------------------
    void query_window_size(std::size_t s) {
        queryWindowSize_ = s;
    }
    //-----------------------------------------------------
    size_t query_window_size() const noexcept {
        return queryWindowSize_;
    }
    //---------------------------------------------------------------
    void query_window_stride(std::size_t s) {
        if(s < 1) s = 1;
        queryWindowStride_ = s;
    }
    //-----------------------------------------------------
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
    void erase_features_with_more_locations_than(bucket_size_type n)
    {
        for(auto i = features_.begin(), e = features_.end(); i != e; ++i) {
            if(i->size() > n) features_.clear(i);
        }
    }


    //---------------------------------------------------------------
    void erase_non_unique_features_on_rank(taxon_rank r)
    {
        if(taxa_.empty()) {
            throw std::runtime_error{"no taxonomy available!"};
        }

        for(auto i = features_.begin(), e = features_.end(); i != e; ++i) {
            //check if any two locations refer to different taxa on rank r
            if(!i->empty()) {
                auto j = i->begin();
                taxon_id taxid = targets_[j->tgt].ranks[int(r)];
                for(++j; j != i->end(); ++j) {
                    if(taxid != targets_[j->tgt].ranks[int(r)]) {
                        features_.clear(i);
                        break;
                    }
                }
            }
        }
    }


    //---------------------------------------------------------------
    bool add_target(const sequence& seq, sequence_id sid,
                    taxon_id taxid = 0,
                    const sequence_origin& origin = sequence_origin{})
    {
        using std::begin;
        using std::end;

        //reached hard limit for number of targets
        if(nextTargetId_ >= max_target_count()) {
            throw target_limit_exceeded_error{};
        }

        using iter_t = typename sequence::const_iterator;

        if(seq.empty()) return false;

        //don't allow non-unique sequence ids
        auto it = sid2gid_.find(sid);
        if(it != sid2gid_.end()) return false;

        sid2gid_.insert({sid, nextTargetId_});

        targets_.emplace_back(std::move(sid), taxid, origin);
        if(taxid > 0 && !taxa_.empty()) {
            rank_target(nextTargetId_, taxid);
        }

        window_id win = 0;
        for_each_window(seq, targetWindowSize_, targetWindowStride_,
            [this, &win] (iter_t b, iter_t e) {
                auto sk = targetSketcher_(b,e);
                //insert features from sketch into database
                for(const auto& f : sk) {
                    auto it = features_.insert(f, location{nextTargetId_, win});
                    if(it->size() > maxLocsPerFeature_) {
                        features_.shrink(it, maxLocsPerFeature_);
                    }
                }
                ++win;
            });

        ++nextTargetId_;

        return true;
    }

    //---------------------------------------------------------------
    bool
    is_valid(target_id tid) const noexcept {
        return tid < nextTargetId_;
    }
    static constexpr target_id
    invalid_target_id() noexcept {
        return max_target_count();
    }
    std::uint64_t
    target_count() const noexcept {
        return targets_.size();
    }
    static constexpr std::uint64_t
    max_target_count() noexcept {
        return std::numeric_limits<target_id>::max();
    }

    //-----------------------------------------------------
    bool empty() const noexcept {
        return features_.empty();
    }


    //---------------------------------------------------------------
    void clear() {
        targets_.clear();
        sid2gid_.clear();
        features_.clear();
    }

    //-----------------------------------------------------
    /**
     * @brief very dangerous! clears feature map without memory deallocation
     */
    void clear_without_deallocation() {
        targets_.clear();
        sid2gid_.clear();
        features_.clear_without_deallocation();
    }


    //---------------------------------------------------------------
    const sequence_id&
    sequence_id_of_target(target_id tid) const noexcept {
        return targets_[tid].id;
    }
    //-----------------------------------------------------
    target_id
    target_id_of_sequence(const sequence_id& sid) const noexcept {
        auto it = sid2gid_.find(sid);
        return (it != sid2gid_.end()) ? it->second : nextTargetId_;
    }


    //---------------------------------------------------------------
    void apply_taxonomy(const taxonomy& tax) {
        taxa_ = tax;
        for(auto& g : targets_) update_lineages(g);
    }
    //-----------------------------------------------------
    void apply_taxonomy(taxonomy&& tax) {
        taxa_ = std::move(tax);
        for(auto& gp : targets_) update_lineages(gp);
    }


    //---------------------------------------------------------------
    std::uint64_t
    taxon_count() const noexcept {
        return taxa_.taxon_count();
    }
    //-----------------------------------------------------
    const taxon&
    taxon_with_id(taxon_id id) const noexcept {
        return taxa_[id];
    }


    //---------------------------------------------------------------
    void rank_target(target_id tid, taxon_id taxid)
    {
        targets_[tid].taxonId = taxid;
        update_lineages(targets_[tid]);
    }
    //-----------------------------------------------------
    taxon_id
    taxon_id_of_target(target_id id) const noexcept {
        return taxa_[targets_[id].taxonId].id;
    }
    //-----------------------------------------------------
    const taxon&
    taxon_of_target(target_id id) const noexcept {
        return taxa_[targets_[id].taxonId];
    }

    //-----------------------------------------------------
    const sequence_origin&
    origin_of_target(target_id id) const noexcept {
        return targets_[id].origin;
    }

    //-----------------------------------------------------
    full_lineage
    lineage(const taxon& tax) const noexcept {
        return taxa_.lineage(tax.id);
    }
    //-----------------------------------------------------
    const full_lineage&
    lineage_of_target(target_id id) const noexcept {
        return targets_[id].lineage;
    }

    //-----------------------------------------------------
    ranked_lineage
    ranks(const taxon& tax) const noexcept {
        return taxa_.ranks(tax.id);
    }
    //-----------------------------------------------------
    const ranked_lineage&
    ranks_of_target(target_id id) const noexcept {
        return targets_[id].ranks;
    }

    //---------------------------------------------------------------
    const taxon&
    lca(const taxon& ta, const taxon& tb) const
    {
        return taxa_.lca(ta,tb);
    }
    //---------------------------------------------------------------
    const taxon&
    lca_of_targets(target_id ta, target_id tb) const
    {
        return taxa_.lca(lineage_of_target(ta), lineage_of_target(tb));
    }

    //---------------------------------------------------------------
    const taxon&
    ranked_lca(const taxon& ta, const taxon& tb) const
    {
        return taxa_.ranked_lca(ta,tb);
    }
    //---------------------------------------------------------------
    const taxon&
    ranked_lca_of_targets(target_id ta, target_id tb) const
    {
        return taxa_.lca(ranks_of_target(ta), ranks_of_target(tb));
    }

    //---------------------------------------------------------------
    taxon_rank
    lowest_rank(const taxon& tax) const {
        return taxa_.lowest_rank(tax.id);
    }
    //-----------------------------------------------------
    taxon_rank
    lowest_rank_of_taxon(taxon_id id) const {
        return taxa_.lowest_rank(id);
    }

    //---------------------------------------------------------------
    /**
     * @return number of times a taxon is covered by any target in the DB
     */
    int covering(const taxon& t) const {
        return covering_of_taxon(t.id);
    }
    int covering_of_taxon(taxon_id id) const {
        int cover = 0;

        for(const auto& g : targets_) {
            for(taxon_id taxid : g.lineage) {
                if(taxid == id) ++cover;
            }
        }
        return cover;
    }

    //-----------------------------------------------------
    /**
     * @return true, if taxon is covered by any target in the DB
     */
    bool covers(const taxon& t) const {
        return covers_taxon(t.id);
    }
    bool covers_taxon(taxon_id id) const {
        for(const auto& g : targets_) {
            for(taxon_id taxid : g.lineage) {
                if(taxid == id) return true;
            }
        }
        return false;
    }


    //---------------------------------------------------------------
    template<class Consumer>
    void
    for_each_match(const sequence& query, Consumer&& consume) const
    {
        using iter_t = typename sequence::const_iterator;

        if(query.empty()) return;

        for_each_window(query, queryWindowSize_, queryWindowStride_,
            [this, &consume] (iter_t b, iter_t e) {
                auto sk = querySketcher_(b,e);
                for(auto f : sk) {
                    auto sit = features_.find(f);
                    if(sit != features_.end()) {
                        for(const auto& pos : *sit) {
                            consume(pos);
                        }
                    }
                }
            });
    }


    //---------------------------------------------------------------
    match_result
    matches(const sequence& query) const
    {
        auto res = match_result{};
        accumulate_matches(query, res);
        return res;
    }
    //---------------------------------------------------------------
    void
    accumulate_matches(const sequence& query,
                       match_result& res) const
    {
        for_each_match(query,
            [this, &res] (const location& t) {
                auto it = res.find(t);
                if(it != res.end()) {
                    ++(it->second);
                }
                else {
                    res.insert(it, {t, 1});
                }
            });
    }


    //-------------------------------------------------------------------
    void max_load_factor(float lf) {
        features_.max_load_factor(lf);
    }
    //-----------------------------------------------------
    float max_load_factor() const noexcept {
        return features_.max_load_factor();
    }


    //-------------------------------------------------------------------
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
        uint64_t dbVer = 0;
        read_binary(is, dbVer);

        if(MC_DB_VERSION != dbVer) {
            throw file_read_error{
                "Database " + filename +
                " is incompatible with incompatible version of MetaCache" };
        }

        //data type info
        std::uint8_t featureSize = 0;
        read_binary(is, featureSize);

        std::uint8_t targetSize = 0;
        read_binary(is, targetSize);

        std::uint8_t windowSize = 0;
        read_binary(is, windowSize);

        if( (sizeof(feature) != featureSize) ||
            (sizeof(target_id) != targetSize) ||
            (sizeof(window_id) != windowSize) )
        {
            throw file_read_error{
                "Database " + filename +
                " is incompatible with this variant of MetaCache" +
                " due to different data type widths"};
        }

        clear();

        //sketching parameters
        read_binary(is, targetSketcher_);
        read_binary(is, targetWindowSize_);
        read_binary(is, targetWindowStride_);

        read_binary(is, querySketcher_);
        read_binary(is, queryWindowSize_);
        read_binary(is, queryWindowStride_);

        read_binary(is, maxLocsPerFeature_);

        //taxon metadata
        read_binary(is, taxa_);

        read_binary(is, nextTargetId_);
        read_binary(is, targets_);
        if(targets_.size() < 1) return;

        sid2gid_.clear();
        for(target_id i = 0; i < target_id(targets_.size()); ++i) {
            sid2gid_.insert({targets_[i].id, i});
        }

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
        std::ofstream os{filename, std::ios::out | std::ios::binary};

        if(!os.good()) {
            throw file_access_error{"can't open file " + filename};
        }

        //database version info
        write_binary(os, std::uint64_t(MC_DB_VERSION));

        //data type info
        write_binary(os, std::uint8_t(sizeof(feature)));
        write_binary(os, std::uint8_t(sizeof(target_id)));
        write_binary(os, std::uint8_t(sizeof(window_id)));

        //sketching parameters
        write_binary(os, targetSketcher_);
        write_binary(os, targetWindowSize_);
        write_binary(os, targetWindowStride_);

        write_binary(os, querySketcher_);
        write_binary(os, queryWindowSize_);
        write_binary(os, queryWindowStride_);

        write_binary(os, maxLocsPerFeature_);

        //taxon metadata
        write_binary(os, taxa_);

        //target metainformation
        write_binary(os, nextTargetId_);
        write_binary(os, targets_);

        //hash table
        write_binary(os, features_);
    }


    //---------------------------------------------------------------
    std::uint64_t bucket_count() const noexcept {
        return features_.bucket_count();
    }
    //---------------------------------------------------------------
    std::uint64_t feature_count() const noexcept {
        return features_.non_empty_bucket_count();
    }
    //---------------------------------------------------------------
    std::uint64_t location_count() const noexcept {
        return features_.value_count();
    }


    //---------------------------------------------------------------
    variance_accumulator<double>
    bucket_size_statistics() const {
        auto priSize = variance_accumulator<double>{};

        for(const auto& bucket : features_) {
            priSize += bucket.size();
        }

        return priSize;
    }


    //---------------------------------------------------------------
    void print_feature_map(std::ostream& os) const {
        for(const auto& bucket : features_) {
            if(!bucket.empty()) {
                os << bucket.key() << " -> ";
                for(location p : bucket) {
                    os << '(' << p.tgt << ',' << p.win << ')';
                }
                os << '\n';
            }
        }
    }


private:

    //---------------------------------------------------------------
    void update_lineages(target_property& gp)
    {
        if(gp.taxonId > 0) {
            gp.lineage = taxa_.lineage(gp.taxonId);
            gp.ranks = taxa_.ranks(gp.taxonId);
        }
    }


    //---------------------------------------------------------------
    sketcher targetSketcher_;
    sketcher querySketcher_;
    std::uint64_t targetWindowSize_;
    std::uint64_t targetWindowStride_;
    std::uint64_t queryWindowSize_;
    std::uint64_t queryWindowStride_;
    std::uint64_t maxLocsPerFeature_;
    target_id nextTargetId_;
    std::vector<target_property> targets_;
    feature_store features_;
    std::map<sequence_id,target_id> sid2gid_;
    taxonomy taxa_;
};






/*****************************************************************************
 *
 * @brief
 *
 *****************************************************************************/
template<class ValueT>
struct index_range
{
    using value_type = ValueT;

    constexpr
    index_range() noexcept :
        beg(0), end(0)
    {}
    constexpr
    index_range(value_type first, value_type last) noexcept :
        beg{first}, end{last}
    {}

    constexpr value_type size()  const noexcept { return beg - end; }
    constexpr value_type empty() const noexcept { return beg == end; }

    value_type beg;
    value_type end;
};



/*****************************************************************************
*
* @brief processes a database hit list and
*        stores the top (in terms of accumulated hits) 'maxNo' contiguous
*        window ranges
*
*****************************************************************************/
template<int maxNo, class TgtId>
class matches_in_contiguous_window_range_top
{
    static_assert(maxNo > 1, "no must be > 1");


public:

    //---------------------------------------------------------------
    using tid_t = TgtId;
    using hit_t = std::uint_least64_t;
    using win_t = std::uint_least64_t;
    using window_range = index_range<win_t>;

    //---------------------------------------------------------------
    static constexpr int max_count() noexcept { return maxNo; }

    static constexpr tid_t invalid_tgt() noexcept {
        return std::numeric_limits<tid_t>::max();
    }


    /****************************************************************
     * @pre matches must be sorted by target (first) and window (second)
     */
    template<class MatchResult>
    matches_in_contiguous_window_range_top(
        const MatchResult& matches, win_t numWindows = 3)
    :
        tgt_{}, hits_{}, pos_{}
    {
        using std::begin;
        using std::end;

        for(int i = 0; i < maxNo; ++i) {
            tgt_[i] = invalid_tgt();
            hits_[i] = 0;
        }

        tid_t tgt = invalid_tgt();
        hit_t hits = 0;
        hit_t maxHits = 0;
        hit_t win = 0;
        win_t maxWinBeg = 0;
        win_t maxWinEnd = 0;

        //check hits per query sequence
        auto fst = begin(matches);
        auto lst = fst;
        while(lst != end(matches)) {
            //look for neighboring windows with the highest total hit count
            //as long as we are in the same target and the windows are in a
            //contiguous range
            if(lst->first.tgt == tgt) {
                //add new hits to the right
                hits += lst->second;
                //subtract hits to the left that fall out of range
                while(fst != lst &&
                     (lst->first.win - fst->first.win) >= numWindows)
                {
                    hits -= fst->second;
                    //move left side of range
                    ++fst;
                    win = fst->first.win;
                }
                //track best of the local sub-ranges
                if(hits > maxHits) {
                    maxHits = hits;
                    maxWinBeg = win;
                    maxWinEnd = win + distance(fst,lst);
                }
            }
            else {
                //reset to new target
                win = lst->first.win;
                tgt = lst->first.tgt;
                hits = lst->second;
                maxHits = hits;
                maxWinBeg = win;
                maxWinEnd = win;
                fst = lst;
            }
            //keep track of 'maxNo' largest
            //TODO binary search for large maxNo?
            for(int i = 0; i < maxNo; ++i) {
                if(maxHits >= hits_[i]) {
                    //shift to the right
                    for(int j = maxNo-1; j > i; --j) {
                        hits_[j] = hits_[j-1];
                        tgt_[j] = tgt_[j-1];
                        pos_[j] = pos_[j-1];
                    }
                    //set hits & associated sequence (position)
                    hits_[i] = maxHits;
                    tgt_[i] = tgt;
                    pos_[i].beg = maxWinBeg;
                    pos_[i].end = maxWinEnd;
                    break;
                }
            }
            ++lst;
        }
    }

    int count() const noexcept {
        for(int i = 0; i < maxNo; ++i) {
            if(hits_[i] < 1) return i;
        }
        return maxNo;
    }

    tid_t target_id(int rank)   const noexcept { return tgt_[rank];  }
    hit_t hits(int rank) const noexcept { return hits_[rank]; }

    hit_t total_hits() const noexcept {
        int h = 0;
        for(int i = 0; i < maxNo; ++i) h += hits_[i];
        return h;
    }

    const window_range&
    window(int rank) const noexcept { return pos_[rank]; }

private:
    tid_t tgt_[maxNo];
    hit_t hits_[maxNo];
    window_range pos_[maxNo];
};




/*****************************************************************************
 *
 * @brief reads database from file
 *
 *****************************************************************************/
template<class Database>
Database
make_database(const std::string& filename,
              typename Database::scope what = Database::scope::everything)
{
    Database db;

    std::cout << "Reading database from file '"
              << filename << "' ... " << std::flush;
    try {
        db.read(filename, what);
        std::cout << "done." << std::endl;

        return db;
    }
    catch(const file_access_error& e) {
        std::cout << "FAIL" << std::endl;
        throw file_access_error{"Could not read database file '" + filename + "'"};
    }

    return db;
}

//-------------------------------------------------------------------
template<class Database>
Database
make_database_metadata_only(const std::string& filename)
{
    return make_database<Database>(filename, Database::scope::metadata_only);
}




/*****************************************************************************
 *
 * @brief writes database to file
 *
 *****************************************************************************/
template<class S, class K, class G, class W>
void
write_database(const sketch_database<S,K,G,W>& db, const std::string& filename)
{
    std::cout << "Writing database to file'"
              << filename << "' ... " << std::flush;
    try {
        db.write(filename);
        std::cout << "done." << std::endl;

    }
    catch(const file_access_error&) {
        std::cout << "FAIL" << std::endl;
        throw file_access_error{"Could not write database file '" + filename + "'"};
    }
}




/*****************************************************************************
 *
 * @brief prints database properties to stdout
 *
 *****************************************************************************/
template<class S, class K, class G, class W>
void print_config(const sketch_database<S,K,G,W>& db)
{
    using db_t = sketch_database<S,K,G,W>;
    using target_id = typename db_t::target_id;
    using window_id = typename db_t::window_id;
    using feature_t = typename db_t::feature;

    int numRankedTargets = 0;
    for(target_id i = 0; i < db.target_count(); ++i) {
        if(!db.taxon_of_target(i).none()) ++numRankedTargets;
    }

    std::cout
        << "database format  " << MC_DB_VERSION << '\n'
        << "sequence type    " << typeid(typename db_t::sequence).name() << '\n'
        << "target id type:  " << typeid(target_id).name() << " " << (sizeof(target_id)*8) << " bits\n"
        << "window id type:  " << typeid(window_id).name() << " " << (sizeof(window_id)*8) << " bits\n"
        << "window length:   " << db.target_window_size() << '\n'
        << "window stride:   " << db.target_window_stride() << '\n'
        << "sketcher type    " << typeid(typename db_t::sketcher).name() << '\n'
        << "feature type:    " << typeid(feature_t).name() << " " << (sizeof(feature_t)*8) << " bits\n"
        << "kmer size:       " << std::uint64_t(db.target_sketcher().kmer_size()) << '\n'
        << "sketch size:     " << db.target_sketcher().sketch_size() << '\n'
        << "location limit:  " << std::uint64_t(db.max_locations_per_feature()) << '\n';
}


//-------------------------------------------------------------------
template<class S, class K, class G, class W>
void print_data_properties(const sketch_database<S,K,G,W>& db)
{
    using db_t = sketch_database<S,K,G,W>;
    using target_id = typename db_t::target_id;

    std::uint64_t numRankedTargets = 0;
    for(target_id i = 0; i < db.target_count(); ++i) {
        if(!db.taxon_of_target(i).none()) ++numRankedTargets;
    }

    std::cout
        << "targets:         " << db.target_count() << '\n'
        << "ranked targets:  " << numRankedTargets << '\n'
        << "taxa in tree:    " << db.taxon_count() << '\n';
}


//-------------------------------------------------------------------
template<class S, class K, class G, class W>
void print_statistics(const sketch_database<S,K,G,W>& db)
{
    auto hbs = db.bucket_size_statistics();

    std::cout
        << "buckets:         " << db.bucket_count() << '\n'
        << "bucket size:     " << hbs.mean() << " +/- " << hbs.stddev() << '\n'
        << "features:        " << db.feature_count() << '\n'
        << "locations:       " << db.location_count() << '\n';
}


} // namespace mc


#undef MC_DB_VERSION

#endif

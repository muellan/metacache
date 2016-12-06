/*****************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * version 1.0
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

#include "io_error.h"
#include "stat_moments.h"
#include "taxonomy.h"
#include "hash_multimap.h"


namespace mc {

//helps to distinguish potentially incompatible database versions
#define METACACHE_DB_VERSION 20160912


/*****************************************************************************
 *
 * @brief  maps 'features' (e.g. hash values obtained by min-hashing)
 *         to 'targets' = positions in reference genomes
 *
 * @details
 *   terminology
 *   genome:      reference sequence whose sketches are stored in the DB
 *
 *   query:       sequence (usually short reads) that shall be matched against
 *                the reference genomes
 *
 *   genome_id:   numeric database identifier of reference genomes
 *                in the range [0,n-1] where n is the number of genomes in the DB
 *
 *   window_id:   window index (starting with 0) within a reference genome
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
 *   GenomeId:      type for reference sequence identification;
 *
 *                  may have heavy impact on memory footprint of database
 *
 *****************************************************************************/
template<
    class SequenceType,
    class Sketcher,
    class GenomeId = std::uint16_t
>
class sketch_database
{
public:
    //---------------------------------------------------------------
    using sequence = SequenceType;
    using sketcher = Sketcher;
    //-----------------------------------------------------
    using genome_id   = GenomeId;
    using window_id   = std::uint16_t;
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
    class genome_limit_exceeded_error : public std::runtime_error {
    public:
        genome_limit_exceeded_error():
            std::runtime_error{"genome count limit exceeded"}
        {}
    };


private:
    //-----------------------------------------------------
    /**
     * @brief genome meta-information
     */
    class genome_property
    {
    public:
        using taxon_id = sketch_database::taxon_id;

        explicit
        genome_property(std::string identifier = "",
                        taxon_id tax = 0,
                        sequence_origin origin = sequence_origin{})
        :
            id(std::move(identifier)),
            taxonId(tax),
            ranks{}, lineage{},
            origin{std::move(origin)}
        {}

        friend void
        read_binary(std::istream& is, genome_property& p) {
            read_binary(is, p.id);
            read_binary(is, p.taxonId);
            read_binary(is, p.ranks);
            read_binary(is, p.lineage);
            read_binary(is, p.origin.filename);
            read_binary(is, p.origin.index);
        }

        friend void
        write_binary(std::ostream& os, const genome_property& p) {
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
    struct target
    {
        constexpr
        target(genome_id g = 0, window_id w = 0) noexcept :
            gid{g}, win{w}
        {}

        genome_id gid;
        window_id win;

        friend bool
        operator < (const target& a, const target& b) noexcept {
            if(a.gid < b.gid) return true;
            if(a.gid > b.gid) return false;
            return (a.win < b.win);
        }

        friend void read_binary(std::istream& is, target& p) {
            read_binary(is, p.gid);
            read_binary(is, p.win);
        }
        friend void write_binary(std::ostream& os, const target& p) {
            write_binary(os, p.gid);
            write_binary(os, p.win);
        }
    };


private:
    //-----------------------------------------------------
    using feature_store = hash_multimap<feature,target>;


public:
    //-------------------------------------------------------------------
    using bucket_size_type  = typename feature_store::bucket_size_type;
    using match_result = std::map<target,std::uint16_t>; //features per target


    //---------------------------------------------------------------
    explicit
    sketch_database(sketcher genomeSketcher = sketcher{}) :
        genomeSketcher_{std::move(genomeSketcher)},
        querySketcher_{genomeSketcher_},
        genomeWindowSize_(128),
        genomeWindowStride_(128-15),
        queryWindowSize_(genomeWindowSize_),
        queryWindowStride_(genomeWindowStride_),
        maxRefsPerFeature_(features_.max_bucket_size()-1),
        nextGenomeId_(0),
        genomes_{},
        features_{},
        sid2gid_{},
        taxa_{}
    {}
    //-----------------------------------------------------
    explicit
    sketch_database(sketcher genomeSketcher, sketcher querySketcher) :
        genomeSketcher_{std::move(genomeSketcher)},
        querySketcher_{std::move(querySketcher)},
        genomeWindowSize_(128),
        genomeWindowStride_(128-15),
        queryWindowSize_(genomeWindowSize_),
        queryWindowStride_(genomeWindowStride_),
        maxRefsPerFeature_(features_.max_bucket_size()-1),
        nextGenomeId_(0),
        genomes_{},
        features_{},
        sid2gid_{},
        taxa_{}
    {}

    sketch_database(const sketch_database&) = delete;
    sketch_database(sketch_database&&)      = default;

    sketch_database& operator = (const sketch_database&) = delete;
    sketch_database& operator = (sketch_database&&)      = default;


    //---------------------------------------------------------------
    const sketcher&
    genome_sketcher() const noexcept {
        return genomeSketcher_;
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
    void genome_window_size(std::size_t s) {
        genomeWindowSize_ = s;
    }
    //-----------------------------------------------------
    size_t genome_window_size() const noexcept {
        return genomeWindowSize_;
    }
    //---------------------------------------------------------------
    void genome_window_stride(std::size_t s) {
        if(s < 1) s = 1;
        genomeWindowStride_ = s;
    }
    //-----------------------------------------------------
    size_t genome_window_stride() const noexcept {
        return genomeWindowStride_;
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
    void max_genomes_per_feature(bucket_size_type n)
    {
        if(n < 1) n = 1;
        if(n >= features_.max_bucket_size()) {
            n = features_.max_bucket_size() - 1;
        }
        else if(n < maxRefsPerFeature_) {
            for(auto i = features_.begin(), e = features_.end(); i != e; ++i) {
                if(i->size() > n) features_.shrink(i, n);
            }
        }
        maxRefsPerFeature_ = n;
    }
    //-----------------------------------------------------
    bucket_size_type max_genomes_per_feature() const noexcept {
        return maxRefsPerFeature_;
    }

    //-----------------------------------------------------
    void erase_features_with_more_genomes_than(bucket_size_type n)
    {
        for(auto i = features_.begin(), e = features_.end(); i != e; ++i) {
            if(i->size() > n) features_.erase(i);
        }
    }


    //---------------------------------------------------------------
    bool add_genome(const sequence& seq,
                    sequence_id sid,
                    taxon_id taxonId = 0,
                    const sequence_origin& origin = sequence_origin{})
    {
        using std::begin;
        using std::end;

        //reached hard limit for number of genomes
        if(nextGenomeId_ >= max_genome_count()) {
            throw genome_limit_exceeded_error{};
        }

        using iter_t = typename sequence::const_iterator;

        if(seq.empty()) return false;

        //don't allow non-unique sequence ids
        auto it = sid2gid_.find(sid);
        if(it != sid2gid_.end()) return false;

        sid2gid_.insert({sid, nextGenomeId_});

        genomes_.emplace_back(std::move(sid), taxonId, origin);
        if(taxonId > 0 && !taxa_.empty()) {
            rank_genome(nextGenomeId_, taxonId);
        }

        window_id win = 0;
        for_each_window(seq, genomeWindowSize_, genomeWindowStride_,
            [this, &win] (iter_t b, iter_t e) {
                auto sk = genomeSketcher_(b,e);
                //insert features from sketch into database
                for(const auto& f : sk) {
                    auto it = features_.insert(f, target{nextGenomeId_, win});
                    if(it->size() > maxRefsPerFeature_) {
                        features_.shrink(it, maxRefsPerFeature_);
                    }
                }
                ++win;
            });

        ++nextGenomeId_;

        return true;
    }

    //---------------------------------------------------------------
    bool
    is_valid(genome_id gid) const noexcept {
        return gid < nextGenomeId_;
    }
    static constexpr genome_id
    invalid_genome_id() noexcept {
        return max_genome_count();
    }
    std::uint64_t
    genome_count() const noexcept {
        return genomes_.size();
    }
    static constexpr std::uint64_t
    max_genome_count() noexcept {
        return std::numeric_limits<genome_id>::max();
    }

    //-----------------------------------------------------
    bool empty() const noexcept {
        return features_.empty();
    }


    //---------------------------------------------------------------
    void clear() {
        genomes_.clear();
        features_.clear();
        sid2gid_.clear();
    }


    //---------------------------------------------------------------
    const sequence_id&
    sequence_id_of_genome(genome_id gid) const noexcept {
        return genomes_[gid].id;
    }
    //-----------------------------------------------------
    genome_id
    genome_id_of_sequence(const sequence_id& sid) const noexcept {
        auto it = sid2gid_.find(sid);
        return (it != sid2gid_.end()) ? it->second : nextGenomeId_;
    }


    //---------------------------------------------------------------
    void apply_taxonomy(const taxonomy& tax) {
        taxa_ = tax;
        for(auto& g : genomes_) update_lineages(g);
    }
    //-----------------------------------------------------
    void apply_taxonomy(taxonomy&& tax) {
        taxa_ = std::move(tax);
        for(auto& gp : genomes_) update_lineages(gp);
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
    void rank_genome(genome_id gid, taxon_id taxid)
    {
        genomes_[gid].taxonId = taxid;
        update_lineages(genomes_[gid]);
    }
    //-----------------------------------------------------
    taxon_id
    taxon_id_of_genome(genome_id id) const noexcept {
        return taxa_[genomes_[id].taxonId].id;
    }
    //-----------------------------------------------------
    const taxon&
    taxon_of_genome(genome_id id) const noexcept {
        return taxa_[genomes_[id].taxonId];
    }

    //-----------------------------------------------------
    const sequence_origin&
    origin_of_genome(genome_id id) const noexcept {
        return genomes_[id].origin;
    }

    //-----------------------------------------------------
    full_lineage
    lineage(const taxon& tax) const noexcept {
        return taxa_.lineage(tax.id);
    }
    //-----------------------------------------------------
    const full_lineage&
    lineage_of_genome(genome_id id) const noexcept {
        return genomes_[id].lineage;
    }

    //-----------------------------------------------------
    ranked_lineage
    ranks(const taxon& tax) const noexcept {
        return taxa_.ranks(tax.id);
    }
    //-----------------------------------------------------
    const ranked_lineage&
    ranks_of_genome(genome_id id) const noexcept {
        return genomes_[id].ranks;
    }

    //---------------------------------------------------------------
    const taxon&
    lca(const taxon& ta, const taxon& tb) const
    {
        return taxa_.lca(ta,tb);
    }
    //---------------------------------------------------------------
    const taxon&
    lca_of_genomes(genome_id ga, genome_id gb) const
    {
        return taxa_[taxonomy::lca_id(lineage_of_genome(ga),
                                      lineage_of_genome(gb))];
    }

    //---------------------------------------------------------------
    const taxon&
    ranked_lca(const taxon& ta, const taxon& tb) const
    {
        return taxa_.ranked_lca(ta,tb);
    }
    //---------------------------------------------------------------
    const taxon&
    ranked_lca_of_genomes(genome_id ga, genome_id gb) const
    {
        return taxa_[taxonomy::ranked_lca_id(ranks_of_genome(ga),
                                             ranks_of_genome(gb))];
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
     * @return number of times a taxon is covered by any genome in the DB
     */
    int covering(const taxon& t) const {
        return covering_of_taxon(t.id);
    }
    int covering_of_taxon(taxon_id id) const {
        int cover = 0;

        for(const auto& g : genomes_) {
            for(taxon_id tid : g.lineage) {
                if(tid == id) ++cover;
            }
        }
        return cover;
    }

    //-----------------------------------------------------
    /**
     * @return true, if taxon is covered by any genome in the DB
     */
    bool covers(const taxon& t) const {
        return covers_taxon(t.id);
    }
    bool covers_taxon(taxon_id id) const {
        for(const auto& g : genomes_) {
            for(taxon_id tid : g.lineage) {
                if(tid == id) return true;
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
            [this, &res] (const target& t) {
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

        if(METACACHE_DB_VERSION != dbVer) {
            throw file_read_error{
                "Database " + filename +
                " was built with an incompatible version of MetaCache\n"
            };
        }

        clear();

        //sketching parameters
        read_binary(is, genomeSketcher_);
        read_binary(is, genomeWindowSize_);
        read_binary(is, genomeWindowStride_);

        read_binary(is, querySketcher_);
        read_binary(is, queryWindowSize_);
        read_binary(is, queryWindowStride_);

        read_binary(is, maxRefsPerFeature_);

        //taxon metadata
        read_binary(is, taxa_);

        read_binary(is, nextGenomeId_);
        read_binary(is, genomes_);
        if(genomes_.size() < 1) return;

        sid2gid_.clear();
        for(genome_id i = 0; i < genome_id(genomes_.size()); ++i) {
            sid2gid_.insert({genomes_[i].id, i});
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
        write_binary(os, std::uint64_t(METACACHE_DB_VERSION));

        //sketching parameters
        write_binary(os, genomeSketcher_);
        write_binary(os, genomeWindowSize_);
        write_binary(os, genomeWindowStride_);

        write_binary(os, querySketcher_);
        write_binary(os, queryWindowSize_);
        write_binary(os, queryWindowStride_);

        write_binary(os, maxRefsPerFeature_);

        //taxon metadata
        write_binary(os, taxa_);

        //genome metainformation
        write_binary(os, nextGenomeId_);
        write_binary(os, genomes_);

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
    std::uint64_t reference_count() const noexcept {
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
            os << bucket.key() << " -> ";
            for(target p : bucket) {
                os << '(' << p.gid << ',' << p.win << ')';
            }
            os << '\n';
        }
    }


private:

    //---------------------------------------------------------------
    void update_lineages(genome_property& gp)
    {
        if(gp.taxonId > 0) {
            gp.lineage = taxa_.lineage(gp.taxonId);
            gp.ranks = taxa_.ranks(gp.taxonId);
        }
    }


    //---------------------------------------------------------------
    sketcher genomeSketcher_;
    sketcher querySketcher_;
    std::uint64_t genomeWindowSize_;
    std::uint64_t genomeWindowStride_;
    std::uint64_t queryWindowSize_;
    std::uint64_t queryWindowStride_;
    std::uint64_t maxRefsPerFeature_;
    genome_id nextGenomeId_;
    std::vector<genome_property> genomes_;
    feature_store features_;
    std::map<sequence_id,genome_id> sid2gid_;
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
template<int maxNo, class GidT = std::uint_least64_t>
class matches_in_contiguous_window_range_top
{
    static_assert(maxNo > 1, "no must be > 1");


public:

    //---------------------------------------------------------------
    using window_range = index_range<int>;
    using hit_t = std::uint_least64_t;
    using gid_t = GidT;

    //---------------------------------------------------------------
    static constexpr int max_count() noexcept { return maxNo; }

    static constexpr gid_t invalid_gid() noexcept {
        return std::numeric_limits<gid_t>::max();
    }


    /****************************************************************
     * @pre matches must be sorted by genome (first) and window (second)
     */
    template<class MatchResult>
    matches_in_contiguous_window_range_top(
        const MatchResult& matches, int numWindows = 3)
    :
        gid_{}, hits_{}, pos_{}
    {
        using std::begin;
        using std::end;

        for(int i = 0; i < maxNo; ++i) {
            gid_[i] = invalid_gid();
            hits_[i] = 0;
        }

        gid_t gid = invalid_gid();
        hit_t hits = 0;
        hit_t maxHits = 0;
        hit_t win = 0;
        hit_t maxWinBeg = 0;
        hit_t maxWinEnd = 0;

        //check hits per query sequence
        auto fst = begin(matches);
        auto lst = fst;
        while(lst != end(matches)) {
            //look for neighboring windows with the highest total hit count
            //as long as we are in the same genome and the windows are in a
            //contiguous range
            if(lst->first.gid == gid) {
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
                //reset to new genome
                win = lst->first.win;
                gid = lst->first.gid;
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
                        gid_[j] = gid_[j-1];
                        pos_[j] = pos_[j-1];
                    }
                    //set hits & associated sequence (position)
                    hits_[i] = maxHits;
                    gid_[i] = gid;
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

    gid_t genome_id(int rank)   const noexcept { return gid_[rank];  }
    hit_t hits(int rank) const noexcept { return hits_[rank]; }

    hit_t total_hits() const noexcept {
        int h = 0;
        for(int i = 0; i < maxNo; ++i) h += hits_[i];
        return h;
    }

    const window_range&
    window(int rank) const noexcept { return pos_[rank]; }

private:
    gid_t gid_[maxNo];
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
template<class S, class K, class G>
void
write_database(const sketch_database<S,K,G>& db, const std::string& filename)
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
template<class S, class K, class G>
void print_statistics(const sketch_database<S,K,G>& db)
{
    using genome_id = typename sketch_database<S,K,G>::genome_id;
    int numRankedGenomes = 0;
    for(genome_id i = 0; i < db.genome_count(); ++i) {
        if(!db.taxon_of_genome(i).none()) ++numRankedGenomes;
    }

    std::cout << "sequences:        " << db.genome_count() << '\n'
              << "ranked sequences: " << numRankedGenomes << '\n'
              << "window length:    " << db.genome_window_size() << '\n'
              << "window stride:    " << db.genome_window_stride() << '\n'
              << "kmer size:        " << int(db.genome_sketcher().kmer_size()) << '\n'
              << "sketch size:      " << db.genome_sketcher().sketch_size() << '\n'
              << "taxa in tree:     " << db.taxon_count() << '\n';

    auto hbs = db.bucket_size_statistics();

    std::cout << "buckets:          " << db.bucket_count() << '\n'
              << "bucket size:      " << hbs.mean() << " +/- " << hbs.stddev() << '\n'
              << "features:         " << db.feature_count() << '\n'
              << "references:       " << db.reference_count() << '\n'
              << std::endl;
}


} // namespace mc


#undef METACACHE_DB_VERSION

#endif

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

#include "generic.h"
#include "io_error.h"
#include "stat_moments.h"
#include "taxonomy.h"
#include "hash_multimap.h"


namespace mc {

#define METACACHE_DB_VERSION 20160606


/*****************************************************************************
 *
 * @brief
 *
 * @details
 *   terminology:
 *   genome      reference sequence whose sketches are stored in the DB
 *
 *   query       sequence (usually short reads) that shall be classified
 *               (ideally a slightly modified version of a genome's sub-sequence)
 *
 *   genome_id   numeric database identifier of reference sequence
 *               in the range [0,n-1] where n is the number of genomes in the DB
 *
 *   sequence_id alphanumeric sequence identifier
 *   taxon_id    numeric taxon identifier
 *
 *****************************************************************************/
template<
    class SequenceType,
    class Sketcher
>
class sketch_database
{
public:
    //---------------------------------------------------------------
    using sequence_type = SequenceType;
    using sketcher_type = Sketcher;
    //-----------------------------------------------------
    using genome_id   = std::uint16_t;
    using window_id   = std::uint16_t;
    using sequence_id = std::string;
    //-----------------------------------------------------
    using taxon_id   = taxonomy::taxon_id;
    using taxon_type = taxonomy::taxon;
    using taxon_rank = taxonomy::taxon_rank;
    //-----------------------------------------------------
    using ranked_lineage_type = taxonomy::ranked_lineage_type;
    using full_lineage_type   = taxonomy::full_lineage_type;


    //---------------------------------------------------------------
    struct sequence_origin {
        std::string filename = "";
        std::size_t index = 0;  //if file contains more than one sequence
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
            rankedLineage{}, fullLineage{},
            origin{std::move(origin)}
        {}

        void read(std::istream& is) {
            //id
            std::string::size_type l = 0;
            is.read(reinterpret_cast<char*>(&l), sizeof(l));
            id.clear();
            id.resize(l);
            is.read(reinterpret_cast<char*>(&(*id.begin())), l);
            //tax id
            is.read(reinterpret_cast<char*>(&taxonId), sizeof(taxonId));
            //ranked lineage
            for(auto& x : rankedLineage) {
                is.read(reinterpret_cast<char*>(&x), sizeof(x));
            }
            //full lineage
            std::size_t linSize = 0;
            is.read(reinterpret_cast<char*>(&linSize), sizeof(linSize));
            fullLineage.clear();
            fullLineage.resize(linSize);
            for(auto& x : fullLineage) {
                is.read(reinterpret_cast<char*>(&x), sizeof(x));
            }
            //sequence filename
            l = 0;
            is.read(reinterpret_cast<char*>(&l), sizeof(l));
            origin.filename.clear();
            origin.filename.resize(l);
            is.read(reinterpret_cast<char*>(&(*origin.filename.begin())), l);
            //index in file
            is.read(reinterpret_cast<char*>(&origin.index), sizeof(origin.index));
        }

        void write(std::ostream& os) const {
            //id
            auto l = id.size();
            os.write(reinterpret_cast<const char*>(&l), sizeof(l));
            os.write(id.c_str(), id.size());
            //tax id
            os.write(reinterpret_cast<const char*>(&taxonId), sizeof(taxonId));
            //ranked lineage
            for(auto x : rankedLineage) {
                os.write(reinterpret_cast<char*>(&x), sizeof(x));
            }
            //full lineage
            std::size_t linSize = fullLineage.size();
            os.write(reinterpret_cast<char*>(&linSize), sizeof(linSize));
            for(auto x : fullLineage) {
                os.write(reinterpret_cast<char*>(&x), sizeof(x));
            }
            //sequence filename
            l = origin.filename.size();
            os.write(reinterpret_cast<const char*>(&l), sizeof(l));
            os.write(origin.filename.c_str(), origin.filename.size());
            //index in file
            os.write(reinterpret_cast<const char*>(&origin.index), sizeof(origin.index));
        }

        std::string id;
        taxon_id taxonId;
        ranked_lineage_type rankedLineage;
        full_lineage_type fullLineage;
        sequence_origin origin;
    };


public:
    //-----------------------------------------------------
    using sketch_type  = typename sketcher_type::result_type;
    using sketch_value = typename sketch_type::value_type;

    //-----------------------------------------------------
    struct reference_pos {
        reference_pos() = default;

        explicit
        reference_pos(genome_id s, window_id w) noexcept :
            gid{s}, win{w}
        {}

        genome_id gid = 0;
        window_id win = 0;

        friend bool
        operator < (const reference_pos& a, const reference_pos& b) noexcept {
            if(a.gid < b.gid) return true;
            if(a.gid > b.gid) return false;
            return (a.win < b.win);
        }
    };


private:
    //-----------------------------------------------------
    using store_t = hash_multimap<sketch_value,reference_pos>;


public:
    //-------------------------------------------------------------------
    using bucket_size_type  = typename store_t::bucket_size_type;
    using match_result_type = std::map<reference_pos,std::uint16_t>;


    //---------------------------------------------------------------
    explicit
    sketch_database(sketcher_type sketcher = sketcher_type{}) :
        refSketcher_{std::move(sketcher)},
        querySketcher_{refSketcher_},
        genomeWindowSize_(128),
        genomeWindowStride_(128 - querySketcher_.kmer_size()),
        queryWindowSize_(genomeWindowSize_),
        queryWindowStride_(genomeWindowStride_),
        maxRefsPerSketchVal_(1024),
        numSeq_(0),
        sketchVals_{},
        genomes_{}, noGenome_{},
        sid2gid_{},
        taxa_{}
    { }


    //---------------------------------------------------------------
    size_t kmer_size() const noexcept {
        return refSketcher_.kmer_size();
    }


    //---------------------------------------------------------------
    void genome_window_size(std::size_t s) {
        if(s < 2*kmer_size()) s = 2*kmer_size();
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
    size_t genome_sketch_size() const noexcept {
        return refSketcher_.sketch_size();
    }


    //---------------------------------------------------------------
    void query_window_size(std::size_t s) {
        if(s < 2*kmer_size()) s = 2*kmer_size();
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
    size_t query_sketch_size() const noexcept {
        return querySketcher_.sketch_size();
    }
    //---------------------------------------------------------------
    void query_sketch_size(size_t s) {
        querySketcher_.sketch_size(s);
    }


    //---------------------------------------------------------------
    void max_genomes_per_sketch_value(bucket_size_type n)
    {
        if(n < 1) n = 1;
        maxRefsPerSketchVal_ = n;

        for(auto i = sketchVals_.begin(), e = sketchVals_.end(); i != e; ++i) {
            if(i->size() > n) sketchVals_.limit(i, n);
        }
    }
    //-----------------------------------------------------
    bucket_size_type max_genomes_per_sketch_value() const noexcept {
        return maxRefsPerSketchVal_;
    }

    //-----------------------------------------------------
    void erase_sketch_values_with_more_genomes_than(bucket_size_type n)
    {
        for(auto i = sketchVals_.begin(), e = sketchVals_.end(); i != e; ++i) {
            if(i->size() > n) sketchVals_.erase(i);
        }
    }


    //---------------------------------------------------------------
    bool add_genome(const sequence_type& seq,
                    sequence_id sid,
                    taxon_id taxonId = 0,
                    const sequence_origin& origin = sequence_origin{})
    {
        using std::begin;
        using std::end;

        using iter_t = typename sequence_type::const_iterator;

        if(seq.size() < kmer_size()) return false;

        //don't allow non-unique sequence ids
        auto it = sid2gid_.find(sid);
        if(it != sid2gid_.end()) return false;

        sid2gid_.insert({sid, numSeq_});

        genomes_.emplace_back(std::move(sid), taxonId, origin);
        if(taxonId > 0 && !taxa_.empty()) {
            rank_genome(numSeq_, taxonId);
        }

        window_id win = 0;
        for_each_window(genomeWindowSize_, genomeWindowStride_, seq,
            [this, &win] (iter_t b, iter_t e) {
                auto sketch = refSketcher_(b,e);
                //insert values from sketch into database
                for(auto s : sketch) {
                    auto it = sketchVals_.insert(s, reference_pos{numSeq_, win});
                    if(it->size() > maxRefsPerSketchVal_) {
                        sketchVals_.limit(it, maxRefsPerSketchVal_);
                    }
                }
                ++win;
            });

        ++numSeq_;

        return true;
    }

    //---------------------------------------------------------------
    std::size_t genome_count() const noexcept {
        return numSeq_;
    }

    //-----------------------------------------------------
    bool empty() const noexcept {
        return sketchVals_.empty();
    }


    //---------------------------------------------------------------
    void clear() {
        genomes_.clear();
        sketchVals_.clear();
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
        return (it != sid2gid_.end()) ? it->second : numSeq_;
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
    std::size_t
    taxon_count() const noexcept {
        return taxa_.taxon_count();
    }
    //-----------------------------------------------------
    const taxon_type&
    taxon(taxon_id id) const noexcept {
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
    const taxon_type&
    taxon_of_genome(genome_id id) const noexcept {
        return taxa_[genomes_[id].taxonId];
    }

    //-----------------------------------------------------
    const sequence_origin&
    origin_of_genome(genome_id id) const noexcept {
        return genomes_[id].origin;
    }

    //-----------------------------------------------------
    full_lineage_type
    full_lineage(const taxon_type& tax) const noexcept {
        return taxa_.full_lineage(tax.id);
    }
    //-----------------------------------------------------
    const full_lineage_type&
    full_lineage_of_genome(genome_id id) const noexcept {
        return genomes_[id].fullLineage;
    }

    //-----------------------------------------------------
    ranked_lineage_type
    ranked_lineage(const taxon_type& tax) const noexcept {
        return taxa_.ranked_lineage(tax.id);
    }
    //-----------------------------------------------------
    const ranked_lineage_type&
    ranked_lineage_of_genome(genome_id id) const noexcept {
        return genomes_[id].rankedLineage;
    }

    //---------------------------------------------------------------
    const taxon_type&
    lca(const taxon_type& ta, const taxon_type& tb) const
    {
        return taxa_.lca(ta,tb);
    }
    //---------------------------------------------------------------
    const taxon_type&
    lca_of_genomes(genome_id ga, genome_id gb) const
    {
        return taxa_[taxonomy::lca_id(full_lineage_of_genome(ga),
                                      full_lineage_of_genome(gb))];
    }

    //---------------------------------------------------------------
    const taxon_type&
    ranked_lca(const taxon_type& ta, const taxon_type& tb) const
    {
        return taxa_.ranked_lca(ta,tb);
    }
    //---------------------------------------------------------------
    const taxon_type&
    ranked_lca_of_genomes(genome_id ga, genome_id gb) const
    {
        return taxa_[taxonomy::ranked_lca_id(ranked_lineage_of_genome(ga),
                                             ranked_lineage_of_genome(gb))];
    }

    //---------------------------------------------------------------
    taxon_rank
    lowest_rank(const taxon_type& tax) const {
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
    int covering(const taxon_type& t) const {
        return covering_of_taxon(t.id);
    }
    int covering_of_taxon(taxon_id id) const {
        int cover = 0;

        for(const auto& g : genomes_) {
            for(taxon_id tid : g.fullLineage) {
                if(tid == id) ++cover;
            }
        }
        return cover;
    }

    //-----------------------------------------------------
    /**
     * @return true, if taxon is covered by any genome in the DB
     */
    bool covers(const taxon_type& t) const {
        return covers_taxon(t.id);
    }
    bool covers_taxon(taxon_id id) const {
        for(const auto& g : genomes_) {
            for(taxon_id tid : g.fullLineage) {
                if(tid == id) return true;
            }
        }
        return false;
    }


    //---------------------------------------------------------------
    template<class Consumer>
    void
    for_each_match(const sequence_type& query, Consumer&& consume) const
    {
        using iter_t = typename sequence_type::const_iterator;

        if(query.empty()) return;

        for_each_window(queryWindowSize_, queryWindowStride_, query,
            [this, &consume] (iter_t b, iter_t e) {
                auto sketch = querySketcher_(b,e);
                for(auto s : sketch) {
                    auto sit = sketchVals_.find(s);
                    if(sit != sketchVals_.end()) {
                        for(const auto& pos : *sit) {
                            consume(pos);
                        }
                    }
                }
            });
    }


    //---------------------------------------------------------------
    match_result_type
    matches(const sequence_type& query) const
    {
        auto res = match_result_type{};
        accumulate_matches(query, res);
        return res;
    }
    //---------------------------------------------------------------
    void
    accumulate_matches(const sequence_type& query,
                       match_result_type& res) const
    {
        for_each_match(query,
            [this, &res] (const reference_pos& pos) {
                auto it = res.find(pos);
                if(it != res.end()) {
                    ++(it->second);
                }
                else {
                    res.insert(it, {pos, 1});
                }
            });
    }


    //-------------------------------------------------------------------
    void max_load_factor(float lf) {
        sketchVals_.max_load_factor(lf);
    }
    //-----------------------------------------------------
    float max_load_factor() const noexcept {
        return sketchVals_.max_load_factor();
    }


    //-------------------------------------------------------------------
    /**
     * @brief   read database from binary file
     * @details Note that the map is not just de-serialized but
     *          rebuilt by inserting individual keys and values
     *          This should make DB files more robust against changes in the
     *          internal mapping structure.
     */
    void read(const std::string& filename)
    {
        std::ifstream is{filename, std::ios::in | std::ios::binary};

        if(!is.good()) {
            throw file_access_error{"can't open file " + filename};
        }

        //database version info
        std::uint32_t dbVer = 0;
        is.read(reinterpret_cast<char*>(&dbVer), sizeof(dbVer));

        if(METACACHE_DB_VERSION < dbVer) {
            throw file_read_error{
                "Database " + filename +
                " was built with a newer version of MetaCache\n" +
                "that is incompatible with this version!"
            };
        }

        clear();

        //read basic properties
        auto kmersize = refSketcher_.kmer_size();
        is.read(reinterpret_cast<char*>(&kmersize), sizeof(kmersize));
        refSketcher_.kmer_size(kmersize);
        querySketcher_.kmer_size(kmersize);

        auto refsketchsize = refSketcher_.sketch_size();
        is.read(reinterpret_cast<char*>(&refsketchsize), sizeof(refsketchsize));
        refSketcher_.sketch_size(refsketchsize);
        is.read(reinterpret_cast<char*>(&genomeWindowSize_), sizeof(genomeWindowSize_));
        is.read(reinterpret_cast<char*>(&genomeWindowStride_), sizeof(genomeWindowStride_));

        auto querysketchsize = refSketcher_.sketch_size();
        is.read(reinterpret_cast<char*>(&querysketchsize), sizeof(querysketchsize));
        querySketcher_.sketch_size(querysketchsize);
        is.read(reinterpret_cast<char*>(&queryWindowSize_), sizeof(queryWindowSize_));
        is.read(reinterpret_cast<char*>(&queryWindowStride_), sizeof(queryWindowStride_));

        is.read(reinterpret_cast<char*>(&maxRefsPerSketchVal_), sizeof(maxRefsPerSketchVal_));

        taxa_.read(is);

        //number of sequences
        is.read(reinterpret_cast<char*>(&numSeq_), sizeof(numSeq_));
        if(numSeq_ < 1) return;

        //read genome metainformation
        genomes_.clear();
        genomes_.reserve(numSeq_);
        for(size_t i = 0; i < numSeq_; ++i) {
            genome_property p;
            p.read(is);
            genomes_.push_back(std::move(p));
        }

        sid2gid_.clear();
        for(genome_id i = 0; i < genome_id(genomes_.size()); ++i) {
            sid2gid_.insert({genomes_[i].id, i});
        }

        //read sketch values (keys) and their origin references
        size_t nkeys = 0;
        is.read(reinterpret_cast<char*>(&nkeys), sizeof(nkeys));
        if(nkeys < 1) return;
        size_t nvalues = 0;
        is.read(reinterpret_cast<char*>(&nvalues), sizeof(nvalues));
        if(nvalues < 1) return;

        sketchVals_.clear();
        sketchVals_.rehash( typename store_t::size_type(
            10 + (1/sketchVals_.max_load_factor() * nkeys)));

        for(size_t i = 0; i < nkeys; ++i) {
            typename store_t::key_type key;
            typename store_t::bucket_size_type nvals;
            is.read(reinterpret_cast<char*>(&key), sizeof(key));
            is.read(reinterpret_cast<char*>(&nvals), sizeof(nvals));

            auto vals = std::unique_ptr<reference_pos[]>{new reference_pos[nvals]};

            for(auto v = vals.get(), e = v+nvals; v < e; ++v) {
                is.read(reinterpret_cast<char*>(&(v->gid)), sizeof(v->gid));
                is.read(reinterpret_cast<char*>(&(v->win)), sizeof(v->win));
            }

            sketchVals_.replace(key, nvals, std::move(vals));
        }

//        std::cout << '\n'
//            << "keys:        " << sketchVals_.key_count() << " <> " << nkeys << '\n'
//            << "values:      " << sketchVals_.value_count() << " <> " << nvalues << '\n'
//            << "buckets:     " << sketchVals_.bucket_count() << '\n'
//            << "bucket size: " << (8*sizeof(sketchVals_.bucket(0))) << '\n'
//            << "key size:    " << (8*sizeof(sketchVals_.bucket(0).key())) << '\n'
//            << "value size:  " << (8*sizeof(sketchVals_.bucket(0)[0])) << '\n'
//            << "gid size:    " << (8*sizeof(sketchVals_.bucket(0)[0].gid)) << '\n'
//            << "win size:    " << (8*sizeof(sketchVals_.bucket(0)[0].win)) << '\n'
//            << "addr0      : " << &sketchVals_.bucket(0)[0] << '\n'
//            << "addr1      : " << &sketchVals_.bucket(0)[1] << '\n'
//            << std::flush;
//        char c;
//        std::cin >> c;
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
        std::uint32_t dbVer = METACACHE_DB_VERSION;
        os.write(reinterpret_cast<const char*>(&dbVer), sizeof(dbVer));

        //write basic properties
        auto kmersize = refSketcher_.kmer_size();
        os.write(reinterpret_cast<const char*>(&kmersize), sizeof(kmersize));

        auto refsketchsize = refSketcher_.sketch_size();
        os.write(reinterpret_cast<const char*>(&refsketchsize), sizeof(refsketchsize));
        os.write(reinterpret_cast<const char*>(&genomeWindowSize_), sizeof(genomeWindowSize_));
        os.write(reinterpret_cast<const char*>(&genomeWindowStride_), sizeof(genomeWindowStride_));

        auto querysketchsize = querySketcher_.sketch_size();
        os.write(reinterpret_cast<const char*>(&querysketchsize), sizeof(querysketchsize));
        os.write(reinterpret_cast<const char*>(&queryWindowSize_), sizeof(queryWindowSize_));
        os.write(reinterpret_cast<const char*>(&queryWindowStride_), sizeof(queryWindowStride_));

        os.write(reinterpret_cast<const char*>(&maxRefsPerSketchVal_), sizeof(maxRefsPerSketchVal_));

        taxa_.write(os);

        //number of genomes
        os.write(reinterpret_cast<const char*>(&numSeq_), sizeof(numSeq_));

        //write genome metainformation
        for(const auto& a : genomes_) {
            a.write(os);
        }

        //write sketch values (keys) and their origin references
        auto nkeys = sketchVals_.key_count();
        os.write(reinterpret_cast<char*>(&nkeys), sizeof(nkeys));
        auto nvalues = sketchVals_.value_count();
        os.write(reinterpret_cast<char*>(&nvalues), sizeof(nvalues));

        for(const auto& bucket : sketchVals_) {
            //store non-empty buckets only
            if(!bucket.empty()) {
                auto key = bucket.key();
                os.write(reinterpret_cast<const char*>(&key), sizeof(key));

                auto nvals = bucket.size();
                os.write(reinterpret_cast<const char*>(&nvals), sizeof(nvals));

                for(const auto& ref : bucket) {
                    os.write(reinterpret_cast<const char*>(&ref.gid), sizeof(ref.gid));
                    os.write(reinterpret_cast<const char*>(&ref.win), sizeof(ref.win));
                }
            }
        }
    }

    //---------------------------------------------------------------
    size_t bucket_count() const noexcept {
        return sketchVals_.bucket_count();
    }


    //---------------------------------------------------------------
    variance_accumulator<double>
    bucket_sizes() const {
        auto priSize = variance_accumulator<double>{};

        for(const auto& bucket : sketchVals_) {
            priSize += bucket.size();
        }

        return priSize;
    }


    //---------------------------------------------------------------
    timer hash_table_query_benchmark(const std::vector<sketch_value>& keys,
                                     std::vector<reference_pos>& result)
    {
        timer time;
        time.start();

        for(auto k : keys) {
            auto it = sketchVals_.find(k);
            if(it != sketchVals_.end()) {
                for(const auto& v : *it) {
                    result.push_back(v);
                }
            }
        }

        time.stop();

        return time;
    }


private:
    //---------------------------------------------------------------
    void update_lineages(genome_property& gp)
    {
        if(gp.taxonId > 0) {
            gp.fullLineage = taxa_.full_lineage(gp.taxonId);
            gp.rankedLineage = taxa_.ranked_lineage(gp.taxonId);
        }
    }

    
    //---------------------------------------------------------------
    sketcher_type refSketcher_;
    sketcher_type querySketcher_;
    size_t genomeWindowSize_;
    size_t genomeWindowStride_;
    size_t queryWindowSize_;
    size_t queryWindowStride_;
    size_t maxRefsPerSketchVal_;
    genome_id numSeq_;
    store_t sketchVals_;
    std::vector<genome_property> genomes_;
    genome_property noGenome_;
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
* @brief
*
*****************************************************************************/
template<int maxNo>
class matches_in_contiguous_window_range_top
{
public:
    static_assert(maxNo > 1, "no must be > 1");

    using window_range = index_range<int>;

    static constexpr int max_count() noexcept { return maxNo; }

    /**
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
            gid_[i] = -1;
            hits_[i] = 0;
        }

        int gid = -1;
        int hits = 0;
        int win = -1;

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
            }
            else {
                //reset to new genome
                win = lst->first.win;
                gid = lst->first.gid;
                hits = lst->second;
                fst = lst;
            }
            //keep track of 'maxNo' largest
            //TODO binary search?
            for(int i = 0; i < maxNo; ++i) {
                if(hits >= hits_[i]) {
                    //shift to the right
                    for(int j = maxNo-1; j > i; --j) {
                        hits_[j] = hits_[j-1];
                        gid_[j] = gid_[j-1];
                        pos_[j] = pos_[j-1];
                    }
                    //set hits & associated sequence (position)
                    hits_[i] = hits;
                    gid_[i] = gid;
                    pos_[i].beg = win;
                    pos_[i].end = win + distance(fst,lst) + 1;
                    break;
                }
            }
            ++lst;
        }
    }

    int count() const noexcept {
        for(int i = 0; i < maxNo; ++i) {
            if(gid_[i] < 0 || hits_[i] < 1) return i;
        }
        return maxNo;
    }

    int genome_id(int rank)   const noexcept { return gid_[rank];  }
    int hits(int rank) const noexcept { return hits_[rank]; }

    int total_hits() const noexcept {
        int h = 0;
        for(int i = 0; i < maxNo; ++i) h += hits_[i];
        return h;
    }

    const window_range&
    window(int rank) const noexcept { return pos_[rank]; }

private:
    int gid_[maxNo];
    int hits_[maxNo];
    window_range pos_[maxNo];
};




/*****************************************************************************
 *
 * @brief reads database from file
 *
 *****************************************************************************/
template<class Database>
Database
make_database(const std::string& filename)
{
    Database db;

    std::cout << "Reading database from file '"
              << filename << "' ... " << std::flush;
    try {
        db.read(filename);
        std::cout << "done." << std::endl;

        return db;
    }
    catch(const file_access_error& e) {
        std::cout << "FAIL" << std::endl;
        throw file_access_error{"Could not read database file '" + filename + "'"};
    }

    return db;
}




/*****************************************************************************
 *
 * @brief writes database to file
 *
 *****************************************************************************/
template<class S, class K>
void
write_database(const sketch_database<S,K>& db, const std::string& filename)
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
template<class S, class K>
void print_info(const sketch_database<S,K>& db)
{
    using genome_id = typename sketch_database<S,K>::genome_id;
    int numRankedGenomes = 0;
    for(genome_id i = 0; i < db.genome_count(); ++i) {
        if(!db.taxon_of_genome(i).none()) ++numRankedGenomes;
    }

    std::cout << "sequences:        " << db.genome_count() << '\n'
              << "ranked sequences: " << numRankedGenomes << '\n'
              << "window length:    " << db.genome_window_size() << '\n'
              << "window stride:    " << db.genome_window_stride() << '\n'
              << "kmer size:        " << db.kmer_size() << '\n'
              << "sketch size:      " << db.genome_sketch_size() << '\n'
              << "taxa in tree:     " << db.taxon_count() << '\n';

    auto hbs = db.bucket_sizes();

    std::cout << "buckets:          " << db.bucket_count() << '\n'
              << "bucket size:      " << hbs.mean() << " +/- " << hbs.stddev()
              << std::endl;
}


} // namespace mc


#undef METACACHE_DB_VERSION

#endif

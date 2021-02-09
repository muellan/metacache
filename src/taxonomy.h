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

#ifndef MC_TAXONOMY_H_
#define MC_TAXONOMY_H_


#include "io_serialize.h"
#include "config.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <mutex>
#include <set>
#include <shared_mutex>
#include <unordered_map>
#include <vector>


namespace mc {


/*************************************************************************//**
 *
 * @brief id-based taxonomy
 *
 * @details the directed graph is stored implicitly
 *
 *****************************************************************************/
class taxonomy
{
public:
    //---------------------------------------------------------------
    using taxon_id   = std::int_least64_t;
    using taxon_name = std::string;

    static constexpr taxon_id none_id() noexcept { return 0; }


    //-----------------------------------------------------
    /**
     * @brief encodes taxnomic ranks
     */
    enum class rank : std::uint8_t {
        Sequence,
                Form,
                Variety,
                    subSpecies,
            Species,
                    subGenus,
            Genus,
                    subTribe,
                Tribe,
                    subFamily,
            Family,
                    subOrder,
            Order,
                    subClass,
            Class,
                    subPhylum,
            Phylum,
                    subKingdom,
            Kingdom,
            Domain,
        root,
        none
    };

    inline friend void write_binary(std::ostream& os, rank r) {
        write_binary(os, std::uint8_t(r));
    }
    inline friend void read_binary(std::istream& is, rank& r) {
        std::uint8_t n = 0;
        read_binary(is, n);
        r = rank(n);
    }

    //---------------------------------------------------------------
    static constexpr int num_ranks = static_cast<int>(rank::none);

    //---------------------------------------------------------------
    static rank
    next_main_rank(rank r) noexcept {
        switch(r) {
            case rank::Sequence:     return rank::Species;
            case rank::Form:         return rank::Species;
            case rank::Variety:      return rank::Species;
            case rank::subSpecies:   return rank::Species;
            case rank::Species:      return rank::Genus;
            case rank::subGenus:     return rank::Genus;
            case rank::Genus:        return rank::Family;
            case rank::subTribe:     return rank::Family;
            case rank::Tribe:        return rank::Family;
            case rank::subFamily:    return rank::Family;
            case rank::Family:       return rank::Order;
            case rank::subOrder:     return rank::Order;
            case rank::Order:        return rank::Class;
            case rank::subClass:     return rank::Class;
            case rank::Class:        return rank::Phylum;
            case rank::subPhylum:    return rank::Phylum;
            case rank::Phylum:       return rank::Kingdom;
            case rank::subKingdom:   return rank::Kingdom;
            case rank::Kingdom:      return rank::Domain;
            case rank::Domain:       return rank::root;
            default:
            case rank::root:
            case rank::none:         return rank::none;
        }
    }
    static rank
    prev_main_rank(rank r) noexcept {
        switch(r) {
            case rank::Sequence:     return rank::none;
            case rank::Form:         return rank::Sequence;
            case rank::Variety:      return rank::Sequence;
            case rank::subSpecies:   return rank::Sequence;
            case rank::Species:      return rank::Sequence;
            case rank::subGenus:     return rank::Species;
            case rank::Genus:        return rank::Species;
            case rank::subTribe:     return rank::Genus;
            case rank::Tribe:        return rank::Genus;
            case rank::subFamily:    return rank::Genus;
            case rank::Family:       return rank::Genus;
            case rank::subOrder:     return rank::Family;
            case rank::Order:        return rank::Family;
            case rank::subClass:     return rank::Order;
            case rank::Class:        return rank::Order;
            case rank::subPhylum:    return rank::Class;
            case rank::Phylum:       return rank::Class;
            case rank::subKingdom:   return rank::Phylum;
            case rank::Kingdom:      return rank::Phylum;
            case rank::Domain:       return rank::Kingdom;
            case rank::root:         return rank::Domain;
            default:
            case rank::none:         return rank::none;
        }
    }


    //---------------------------------------------------------------
    inline friend rank& operator ++ (rank& r) {
        return std::uint8_t(r) < num_ranks ? (r = rank(std::uint8_t(r) + 1)) : r;
    }
    inline friend rank& operator -- (rank& r) {
        return std::uint8_t(r) > 0 ? (r = rank(std::uint8_t(r) - 1)) : r;
    }

    inline friend rank operator + (rank r, std::uint8_t offset) {
        return std::uint8_t(r) < num_ranks ? rank(std::uint8_t(r) + offset) : rank::root;
    }
    inline friend rank operator - (rank r, std::uint8_t offset) {
        return std::uint8_t(r) >= offset ? rank(std::uint8_t(r) - offset) : rank::Sequence;
    }


    //---------------------------------------------------------------
    static rank
    rank_from_name(std::string name) noexcept {
        std::transform(name.begin(), name.end(), name.begin(), ::tolower);
        if(name == "sequence")         return rank::Sequence;
        if(name == "genome")           return rank::Sequence;
        if(name == "form")             return rank::Form;
        if(name == "forma")            return rank::Form;
        if(name == "variety")          return rank::Variety;
        if(name == "varietas")         return rank::Variety;
        if(name == "subspecies")       return rank::subSpecies;
        if(name == "species")          return rank::Species;
        if(name == "species group")    return rank::subGenus;
        if(name == "species subgroup") return rank::subGenus;
        if(name == "subgenus")         return rank::subGenus;
        if(name == "genus")            return rank::Genus;
        if(name == "subtribe")         return rank::subTribe;
        if(name == "tribe")            return rank::Tribe;
        if(name == "subfamily")        return rank::subFamily;
        if(name == "family")           return rank::Family;
        if(name == "superfamily")      return rank::subOrder;
        if(name == "parvorder")        return rank::subOrder;
        if(name == "infraorder")       return rank::subOrder;
        if(name == "suborder")         return rank::subOrder;
        if(name == "order")            return rank::Order;
        if(name == "superorder")       return rank::subClass;
        if(name == "infraclass")       return rank::subClass;
        if(name == "subclass")         return rank::subClass;
        if(name == "class")            return rank::Class;
        if(name == "superclass")       return rank::subPhylum;
        if(name == "subphylum")        return rank::subPhylum;
        if(name == "phylum")           return rank::Phylum;
        if(name == "division")         return rank::Phylum;
        if(name == "superphylum")      return rank::subKingdom;
        if(name == "subkingdom")       return rank::subKingdom;
        if(name == "kingdom")          return rank::Kingdom;
        if(name == "subdomain")        return rank::Kingdom;
        if(name == "superkingdom")     return rank::Domain;
        if(name == "domain")           return rank::Domain;
        if(name == "root")             return rank::root;
        return rank::none;
    }


    //---------------------------------------------------------------
    static const char*
    rank_name(rank r) {
        switch(r) {
            case rank::Sequence:     return "sequence";
            case rank::Form:         return "form";
            case rank::Variety:      return "variety";
            case rank::subSpecies:   return "subspecies";
            case rank::Species:      return "species";
            case rank::subGenus:     return "subgenus";
            case rank::Genus:        return "genus";
            case rank::subTribe:     return "subtribe";
            case rank::Tribe:        return "tribe";
            case rank::subFamily:    return "subfamily";
            case rank::Family:       return "family";
            case rank::subOrder:     return "suborder";
            case rank::Order:        return "order";
            case rank::subClass:     return "subclass";
            case rank::Class:        return "class";
            case rank::subPhylum:    return "subphylum";
            case rank::Phylum:       return "phylum";
            case rank::subKingdom:   return "subkingdom";
            case rank::Kingdom:      return "kingdom";
            case rank::Domain:       return "domain";
            case rank::root:         return "root";
            default:
            case rank::none:         return "none";
        }
    }


    /************************************************************
     *
     * @brief taxonomic node
     *
     ************************************************************/
    class taxon {
        friend class taxonomy;
    public:
        //-----------------------------------------------------
        using rank_type = taxonomy::rank;

        //-----------------------------------------------------
        struct file_source {
            using index_t   = std::uint_least64_t;
            using window_id = std::uint_least64_t;

            explicit
            file_source(std::string filename = "", index_t index = 0,
                        window_id numWindows = 0)
            :
                filename{std::move(filename)}, windows{numWindows}, index{index}
            {}

            std::string filename;
            window_id windows;
            index_t index;
        };

        //default: empty taxon
        explicit
        taxon(taxon_id taxonId = none_id(),
              taxon_id parentId = none_id(),
              std::string taxonName = "--",
              rank_type rk = rank_type::none,
              file_source source = file_source{})
        :
            id_{taxonId}, parent_{parentId},
            name_{std::move(taxonName)},
            source_{std::move(source)},
            rank_{rk}
        {}

        //makes taxa sortable
        inline friend bool
        operator < (const taxon& a, const taxon& b) noexcept {
            return a.id_ < b.id_;
        }

        taxon_id id() const noexcept { return id_; }

        const taxon_name& name() const noexcept { return name_; }

        rank_type rank() const noexcept { return rank_; }

        const file_source& source() const noexcept { return source_; }

        const char* rank_name() const noexcept {
            return taxonomy::rank_name(rank_);
        }

        taxon_id parent_id() const noexcept { return parent_; }

        bool has_parent() const noexcept {
            return parent_ != taxonomy::none_id();
        }


        //-----------------------------------------------------
        friend
        void read_binary(std::istream& is, taxon& t) {
            read_binary(is, t.id_);
            read_binary(is, t.parent_);
            read_binary(is, t.rank_);
            read_binary(is, t.name_);
            read_binary(is, t.source_.filename);
            read_binary(is, t.source_.index);
            read_binary(is, t.source_.windows);
        }

        //-----------------------------------------------------
        friend
        void write_binary(std::ostream& os, const taxon& t) {
            write_binary(os, t.id_);
            write_binary(os, t.parent_);
            write_binary(os, t.rank_);
            write_binary(os, t.name_);
            write_binary(os, t.source_.filename);
            write_binary(os, t.source_.index);
            write_binary(os, t.source_.windows);
        }

        //-----------------------------------------------------
    private:
        taxon_id id_;
        taxon_id parent_;
        taxon_name name_;
        file_source source_;
        rank_type rank_;
    };


private:
    using taxon_store = std::set<taxon>;

public:
    //-----------------------------------------------------
    using full_lineage   = std::vector<const taxon*>;
    using ranked_lineage = std::array<const taxon*,num_ranks>;
    using file_source    = taxon::file_source;
    using size_type      = taxon_store::size_type;
    using const_iterator = taxon_store::const_iterator;
    //-----------------------------------------------------
    class const_range {
    public:
        explicit
        const_range(const_iterator beg, const_iterator end):
            beg_{beg}, end_{end}
        {}
        const_iterator begin() const { return beg_; }
        const_iterator end() const   { return end_; }
    private:
        const_iterator beg_;
        const_iterator end_;
    };


    //-------------------------------------------------------------------
    void clear() {
        taxa_.clear();
    }
    //-----------------------------------------------------
    /**
     * @brief erase all taxa with rank > r (except taxon with id = 0)
     */
    void erase_above(rank r) {
        auto i = taxa_.begin();
        while(i != taxa_.end()) {
            if(i->rank() > r) {
                i = taxa_.erase(i);
            } else {
                ++i;
            }
        }
    }


    //-------------------------------------------------------------------
    const_iterator
    emplace(taxon_id taxonId,
            taxon_id parentId,
            std::string taxonName,
            std::string rankName,
            file_source source = file_source{})
    {
        return emplace(taxonId, parentId, std::move(taxonName),
                       rank_from_name(rankName),
                       std::move(source));
    }
    //-----------------------------------------------------
    const_iterator
    emplace(taxon_id taxonId,
            taxon_id parentId = none_id(),
            std::string taxonName = "",
            rank rank = rank::none,
            file_source source = file_source{})
    {
        if(taxonId == none_id()) return taxa_.end();

        return taxa_.emplace(taxonId, parentId,
                             std::move(taxonName), rank,
                             std::move(source)).first;
    }


    //---------------------------------------------------------------
    const_iterator
    insert(const taxon& t) {
        return taxa_.insert(t).first;
    }
    //-----------------------------------------------------
    const_iterator
    insert(taxon&& t) {
        return taxa_.insert(std::move(t)).first;
    }


    //---------------------------------------------------------------
    const_iterator
    insert_or_replace(const taxon& t)
    {
        auto i = taxa_.find(taxon{t.id()});
        if(i != taxa_.end()) {
            taxa_.erase(i);
            return taxa_.insert(t).first;
        } else {
            return taxa_.insert(t).first;
        }
    }
    //-----------------------------------------------------
    const_iterator
    insert_or_replace(taxon&& t)
    {
        auto i = taxa_.find(taxon{t.id()});
        if(i != taxa_.end()) {
            taxa_.erase(i);
            return taxa_.insert(std::move(t)).first;
        } else {
            return taxa_.insert(t).first;
        }
    }


    //---------------------------------------------------------------
    bool reset_parent(taxon_id id, taxon_id parent)
    {
        if(id == none_id()) return false;
        auto i = taxa_.find(taxon{id});
        if(i == taxa_.end()) return false;
        const_cast<taxon*>(&*i)->parent_ = parent;
        return true;
    }


    //---------------------------------------------------------------
    bool reset_rank(taxon_id id, rank rank)
    {
        if(id == none_id()) return false;
        auto i = taxa_.find(taxon{id});
        if(i == taxa_.end()) return false;
        const_cast<taxon*>(&*i)->rank_ = rank;
        return true;
    }


    //---------------------------------------------------------------
    const_iterator
    find(taxon_id id) const {
        return taxa_.find(taxon{id});
    }
    /**
     * @return taxon to 'id' or nullptr;
     */
    const taxon*
    operator [] (taxon_id id) const {
        auto i = find(id);
        if(i == taxa_.end()) return nullptr;
        return &(*i);
    }


    //---------------------------------------------------------------
    bool
    contains(taxon_id id) const {
        return find(id) != taxa_.end();
    }


    //---------------------------------------------------------------
    const taxon*
    lca(const taxon& a, const taxon& b) const
    {
        return lca(make_lineage(a), make_lineage(b));
    }
    //---------------------------------------------------------------
    const taxon*
    lca(const full_lineage& lina, const full_lineage& linb) const
    {
        for(auto ta : lina) {
            for(auto tb : linb) {
                if(ta && tb && (ta == tb)) return ta;
            }
        }
        return nullptr;
    }


    //---------------------------------------------------------------
    const taxon*
    parent(const taxon& tax) const {
        return operator[](tax.parent_);
    }


    //---------------------------------------------------------------
    /** @brief return next taxon that has a rank assigned */
    const taxon*
    next_ranked_ancestor(const taxon& tax) const {
        if(tax.rank() != rank::none) return &tax;
        return next_ranked_ancestor(tax.id());
    }
    //-----------------------------------------------------
    const taxon*
    next_ranked_ancestor(taxon_id id) const
    {
        while(id != none_id()) {
            auto it = taxa_.find(taxon{id});
            if(it == taxa_.end()) break;
            if(it->rank() != rank::none) {
                return &(*it);
            }
            if(it->parent_ == id) break; //break cycles
            id = it->parent_;
        }
        return nullptr;
    }


    //---------------------------------------------------------------
    ranked_lineage
    make_ranks(const taxon& tax) const {
        return make_ranks(tax.id());
    }
    //-----------------------------------------------------
    ranked_lineage
    make_ranks(taxon_id id) const
    {
        auto lin = ranked_lineage{};
        for(auto& x : lin) x = nullptr;

        while(id != none_id()) {
            auto it = taxa_.find(taxon{id});
            if(it == taxa_.end()) break;
            if(it->rank() != rank::none) {
                lin[static_cast<int>(it->rank())] = &(*it);
            }
            if(it->parent_ == id) break; //break cycles
            id = it->parent_;
        }
        return lin;
    }


    //---------------------------------------------------------------
    full_lineage
    make_lineage(const taxon& tax) const {
        return make_lineage(tax.id());
    }
    //-----------------------------------------------------
    full_lineage
    make_lineage(taxon_id id) const
    {
        auto lin = full_lineage{};

        while(id != none_id()) {
            auto it = taxa_.find(taxon{id});
            if(it != taxa_.end()) {
                lin.push_back(&(*it));
                if(it->parent_ != id) { //break cycles
                    id = it->parent_;
                } else {
                    id = none_id();
                }
            } else {
                id = none_id();
            }
        }
        return lin;
    }


    //---------------------------------------------------------------
    size_type size() const noexcept {
        return taxa_.size();
    }
    //-----------------------------------------------------
    bool empty() const noexcept {
        return taxa_.empty();
    }


    //---------------------------------------------------------------
    const_iterator begin() const { return taxa_.begin(); }
    const_iterator end()   const { return taxa_.end(); }
    //-----------------------------------------------------
    const_range
    full_range() const {
        return const_range{taxa_.begin(), taxa_.end()};
    }
    //-----------------------------------------------------
    const_range
    subrange(taxon_id first, taxon_id last) const {
        return const_range{taxa_.lower_bound(taxon{first}),
                           taxa_.upper_bound(taxon{last})};
    }
    //-----------------------------------------------------
    const_range
    subrange_from(taxon_id first) const {
        return const_range{taxa_.lower_bound(taxon{first}), taxa_.end()};
    }
    //-----------------------------------------------------
    const_range
    subrange_until(taxon_id last) const {
        return const_range{taxa_.begin(), taxa_.upper_bound(taxon{last})};
    }


    //---------------------------------------------------------------
    friend void
    read_binary(std::istream& is, taxonomy& tax)
    {
        tax.taxa_.clear();
        std::uint64_t nt = 0;
        read_binary(is, nt);
        for(std::uint64_t i = 0; i < nt; ++i) {
            taxon t;
            read_binary(is, t);
            tax.taxa_.insert(std::move(t));
        }
    }

    //-----------------------------------------------------
    friend void
    write_binary(std::ostream& os, const taxonomy& tax)
    {
        write_binary(os, std::uint64_t(tax.taxa_.size()));
        for(const auto& t : tax.taxa_) {
            write_binary(os, t);
        }
    }


private:
    //---------------------------------------------------------------
    taxon_store taxa_;
};




/*************************************************************************//**
 *
 * @brief ranked lineage cache
 *
 *****************************************************************************/
class ranked_lineages_cache
{
public:
    //---------------------------------------------------------------
    using taxon          = taxonomy::taxon;
    using ranked_lineage = taxonomy::ranked_lineage;
    using taxon_rank     = taxonomy::rank;

public:
    //---------------------------------------------------------------
    /**
     * @param highestRank  initialize cache for taxa below and at this rank
     *                     and for their ranked ancestors
     */
    explicit
    ranked_lineages_cache(const taxonomy& taxa)
    :
        taxa_(taxa), highestRank_{taxon_rank::Sequence},
        empty_{}, lins_{},
        outdated_(true)
    {
        for(auto& x : empty_) x = nullptr;
    }


    //---------------------------------------------------------------
    ranked_lineages_cache(const ranked_lineages_cache&) = delete;

    ranked_lineages_cache(ranked_lineages_cache&& src):
        taxa_(src.taxa_), highestRank_{src.highestRank_},
        lins_{std::move(src.lins_)},
        outdated_(src.outdated_)
    {}

    ranked_lineages_cache& operator = (const ranked_lineages_cache&) = delete;
    ranked_lineages_cache& operator = (ranked_lineages_cache&&) = delete;

    //---------------------------------------------------------------
    void mark_outdated() {
        outdated_ = true;
    }

    //---------------------------------------------------------------
    void
    init_from_targets(const std::vector<ranked_lineage>& lins) {
        clear();
        highestRank_ = taxon_rank::Sequence;
        for(const auto& lin : lins)
            insert_lineage(lin);
    }

    //---------------------------------------------------------------
    /**
     * @brief  insert lineage for all its ranked ancestors
     */
    void
    insert_lineage(ranked_lineage lin) {
        for(int i = static_cast<int>(taxon_rank::Sequence);
                i < static_cast<int>(taxon_rank::root);
                ++i)
        {
            const taxon* tax = lin[i];
            if(tax) {
                // insert if not present
                if(lins_.find(tax) == lins_.end())
                    lins_.emplace(tax, lin);
                // remove this rank for next ancestor
                lin[i] = nullptr;
            }
        }
    }

    //---------------------------------------------------------------
    /**
     * @brief  insert lineage for taxon and all its ranked ancestors
     */
    void
    insert(const taxon& tax) {
        auto i = lins_.find(&tax);
        if(i == lins_.end()) {
            auto lin = lins_.emplace(&tax, taxa_.make_ranks(tax)).first->second;
            insert_lineage(lin);
        }
    }

    //-----------------------------------------------------
    /**
     * @brief    increase highest rank of cache, then update if outdated
     *
     * @details  if new rank is not higher than previous rank,
     *           only update if already outdated
     *
     * @param highestRank  initialize cache for taxa below and at this rank
     *                     and for their ranked ancestors
     */
    void update(taxon_rank highestRank = taxon_rank::Sequence) {
        if(highestRank_ < highestRank) {
            highestRank_ = highestRank;
            outdated_ = true;
        }

        if(!outdated_) return;

        for(const auto& t : taxa_) {
            if(t.rank() <= highestRank_) {
                insert(t);
            }
        }
        outdated_ = false;
    }

    //-----------------------------------------------------
    void clear() {
        lins_.clear();
        outdated_ = true;
    }

    //-----------------------------------------------------
    /**
     * @brief  rebuild cache using previous highest rank
     */
    void reset() {
        clear();
        update();
    }

    //-----------------------------------------------------
    /**
     * @brief  rebuild cache using new highest rank
     */
    void reset(taxon_rank highestRank) {
        clear();
        highestRank_ = highestRank;
        update();
    }

    //---------------------------------------------------------------
    /// @brief only works if tax is cached - make sure to call update first
    const ranked_lineage&
    operator [](const taxon* tax) const {
        return tax ? operator [](*tax) : empty_;
    }
    //-----------------------------------------------------
    /// @brief only works if tax is cached - make sure to call update first
    const ranked_lineage&
    operator [](const taxon& tax) const {
        assert(outdated_ == false);
        auto i = lins_.find(&tax);
        if(i != lins_.end()) return i->second;
        return empty_;
    }


private:
    //---------------------------------------------------------------
    const taxonomy& taxa_;
    taxon_rank highestRank_;
    ranked_lineage empty_;
    mutable std::unordered_map<const taxon*,ranked_lineage> lins_;
    bool outdated_;
};




/*************************************************************************//**
*
* @brief ranked lineage cache
*
*****************************************************************************/
class ranked_lineages_of_targets
{
public:
    //---------------------------------------------------------------
    using ranked_lineage = taxonomy::ranked_lineage;
    using taxon          = taxonomy::taxon;
    using taxon_id       = taxonomy::taxon_id;
    using taxon_rank     = taxonomy::rank;

    //use negative numbers for sequence level taxon ids
    static constexpr taxon_id
    taxon_id_of_target(target_id id) noexcept { return -taxon_id(id)-1; }

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
    size_t target_count() const noexcept {
        return lins_.size();
    }


    //---------------------------------------------------------------
    void mark_outdated() {
        outdated_ = true;
    }

    //---------------------------------------------------------------
    /**
     * @brief  update cache if target count increased or already outdated
     */
    void update(target_id numTargets = 0) {
        if(numTargets > lins_.size()) {
            outdated_ = true;
        }

        if(!outdated_) return;

        // add new target ranks
        const target_id oldSize = lins_.size();
        lins_.resize(numTargets);


        for(target_id tgt = oldSize; tgt < numTargets; ++tgt) {
            lins_[tgt] = taxa_.make_ranks(taxon_id_of_target(tgt));
        }
        outdated_ = false;
    }

    //---------------------------------------------------------------
    /**
     * @brief  rebuild cache for given target count
     */
    void reset(target_id numTargets = 0) {
        if(numTargets == 0) numTargets = lins_.size();
        clear();
        update(numTargets);
    }

    //-----------------------------------------------------
    void clear() {
        lins_.clear();
        outdated_ = true;
    }

    //---------------------------------------------------------------
    const ranked_lineage&
    operator [](target_id tgt) const {
        assert(outdated_ == false);
        return lins_[tgt];
    }

    //---------------------------------------------------------------
    const std::vector<ranked_lineage>& lineages() const {
        return lins_;
    }


private:
    //---------------------------------------------------------------
    const taxonomy& taxa_;
    std::vector<ranked_lineage> lins_;
    bool outdated_;
};




/*************************************************************************//**
*
* @brief wrapper for taxonomy and lineage caches
*
*****************************************************************************/
class taxonomy_cache
{
public:
    using taxon          = taxonomy::taxon;
    using taxon_id       = taxonomy::taxon_id;
    using taxon_name     = taxonomy::taxon_name;
    using taxon_rank     = taxonomy::rank;
    using taxon_range    = taxonomy::const_range;
    using ranked_lineage = taxonomy::ranked_lineage;
    using full_lineage   = taxonomy::full_lineage;
    using file_source    = taxonomy::file_source;


    //---------------------------------------------------------------
    taxonomy_cache() :
        taxa_{},
        taxonLineages_{taxa_},
        targetLineages_{taxa_},
        name2tax_{}
    {}

    taxonomy_cache(const taxonomy_cache&) = delete;
    taxonomy_cache(taxonomy_cache&& other) :
        taxa_{std::move(other.taxa_)},
        taxonLineages_{std::move(other.taxonLineages_)},
        targetLineages_{std::move(other.targetLineages_)},
        name2tax_{std::move(other.name2tax_)}
    {}


    //---------------------------------------------------------------
    bool empty() const noexcept {
        return taxa_.empty();
    }

    //---------------------------------------------------------------
    std::uint64_t target_count() const noexcept {
        return targetLineages_.target_count();
    }

    //---------------------------------------------------------------
    const std::vector<ranked_lineage>& target_lineages() const {
        return targetLineages_.lineages();
    }


    //---------------------------------------------------------------
    void clear() {
        taxonLineages_.clear();
        targetLineages_.clear();
        name2tax_.clear();
    }


    //---------------------------------------------------------------
    const taxon* cached_taxon_of_target(target_id id) const {
        return targetLineages_[id][0];
    }
    //-----------------------------------------------------
    const taxon*
    taxon_with_id(taxon_id id) const noexcept {
        return taxa_[id];
    }
    /******************************************************
     * @brief will only find sequence-level taxon names == sequence id
     */
    const taxon*
    taxon_with_name(const taxon_name& name) const noexcept {
        if(name.empty()) return nullptr;
        auto i = name2tax_.find(name);
        if(i == name2tax_.end()) return nullptr;
        return i->second;
    }
    /******************************************************
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


    /****************************************************************
     * @brief concurrency-safe taxon name lookup
     */
    bool contains_name(const taxon_name& sid) {
        std::shared_lock<std::shared_timed_mutex> lock(name2taxMtx);
        return name2tax_.find(sid) != name2tax_.end();
    }
    /****************************************************************
     * @brief concurrency-safe taxon name insert
     */
    void insert_name(taxon_name sid, const taxon* tax) {
        std::unique_lock<std::shared_timed_mutex> lock(name2taxMtx);
        name2tax_.insert({std::move(sid), tax});
    }


    /****************************************************************
     * @brief  concurrency-safe taxon emplacement
     */
    const taxon*
    emplace_taxon(taxon_id taxid, taxon_id parentTaxid, taxon_name sid,
                  taxon_rank rank, file_source source)
    {
        std::lock_guard<std::mutex> lock(taxaMtx);

        auto nit = taxa_.emplace(
            taxid, parentTaxid, std::move(sid),
            rank, std::move(source));

        // should never happen
        if(nit == taxa_.end()) {
            throw std::runtime_error{"target taxon could not be created"};
        }

        return &(*nit);
    }


    //---------------------------------------------------------------
    std::uint64_t
    non_target_taxon_count() const noexcept {
        return taxa_.size() - target_count();
    }
    //---------------------------------------------------------------
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
    full_lineage
    make_lineage(const taxon* tax) const noexcept {
        return tax ? make_lineage(*tax) : full_lineage();
    }
    //-----------------------------------------------------
    full_lineage
    make_lineage(const taxon& tax) const noexcept {
        return taxa_.make_lineage(tax);
    }

    //---------------------------------------------------------------
    const ranked_lineage
    make_ranks(const taxon* tax) const noexcept {
        return tax ? make_ranks(*tax) : ranked_lineage{};
    }
    // -----------------------------------------------------
    const ranked_lineage
    make_ranks(const taxon& tax) const noexcept {
        return taxa_.make_ranks(tax);
    }
    //---------------------------------------------------------------
    const ranked_lineage&
    cached_ranks(const taxon* tax) const noexcept {
        return taxonLineages_[tax];
    }
    // -----------------------------------------------------
    const ranked_lineage&
    cached_ranks(const taxon& tax) const noexcept {
        return taxonLineages_[tax];
    }
    //-----------------------------------------------------
    const ranked_lineage&
    cached_ranks(target_id tgt) const noexcept {
        return targetLineages_[tgt];
    }

    //---------------------------------------------------------------
    const taxon*
    parent(const taxon* tax) const noexcept {
        return tax ? parent(*tax) : nullptr;
    }
    //-----------------------------------------------------
    const taxon*
    parent(const taxon& tax) const noexcept {
        return taxa_.parent(tax);
    }

    //---------------------------------------------------------------
    const taxon*
    cached_taxon_ancestor(const taxon* tax, taxon_rank r) const noexcept {
        return tax ? cached_taxon_ancestor(*tax,r) : nullptr;
    }
    //-----------------------------------------------------
    const taxon*
    cached_taxon_ancestor(const taxon& tax, taxon_rank r) const noexcept {
        return taxonLineages_[tax][int(r)];
    }
    //-----------------------------------------------------
    const taxon*
    cached_ancestor(target_id tgt, taxon_rank r) const noexcept {
        return targetLineages_[tgt][int(r)];
    }

    //---------------------------------------------------------------
    const taxon*
    cached_next_ranked_ancestor(const taxon* tax) const noexcept {
        return tax ? cached_next_ranked_ancestor(*tax) : nullptr;
    }
    //-----------------------------------------------------
    const taxon*
    cached_next_ranked_ancestor(const taxon& tax) const noexcept {
        if(tax.rank() != taxon_rank::none) return &tax;
        for(const taxon* a : taxonLineages_[tax]) {
            if(a) return a;
        }
        return nullptr;
    }

    //---------------------------------------------------------------
    const taxon*
    lowest_ranked_ancestor(target_id tgt, taxon_rank lowest) const {
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
    //-----------------------------------------------------
    const taxon*
    lca(const taxon& ta, const taxon& tb) const {
        return taxa_.lca(ta,tb);
    }

     //-----------------------------------------------------
    const taxon*
    ranked_lca(const ranked_lineage& lina, const ranked_lineage& linb,
               taxon_rank lowest = taxon_rank::Sequence) const
    {
        for(int i  = static_cast<int>(lowest);
                i <= static_cast<int>(taxon_rank::root);
                ++i)
        {
            if(lina[i] && (lina[i] == linb[i])) return lina[i];
        }
        return nullptr;
    }
    //---------------------------------------------------------------
    const taxon*
    make_ranked_lca(const taxon* ta, const taxon* tb) const {
        return (ta && tb) ? make_ranked_lca(*ta,*tb) : nullptr;
    }
    //---------------------------------------------------------------
    const taxon*
    make_ranked_lca(const taxon& ta, const taxon& tb) const
    {
        return ranked_lca(taxa_.make_ranks(ta), taxa_.make_ranks(tb));
    }
    //---------------------------------------------------------------
    const taxon*
    cached_ranked_lca(const taxon* ta, const taxon* tb) const {
        return (ta && tb) ? cached_ranked_lca(*ta,*tb) : nullptr;
    }
    //-----------------------------------------------------
    const taxon*
    cached_ranked_lca(const taxon& ta, const taxon& tb) const {
        return ranked_lca(taxonLineages_[ta], taxonLineages_[tb]);
    }
    //-----------------------------------------------------
    const taxon*
    cached_ranked_lca(target_id ta, target_id tb, taxon_rank lowest) const {
        return ranked_lca(targetLineages_[ta], targetLineages_[tb], lowest);
    }


    /****************************************************************
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
                for(const taxon* t : taxa_.make_lineage(tax)) {
                    if(t == &covered) ++coverage;
                }
            }
        }
        return coverage;
    }

    /****************************************************************
     * @return true, if taxon is covered by any target in the DB
     */
    bool covers(const taxon* covered) const {
        return covered ? covers(*covered) : false;
    }
    //-----------------------------------------------------
    bool covers(const taxon& covered) const {
        for(const auto& tax : taxa_) {
            if(tax.rank() == taxon_rank::Sequence) {
                for(const taxon* t : taxa_.make_lineage(tax)) {
                    if(t == &covered) return true;
                }
            }
        }
        return false;
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
    void reset_taxa_above_sequence_level(taxonomy&& tax) {
        taxa_.erase_above(taxon_rank::Sequence);
        for(auto& t : tax) {
            taxa_.insert_or_replace(std::move(t));
        }
        //re-initialize ranks cache
        targetLineages_.reset();
        taxonLineages_.init_from_targets(targetLineages_.lineages());
    }


    //---------------------------------------------------------------
    void initialize_caches(std::uint64_t targetCount) {
        targetLineages_.reset(targetCount);
        taxonLineages_.init_from_targets(targetLineages_.lineages());

        //sequence id lookup
        for(const auto& t : target_taxa()) {
            name2tax_.insert({t.name(), &t});
        }
    }


    //---------------------------------------------------------------
    void mark_cached_lineages_outdated() {
        taxonLineages_.mark_outdated();
        targetLineages_.mark_outdated();
    }

    /****************************************************************
     * @details if target count is 0, the prvious target count will be used
     */
    void update_cached_lineages(taxon_rank rank, std::uint64_t targetCount = 0) const {
        taxonLineages_.update(rank);
        targetLineages_.update(targetCount);
    }


    friend void
    read_binary(std::istream& is, taxonomy_cache& taxonomyCache) {
        taxonomyCache.clear();
        read_binary(is, taxonomyCache.taxa_);
    }

    friend void
    write_binary(std::ostream& os, const taxonomy_cache& taxonomyCache) {
        write_binary(os, taxonomyCache.taxa_);
    }

private:
    //---------------------------------------------------------------
    taxonomy taxa_;
    mutable ranked_lineages_cache taxonLineages_;
    mutable ranked_lineages_of_targets targetLineages_;
    std::map<taxon_name,const taxon*> name2tax_;

    std::shared_timed_mutex name2taxMtx;
    std::mutex taxaMtx;
};


/*************************************************************************//**
 *
 * @brief pull some types from taxonomy into global namespace
 *
 *****************************************************************************/
using taxon          = taxonomy::taxon;
using taxon_rank     = taxonomy::rank;
using taxon_id       = taxonomy::taxon_id;
using ranked_lineage = taxonomy::ranked_lineage;


} // namespace mc


#endif

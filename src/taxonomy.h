/*****************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2017 André Müller (muellan@uni-mainz.de)
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

#include <array>
#include <vector>
#include <set>
#include <algorithm>
#include <cstdint>
#include <iostream>

#include "io_serialize.h"


namespace mc {


/*****************************************************************************
 *
 * @brief id-based taxonomy
 *
 * @details the directed graph is stored implicitly
 *
 * @invariant has one "empty taxon" at all times
 *
 *****************************************************************************/
class taxonomy
{
public:
    //---------------------------------------------------------------
    using taxon_id  = std::uint64_t;


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
            case rank::none: return rank::none;
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
            case rank::none: return rank::none;
        }
    }


    //---------------------------------------------------------------
    inline friend rank& operator ++ (rank& r) {
        return int(r) < num_ranks ? (r = rank(int(r) + 1)) : r;
    }
    inline friend rank& operator -- (rank& r) {
        return int(r) > 1 ? (r = rank(int(r) - 1)) : r;
    }


    //---------------------------------------------------------------
    static rank
    rank_from_name(const std::string& name) noexcept {
        if(name == "sequence")      return rank::Sequence;
        if(name == "genome")        return rank::Sequence;
        if(name == "form")          return rank::Form;
        if(name == "forma")         return rank::Form;
        if(name == "variety")       return rank::Variety;
        if(name == "varietas")      return rank::Variety;
        if(name == "subspecies")    return rank::subSpecies;
        if(name == "species")       return rank::Species;
        if(name == "subgenus")      return rank::subGenus;
        if(name == "genus")         return rank::Genus;
        if(name == "subtribe")      return rank::subTribe;
        if(name == "tribe")         return rank::Tribe;
        if(name == "subfamily")     return rank::subFamily;
        if(name == "family")        return rank::Family;
        if(name == "superfamily")   return rank::subOrder;
        if(name == "parvorder")     return rank::subOrder;
        if(name == "infraorder")    return rank::subOrder;
        if(name == "suborder")      return rank::subOrder;
        if(name == "order")         return rank::Order;
        if(name == "superorder")    return rank::subClass;
        if(name == "infraclass")    return rank::subClass;
        if(name == "subclass")      return rank::subClass;
        if(name == "class")         return rank::Class;
        if(name == "superclass")    return rank::subPhylum;
        if(name == "subphylum")     return rank::subPhylum;
        if(name == "phylum")        return rank::Phylum;
        if(name == "division")      return rank::Phylum;
        if(name == "superphylum")   return rank::subKingdom;
        if(name == "subkingdom")    return rank::subKingdom;
        if(name == "kingdom")       return rank::Kingdom;
        if(name == "superkingdom")  return rank::Domain;
        if(name == "domain")        return rank::Domain;
        if(name == "root")          return rank::root;
        return rank::none;
    }


    //---------------------------------------------------------------
    static std::string
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


    //-----------------------------------------------------
    /**
     * @brief taxonomic node information
     */
    class taxon {
        friend class taxonomy;
    public:
        using rank_type = taxonomy::rank;

        explicit
        taxon(taxon_id taxonId = 0,
              taxon_id parentId = 0,
              std::string taxonName = "--",
              rank_type rk = rank_type::none)
        :
            id_(taxonId), parent_(parentId),
            rank_(rk),
            name_(std::move(taxonName))
        {}

        //makes taxa sortable
        inline friend bool
        operator < (const taxon& a, const taxon& b) noexcept {
            return a.id_ < b.id_;
        }

        bool none() const noexcept { return id_ < 2; }

        taxon_id id() const noexcept { return id_; }

        std::string name() const noexcept { return name_; }

        rank_type rank() const noexcept { return rank_; }

        std::string rank_name() const noexcept {
            return taxonomy::rank_name(rank_);
        }

        //-----------------------------------------------------
        friend
        void read_binary(std::istream& is, taxon& t) {
            read_binary(is, t.id_);
            read_binary(is, t.parent_);
            read_binary(is, t.rank_);
            read_binary(is, t.name_);
        }

        //-----------------------------------------------------
        friend
        void write_binary(std::ostream& os, const taxon& t) {
            write_binary(os, t.id_);
            write_binary(os, t.parent_);
            write_binary(os, t.rank_);
            write_binary(os, t.name_);
        }

        //-----------------------------------------------------
    private:
        taxon_id id_;
        taxon_id parent_;
        rank_type rank_;
        std::string name_;
    };


private:
    using taxon_store = std::set<taxon>;


public:
    //-----------------------------------------------------
    using size_type      = taxon_store::size_type;
    using const_iterator = taxon_store::const_iterator;
    using full_lineage   = std::vector<taxon_id>;
    using ranked_lineage = std::array<taxon_id,num_ranks>;


    //-------------------------------------------------------------------
    taxonomy():
        taxa_{}
    {
        //insert one "empty taxon"
        taxa_.emplace();
    }


    //-------------------------------------------------------------------
    void
    clear() {
        taxa_.clear();
        //insert one "empty taxon"
        taxa_.emplace();
    }


    //-------------------------------------------------------------------
    bool
    emplace(taxon_id taxonId,
            taxon_id parentId,
            const std::string& taxonName,
            const std::string& rankName)
    {
        return emplace(taxonId, parentId, taxonName, rank_from_name(rankName));
    }
    //-----------------------------------------------------
    bool
    emplace(taxon_id taxonId,
            taxon_id parentId = 0,
            const std::string& taxonName = "",
            rank rank = rank::none)
    {
        if(taxonId < 1) return false;
        taxa_.emplace(taxonId, parentId, taxonName, rank);
        return true;
    }
    //-----------------------------------------------------
    bool
    insert(const taxon& t) {
        if(t.id_ < 1) return false;
        taxa_.insert(t);
        return true;
    }
    //-----------------------------------------------------
    bool
    insert(taxon&& t) {
        if(t.id_ < 1) return false;
        taxa_.insert(std::move(t));
        return true;
    }


    //---------------------------------------------------------------
    /**
     * @return taxon to 'id' or 'empty taxon' with id = 0 otherwise
     */
    const taxon&
    operator [] (taxon_id id) const {
        auto it = taxa_.find(taxon{id});
        return (it != taxa_.end()) ? *it : *taxa_.begin();
    }

    //-----------------------------------------------------
    bool
    contains(taxon_id id) const {
        return taxa_.find(taxon{id}) != taxa_.end();
    }


    //---------------------------------------------------------------
    const taxon&
    lca(const taxon& a, const taxon& b) const
    {
        return lca(a.id_, b.id_);
    }

    //---------------------------------------------------------------
    const taxon&
    ranked_lca(const taxon& a, const taxon& b) const
    {
        return ranked_lca(a.id_, b.id_);
    }


    //---------------------------------------------------------------
    const taxon&
    lca(const full_lineage& lina, const full_lineage& linb) const
    {
        return operator[](lca_id(lina,linb));
    }

    //-----------------------------------------------------
    const taxon&
    lca(const ranked_lineage& lina, const ranked_lineage& linb) const
    {
        return operator[](ranked_lca_id(lina,linb));
    }

    //-----------------------------------------------------
    const taxon&
    parent(const taxon& tax) const noexcept {
        return operator[](tax.parent_);
    }


    //---------------------------------------------------------------
    ranked_lineage
    ranks(const taxon& tax) const {
        return ranks(tax.id_);
    }

    //-----------------------------------------------------
    ranked_lineage
    ranks(taxon_id id) const
    {
        auto lin = ranked_lineage{};
        for(auto& x : lin) x = 0;

        while(id) {
            auto it = taxa_.find(taxon{id});
            if(it != taxa_.end()) {
                if(it->rank_ != rank::none) {
                    lin[static_cast<int>(it->rank_)] = it->id_;
                }
                if(it->parent_ != id) {
                    id = it->parent_;
                } else {
                    id = 0;
                }
            } else {
                id = 0;
            }
        }

        return lin;
    }

    //---------------------------------------------------------------
    full_lineage
    lineage(const taxon& tax) const {
        return lineage(tax.id_);
    }

    //-----------------------------------------------------
    full_lineage
    lineage(taxon_id id) const
    {
        auto lin = full_lineage{};

        while(id) {
            auto it = taxa_.find(taxon{id});
            if(it != taxa_.end()) {
                lin.push_back(id);
                if(it->parent_ != id) {
                    id = it->parent_;
                } else {
                    id = 0;
                }
            } else {
                id = 0;
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
    const taxon&
    lca(taxon_id a, taxon_id b) const
    {
        return operator[](lca_id(lineage(a), lineage(b) ));
    }
    //-----------------------------------------------------
    const taxon&
    ranked_lca(taxon_id a, taxon_id b) const
    {
        return operator[](ranked_lca_id(ranks(a), ranks(b) ));
    }


    //---------------------------------------------------------------
    static taxon_id
    ranked_lca_id(const ranked_lineage& lina,
                  const ranked_lineage& linb)
    {
        for(int i = 0; i < int(rank::root); ++i) {
            if((lina[i] > 0) && (lina[i] == linb[i])) return lina[i];
        }

        return 0;
    }

    //-----------------------------------------------------
    static taxon_id
    lca_id(const full_lineage& lina,
           const full_lineage& linb)
    {
        for(auto ta : lina) {
            for(auto tb : linb) {
                if(ta == tb) return ta;
            }
        }

        return 0;
    }

    //---------------------------------------------------------------
    rank
    lowest_rank(taxon_id id) const
    {
       while(id) {
            auto it = taxa_.find(taxon{id});
            if(it != taxa_.end()) {
                if(it->rank_ != rank::none) {
                    return it->rank_;
                }
                if(it->parent_ != id) {
                    id = it->parent_;
                } else {
                    id = 0;
                }
            } else {
                id = 0;
            }
        }
        return rank::none;
    }
    //-----------------------------------------------------
    rank
    lowest_rank(const taxon& tax) const {
        return lowest_rank(tax.id_);
    }

    //---------------------------------------------------------------
    taxon_store taxa_;
};


} // namespace mc


#endif

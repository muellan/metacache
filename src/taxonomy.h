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

#ifndef MC_TAXONOMY_H_
#define MC_TAXONOMY_H_

#include <array>
#include <vector>
#include <set>
#include <algorithm>
#include <cstdint>
#include <iostream>

#include "confusion.h"
#include "io_serialize.h"


namespace mc {


/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
class taxonomy
{
public:
    //---------------------------------------------------------------
    using taxon_id = std::uint64_t;

    //-----------------------------------------------------
    enum class taxon_rank : std::uint8_t {
        Sequence,
                Form,
                Variety,
                    subSpecies,
            Species,
                    subGenus,
            Genus,
//                    subTribe,
                Tribe,
//                    subFamily,
            Family,
//                    superFamily,
//                    subOrder,
            Order,
//                    superOrder,
//                    subClass,
            Class,
//                    superClass,
//                    subPhylum,
            Phylum,
//                    superPhylum,
//                    subKingdom,
            Kingdom,
            Domain,
        root,
        none

    };
    inline friend void write_binary(std::ostream& os, taxon_rank r) {
        write_binary(os, std::uint8_t(r));
    }
    inline friend void read_binary(std::istream& is, taxon_rank& r) {
        std::uint8_t n = 0;
        read_binary(is, n);
        r = taxon_rank(n);
    }

    //---------------------------------------------------------------
    static constexpr int num_ranks = static_cast<int>(taxon_rank::none);

    //---------------------------------------------------------------
    static taxon_rank
    next_main_rank(taxon_rank r) noexcept {
        switch(r) {
            case taxon_rank::Sequence:     return taxon_rank::Species;
            case taxon_rank::Form:         return taxon_rank::Species;
            case taxon_rank::Variety:      return taxon_rank::Species;
            case taxon_rank::subSpecies:   return taxon_rank::Species;
            case taxon_rank::Species:      return taxon_rank::Genus;
            case taxon_rank::subGenus:     return taxon_rank::Genus;
            case taxon_rank::Genus:        return taxon_rank::Family;
//            case taxon_rank::subTribe:     return taxon_rank::Family;
            case taxon_rank::Tribe:        return taxon_rank::Family;
//            case taxon_rank::subFamily:    return taxon_rank::Family;
            case taxon_rank::Family:       return taxon_rank::Order;
//            case taxon_rank::superFamily:  return taxon_rank::Order;
//            case taxon_rank::subOrder:     return taxon_rank::Order;
            case taxon_rank::Order:        return taxon_rank::Class;
//            case taxon_rank::superOrder:   return taxon_rank::Class;
//            case taxon_rank::subClass:     return taxon_rank::Class;
            case taxon_rank::Class:        return taxon_rank::Phylum;
//            case taxon_rank::superClass:   return taxon_rank::Phylum;
//            case taxon_rank::subPhylum:    return taxon_rank::Phylum;
            case taxon_rank::Phylum:       return taxon_rank::Kingdom;
//            case taxon_rank::superPhylum:  return taxon_rank::Kingdom;
//            case taxon_rank::subKingdom:   return taxon_rank::Kingdom;
            case taxon_rank::Kingdom:      return taxon_rank::Domain;
            case taxon_rank::Domain:       return taxon_rank::root;
            default:
            case taxon_rank::root:
            case taxon_rank::none: return taxon_rank::none;
        }
    }
    static taxon_rank
    prev_main_rank(taxon_rank r) noexcept {
        switch(r) {
            case taxon_rank::Sequence:     return taxon_rank::none;
            case taxon_rank::Form:         return taxon_rank::Sequence;
            case taxon_rank::Variety:      return taxon_rank::Sequence;
            case taxon_rank::subSpecies:   return taxon_rank::Sequence;
            case taxon_rank::Species:      return taxon_rank::Sequence;
            case taxon_rank::subGenus:     return taxon_rank::Species;
            case taxon_rank::Genus:        return taxon_rank::Species;
//            case taxon_rank::subTribe:     return taxon_rank::Genus;
            case taxon_rank::Tribe:        return taxon_rank::Genus;
//            case taxon_rank::subFamily:    return taxon_rank::Genus;
            case taxon_rank::Family:       return taxon_rank::Genus;
//            case taxon_rank::superFamily:  return taxon_rank::Family;
//            case taxon_rank::subOrder:     return taxon_rank::Family;
            case taxon_rank::Order:        return taxon_rank::Family;
//            case taxon_rank::superOrder:   return taxon_rank::Order;
//            case taxon_rank::subClass:     return taxon_rank::Order;
            case taxon_rank::Class:        return taxon_rank::Order;
//            case taxon_rank::superClass:   return taxon_rank::Class;
//            case taxon_rank::subPhylum:    return taxon_rank::Class;
            case taxon_rank::Phylum:       return taxon_rank::Class;
//            case taxon_rank::superPhylum:  return taxon_rank::Phylum;
//            case taxon_rank::subKingdom:   return taxon_rank::Phylum;
            case taxon_rank::Kingdom:      return taxon_rank::Phylum;
            case taxon_rank::Domain:       return taxon_rank::Kingdom;
            case taxon_rank::root:         return taxon_rank::Domain;
            default:
            case taxon_rank::none: return taxon_rank::none;
        }
    }


    //---------------------------------------------------------------
    inline friend taxon_rank& operator ++ (taxon_rank& r) {
        return int(r) < num_ranks ? (r = taxon_rank(int(r) + 1)) : r;
    }
    inline friend taxon_rank& operator -- (taxon_rank& r) {
        return int(r) > 1 ? (r = taxon_rank(int(r) - 1)) : r;
    }


    //---------------------------------------------------------------
    static taxon_rank
    rank_from_name(const std::string& name) noexcept {
        if(name == "sequence")      return taxon_rank::Sequence;
        if(name == "genome")        return taxon_rank::Sequence;
        if(name == "form")          return taxon_rank::Form;
        if(name == "forma")         return taxon_rank::Form;
        if(name == "variety")       return taxon_rank::Variety;
        if(name == "varietas")      return taxon_rank::Variety;
        if(name == "subspecies")    return taxon_rank::subSpecies;
        if(name == "species")       return taxon_rank::Species;
        if(name == "subgenus")      return taxon_rank::subGenus;
        if(name == "genus")         return taxon_rank::Genus;
//        if(name == "subtribe")      return taxon_rank::subTribe;
        if(name == "tribe")         return taxon_rank::Tribe;
//        if(name == "subfamily")     return taxon_rank::subFamily;
        if(name == "family")        return taxon_rank::Family;
//        if(name == "superfamily")   return taxon_rank::superFamily;
//        if(name == "suborder")      return taxon_rank::subOrder;
        if(name == "order")         return taxon_rank::Order;
//        if(name == "superorder")    return taxon_rank::superOrder;
//        if(name == "subclass")      return taxon_rank::subClass;
        if(name == "class")         return taxon_rank::Class;
//        if(name == "superclass")    return taxon_rank::superClass;
//        if(name == "subphylum")     return taxon_rank::subPhylum;
        if(name == "phylum")        return taxon_rank::Phylum;
        if(name == "division")      return taxon_rank::Phylum;
//        if(name == "superphylum")   return taxon_rank::superPhylum;
//        if(name == "subkingdom")    return taxon_rank::subKingdom;
        if(name == "kingdom")       return taxon_rank::Kingdom;
        if(name == "superkingdom")  return taxon_rank::Domain;
        if(name == "domain")        return taxon_rank::Domain;
        if(name == "root")          return taxon_rank::root;
        return taxon_rank::none;
    }


    //---------------------------------------------------------------
    static std::string
    rank_name(taxon_rank r) {
        switch(r) {
            case taxon_rank::Sequence:     return "sequence";
            case taxon_rank::Form:         return "form";
            case taxon_rank::Variety:      return "variety";
            case taxon_rank::subSpecies:   return "subspecies";
            case taxon_rank::Species:      return "species";
            case taxon_rank::subGenus:     return "subgenus";
            case taxon_rank::Genus:        return "genus";
//            case taxon_rank::subTribe:     return "subtribe";
            case taxon_rank::Tribe:        return "tribe";
//            case taxon_rank::subFamily:    return "subfamily";
            case taxon_rank::Family:       return "family";
//            case taxon_rank::superFamily:  return "superfamily";
//            case taxon_rank::subOrder:     return "suborder";
            case taxon_rank::Order:        return "order";
//            case taxon_rank::superOrder:   return "superorder";
//            case taxon_rank::subClass:     return "subclass";
            case taxon_rank::Class:        return "class";
//            case taxon_rank::superClass:   return "superclass";
//            case taxon_rank::subPhylum:    return "supphylum";
            case taxon_rank::Phylum:       return "phylum";
//            case taxon_rank::superPhylum:  return "superphylum";
//            case taxon_rank::subKingdom:   return "subkingdom";
            case taxon_rank::Kingdom:      return "kingdom";
            case taxon_rank::Domain:       return "domain";
            case taxon_rank::root:         return "root";
            default:
            case taxon_rank::none:         return "none";
        }
    }


    //-----------------------------------------------------
    /**
     * @brief taxonomic node information
     */
    class taxon {
    public:
        explicit
        taxon(taxon_id taxonId = 0,
              taxon_id parentId = 0,
              std::string taxonName = "--",
              taxon_rank taxonRank = taxon_rank::none)
        :
            id(taxonId), parent(parentId),
            rank(taxonRank),
            name(std::move(taxonName))
        {}

        inline friend bool
        operator < (const taxon& a, const taxon& b) noexcept {
            return a.id < b.id;
        }

        bool none() const noexcept {
            return id < 2;
        }

        std::string rank_name() const noexcept {
            return taxonomy::rank_name(rank);
        }

        //-----------------------------------------------------
        friend
        void read_binary(std::istream& is, taxon& t) {
            read_binary(is, t.id);
            read_binary(is, t.parent);
            read_binary(is, t.rank);
            read_binary(is, t.name);
        }

        //-----------------------------------------------------
        friend
        void write_binary(std::ostream& os, const taxon& t) {
            write_binary(os, t.id);
            write_binary(os, t.parent);
            write_binary(os, t.rank);
            write_binary(os, t.name);
        }

        //-----------------------------------------------------
        taxon_id id;
        taxon_id parent;
        taxon_rank rank;
        std::string name;
    };


    //-----------------------------------------------------
    using full_lineage   = std::vector<taxon_id>;
    using ranked_lineage = std::array<taxon_id,num_ranks>;


    //-------------------------------------------------------------------
    taxonomy():
        taxa_{}, noTaxon_{}
    {}


    //-------------------------------------------------------------------
    void
    clear() {
        taxa_.clear();
    }


    //-------------------------------------------------------------------
    void
    emplace(taxon_id taxonId,
            taxon_id parentId,
            const std::string& taxonName,
            const std::string& rankName)
    {
        emplace(taxonId, parentId, taxonName, rank_from_name(rankName));
    }
    //-----------------------------------------------------
    void
    emplace(taxon_id taxonId,
            taxon_id parentId = 0,
            const std::string& taxonName = "",
            taxon_rank rank = taxon_rank::none)
    {
        taxa_.emplace(taxonId, parentId, taxonName, rank);
    }
    //-----------------------------------------------------
    void
    insert(const taxon& t) {
        if(t.id > 0) taxa_.insert(t);
    }
    //-----------------------------------------------------
    void
    insert(taxon&& t) {
        if(t.id > 0) taxa_.insert(std::move(t));
    }


    //---------------------------------------------------------------
    /**
     * @brief if no rank set, assign sub<rank> of parent <rank>
     */
    void
    assign_ranks_to_all()
    {
        for(auto& tax : taxa_) {
            if(tax.rank == taxon_rank::none) {
                auto lr = lowest_rank(tax);
                if(lr != taxon_rank::none) {
                    if(lr > taxon_rank::subSpecies) --lr;
                    const_cast<taxon*>(&tax)->rank = lr;
                }
            }
        }
    }


    //---------------------------------------------------------------
    std::size_t
    taxon_count() const noexcept {
        return taxa_.size();
    }
    //-----------------------------------------------------
    bool
    empty() const noexcept {
        return taxa_.empty();
    }

    //-----------------------------------------------------
    const taxon&
    operator [] (taxon_id id) const {
        if(id < 1) return noTaxon_;
        auto it = taxa_.find(taxon{id});
        return (it != taxa_.end()) ? *it : noTaxon_;
    }


    //---------------------------------------------------------------
    const taxon&
    lca(const taxon& a, const taxon& b) const
    {
        return lca(a.id, b.id);
    }

    //-----------------------------------------------------
    const taxon&
    lca(taxon_id a, taxon_id b) const
    {
        return operator[](lca_id(lineage(a), lineage(b) ));
    }

    //---------------------------------------------------------------
    const taxon&
    ranked_lca(const taxon& a, const taxon& b) const
    {
        return ranked_lca(a.id, b.id);
    }
    //-----------------------------------------------------
    const taxon&
    ranked_lca(taxon_id a, taxon_id b) const
    {
        return operator[](ranked_lca_id(ranks(a),
                                        ranks(b) ));
    }


    //---------------------------------------------------------------
    static taxon_id
    ranked_lca_id(const ranked_lineage& lina,
                  const ranked_lineage& linb)
    {
        for(int i = 0; i < int(taxon_rank::root); ++i) {
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
    ranked_lineage
    ranks(const taxon& tax) const {
        return ranks(tax.id);
    }

    //-----------------------------------------------------
    ranked_lineage
    ranks(taxon_id id) const
    {
//        std::cout << "      main rank lineage:" << std::endl;
        auto lin = ranked_lineage{};
        for(auto& x : lin) x = 0;

        while(id) {
            auto it = taxa_.find(taxon{id});
            if(it != taxa_.end()) {

                if(it->rank != taxon_rank::none) {
                    lin[static_cast<int>(it->rank)] = it->id;

//                    std::cout << "        " << id << " (" << it->rank_name()
//                              << " (" << int(it->rank) << ") "
//                              << it->name << ")" << std::endl;

                }

                if(it->parent != id) {
                    id = it->parent;
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
        return lineage(tax.id);
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
                if(it->parent != id) {
                    id = it->parent;
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


    //---------------------------------------------------------------
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
    taxon_rank
    lowest_rank(taxon_id id) const
    {
       while(id) {
            auto it = taxa_.find(taxon{id});
            if(it != taxa_.end()) {
                if(it->rank != taxon_rank::none) {
                    return it->rank;
                }
                if(it->parent != id) {
                    id = it->parent;
                } else {
                    id = 0;
                }
            } else {
                id = 0;
            }
        }
        return taxon_rank::none;
    }
    //-----------------------------------------------------
    taxon_rank
    lowest_rank(const taxon& tax) const {
        return lowest_rank(tax.id);
    }

    //---------------------------------------------------------------
    std::set<taxon> taxa_;
    taxon noTaxon_;
};






/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
class rank_statistics
{
public:
    //---------------------------------------------------------------
    using count_t = std::uint_least64_t;
    using taxon_rank = taxonomy::taxon_rank;

    //---------------------------------------------------------------
    rank_statistics() noexcept :
        classified_{}, known_{}, correct_{}, coverage_{}
    {
        for(auto& x : classified_) x = 0;
        for(auto& x : known_)      x = 0;
        for(auto& x : correct_)    x = 0;
    }

    //---------------------------------------------------------------
    void assign(taxon_rank assigned) noexcept
    {
        if(assigned == taxon_rank::none) {
            ++classified_[int(taxon_rank::none)];
        }
        else {
            for(auto i = int(assigned); i <= int(taxon_rank::root); ++i) {
                ++classified_[i];
            }
        }
    }

    //---------------------------------------------------------------
    void assign_known_correct(taxon_rank assigned,
                              taxon_rank known,
                              taxon_rank correct) noexcept
    {
        assign(assigned);

        if(known == taxon_rank::none) {
            ++known_[int(taxon_rank::none)];
        }
        else {
            for(auto i = int(known); i <= int(taxon_rank::root); ++i) {
                ++known_[i];
            }

            if(correct == taxon_rank::none) {
                ++correct_[int(taxon_rank::none)];
            }
            else {
                for(auto i = int(correct); i <= int(taxon_rank::root); ++i) {
                    ++correct_[i];
                }
            }
        }
    }


    //---------------------------------------------------------------
    void count_coverage_true_pos(taxon_rank r) {
        coverage_[int(r)].count_true_pos();
    }
    void count_coverage_false_pos(taxon_rank r) {
        coverage_[int(r)].count_false_pos();
    }
    void count_coverage_true_neg(taxon_rank r) {
        coverage_[int(r)].count_true_neg();
    }
    void count_coverage_false_neg(taxon_rank r) {
        coverage_[int(r)].count_false_neg();
    }

    //-----------------------------------------------------
    const confusion_statistics&
    coverage(taxon_rank r) const noexcept {
        return coverage_[int(r)];
    }


    //---------------------------------------------------------------
    count_t classified() const noexcept {
        return classified_[int(taxon_rank::root)];
    }
    count_t classified(taxon_rank r) const noexcept {
        return classified_[int(r)];
    }
    count_t unclassified() const noexcept {
        return classified_[int(taxon_rank::none)];
    }

    count_t total() const noexcept {
        return classified() + unclassified();
    }

    count_t known() const noexcept {
        return known_[int(taxon_rank::root)];
    }
    count_t known(taxon_rank r) const noexcept {
        return known_[int(r)];
    }
    count_t unknown() const noexcept {
        return known_[int(taxon_rank::none)];
    }

    count_t correct(taxon_rank r) const noexcept {
        return correct_[int(r)];
    }
    count_t correct() const noexcept {
        return correct_[int(taxon_rank::root)];
    }


    //---------------------------------------------------------------
    double known_rate(taxon_rank r) const noexcept {
        return total() > 0 ?  known(r) / double(total()) : 0;
    }
    double known_rate() const noexcept {
        return total() > 0 ?  known() / double(total()) : 0;
    }
    double unknown_rate() const noexcept {
        return total() > 0 ?  unknown() / double(total()) : 0;
    }

    double classification_rate(taxon_rank r) const noexcept {
        return total() > 0 ? classified(r) / double(total()) : 0;
    }
    double unclassified_rate() const noexcept {
        return total() > 0 ? unclassified() / double(total()) : 0;
    }

    double sensitivity(taxon_rank r) const noexcept {
        return known(r) > 0 ? correct(r) / double(known(r)) : 0;
    }
    double precision(taxon_rank r) const noexcept {
        return classified(r) > 0 ? correct(r) / double(classified(r)) : 0;
    }

private:
    //---------------------------------------------------------------
    count_t classified_[taxonomy::num_ranks+1];
    count_t known_[taxonomy::num_ranks+1];
    count_t correct_[taxonomy::num_ranks+1];
    confusion_statistics coverage_[taxonomy::num_ranks+1];
};



} // namespace mc


#endif

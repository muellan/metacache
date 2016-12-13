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

#ifndef MC_TAXONOMY_H_
#define MC_TAXONOMY_H_

#include <array>
#include <vector>
#include <set>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <atomic>

#include "confusion.h"
#include "io_serialize.h"


namespace mc {


/*****************************************************************************
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
    using taxon_id = std::uint64_t;

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
        if(name == "suborder")      return rank::subOrder;
        if(name == "order")         return rank::Order;
        if(name == "subclass")      return rank::subClass;
        if(name == "class")         return rank::Class;
        if(name == "subphylum")     return rank::subPhylum;
        if(name == "phylum")        return rank::Phylum;
        if(name == "division")      return rank::Phylum;
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
            case rank::subPhylum:    return "supphylum";
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
        using rank_t = taxonomy::rank;
    public:
        explicit
        taxon(taxon_id taxonId = 0,
              taxon_id parentId = 0,
              std::string taxonName = "--",
              rank_t rk = rank_t::none)
        :
            id(taxonId), parent(parentId),
            rank(rk),
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
        rank_t rank;
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
            rank rank = rank::none)
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
    rank_all_unranked()
    {
        for(auto& tax : taxa_) {
            if(tax.rank == rank::none) {
                auto lr = lowest_rank(tax);
                if(lr != rank::none) {
                    if(lr > rank::subSpecies) --lr;
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
        return operator[](ranked_lca_id(ranks(a), ranks(b) ));
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


    //---------------------------------------------------------------
    /**
     * @param id
     * @return
     */
    ranked_lineage
    ranks(taxon_id id) const
    {
        auto lin = ranked_lineage{};
        for(auto& x : lin) x = 0;

        while(id) {
            auto it = taxa_.find(taxon{id});
            if(it != taxa_.end()) {
                if(it->rank != rank::none) {
                    lin[static_cast<int>(it->rank)] = it->id;
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

    //-----------------------------------------------------
    ranked_lineage
    ranks(const taxon& tax) const {
        return ranks(tax.id);
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
                if(it->rank != rank::none) {
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
        return rank::none;
    }
    //-----------------------------------------------------
    rank
    lowest_rank(const taxon& tax) const {
        return lowest_rank(tax.id);
    }

    //---------------------------------------------------------------
    std::set<taxon> taxa_;
    taxon noTaxon_;
};






/*****************************************************************************
 *
 * @brief useful for tracking classification properties (precision, ...)
 *
 *****************************************************************************/
class rank_statistics
{
public:
    //---------------------------------------------------------------
    using count_t = std::uint_least64_t;
    using rank = taxonomy::rank;

    //---------------------------------------------------------------
    rank_statistics() noexcept :
        assigned_{}, known_{}, correct_{}, wrong_{}, coverage_{}
    {
        for(auto& x : assigned_) x = 0;
        for(auto& x : known_)    x = 0;
        for(auto& x : correct_)  x = 0;
        for(auto& x : wrong_)    x = 0;
    }

    //---------------------------------------------------------------
    /**
     * @brief    counts one assignment
     * @details  concurrency-safe
     */
    void assign(rank assigned) noexcept
    {
        if(assigned == rank::none) {
            ++assigned_[int(rank::none)];
        }
        else {
            for(rank r = assigned; r <= rank::root; ++r) ++assigned_[int(r)];
        }
    }

    //---------------------------------------------------------------
    /**
     * @brief   counts one assignment including ground truth and correctness
     *          assessment
     *
     * @details concurrency-safe
     *
     * @param assigned : lowest rank of current assignment
     * @param known    : lowest rank for which ground truth was known
     * @param correct  : lowest rank for which current assignment is correct
     */
    void assign_known_correct(rank assigned, rank known, rank correct) noexcept
    {
        assign(assigned);

        //plausibility check
        if(correct < assigned) correct = assigned;
        if(correct < known)    correct = known;

        //if ground truth known -> count correct and wrong assignments
        if(known == rank::none) {
            ++known_[int(rank::none)];
        }
        else {
            for(rank r = known; r <= rank::root; ++r) ++known_[int(r)];

            if(correct == rank::none) {
                ++correct_[int(rank::none)];
            }
            else {
                for(rank r = correct; r <= rank::root; ++r) ++correct_[int(r)];
            }
            //if ranks above and including the current rank of assignment
            //are wrong => levels below of current assignment must be wrong, too
            if(correct > assigned) {
                for(rank r = known; r < correct; ++r) {
                    ++wrong_[int(r)];
                }
            }
        }
    }


    //---------------------------------------------------------------
    /// @details concurrency-safe
    void count_coverage_true_pos(rank r) {
        coverage_[int(r)].count_true_pos();
    }
    /// @details concurrency-safe
    void count_coverage_false_pos(rank r) {
        coverage_[int(r)].count_false_pos();
    }
    /// @details concurrency-safe
    void count_coverage_true_neg(rank r) {
        coverage_[int(r)].count_true_neg();
    }
    /// @details concurrency-safe
    void count_coverage_false_neg(rank r) {
        coverage_[int(r)].count_false_neg();
    }

    //-----------------------------------------------------
    const confusion_statistics&
    coverage(rank r) const noexcept {
        return coverage_[int(r)];
    }


    count_t assigned() const noexcept {
        return assigned_[int(rank::root)];
    }
    /**
     * @brief number of assignments on a taxonomix rank (and above)
     */
    count_t assigned(rank r) const noexcept {
        return assigned_[int(r)];
    }
    count_t unassigned() const noexcept {
        return assigned_[int(rank::none)];
    }

    //---------------------------------------------------------------
    count_t total() const noexcept {
        return assigned() + unassigned();
    }

    count_t known() const noexcept {
        return known_[int(rank::root)];
    }
    /**
     * @brief number of cases with ground truth known on ranks >= r
     */
    count_t known(rank r) const noexcept {
        return known_[int(r)];
    }
    count_t unknown() const noexcept {
        return known_[int(rank::none)];
    }

    /**
     * @brief number of known correct assignments on ranks >= r
     */
    count_t correct(rank r) const noexcept {
        return correct_[int(r)];
    }
    count_t correct() const noexcept {
        return correct_[int(rank::root)];
    }

    /**
     * @brief number of known wrong assignments on ranks <= r
     */
    count_t wrong(rank r) const noexcept {
        return wrong_[int(r)];
    }
    count_t wrong() const noexcept {
        return wrong_[int(rank::root)];
    }


    //---------------------------------------------------------------
    double known_rate(rank r) const noexcept {
        return total() > 0 ?  known(r) / double(total()) : 0;
    }
    double known_rate() const noexcept {
        return total() > 0 ?  known() / double(total()) : 0;
    }
    double unknown_rate() const noexcept {
        return total() > 0 ?  unknown() / double(total()) : 0;
    }
    double classification_rate(rank r) const noexcept {
        return total() > 0 ? assigned(r) / double(total()) : 0;
    }
    double unclassified_rate() const noexcept {
        return total() > 0 ? unassigned() / double(total()) : 0;
    }

    double sensitivity(rank r) const noexcept {
        return known(r) > 0 ? correct(r) / double(known(r)) : 0;
    }
    double precision(rank r) const noexcept {
        //note that in general tot != assigned(r) and tot != known(r)
        double tot = correct(r) + wrong(r);
        return tot > 0 ? correct(r) / tot : 0;
    }

private:
    //---------------------------------------------------------------
    std::atomic<count_t> assigned_[taxonomy::num_ranks+1];
    std::atomic<count_t> known_[taxonomy::num_ranks+1];
    std::atomic<count_t> correct_[taxonomy::num_ranks+1];
    std::atomic<count_t> wrong_[taxonomy::num_ranks+1];
    confusion_statistics coverage_[taxonomy::num_ranks+1];
};



} // namespace mc


#endif

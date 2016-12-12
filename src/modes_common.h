/*****************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * version 0.1
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

#ifndef MC_MODES_COMMON_H_
#define MC_MODES_COMMON_H_


#include "args_parser.h"

//contains all important typedefs
#include "config.h"


namespace mc {


/*****************************************************************************
 *
 * @brief convenience type aliases from database
 *
 *****************************************************************************/
using taxon_rank     = database::taxon_rank;
using taxon          = database::taxon;
using taxon_id       = database::taxon_id;
using ranked_lineage = database::ranked_lineage;
using match_result   = database::match_result;



/*****************************************************************************
 *
 * @brief needed because the NCBI does not define taxon ids on sequence level
 *
 *****************************************************************************/
struct classification
{

    explicit constexpr
    classification(const taxon* tax = nullptr) noexcept:
        tid_{database::invalid_target_id()}, tax_{tax}
    {}

    explicit constexpr
    classification(target_id tid, const taxon* tax = nullptr) noexcept:
        tid_{tid}, tax_{tax}
    {}

    bool has_taxon() const noexcept {
        return tax_;
    }
    bool sequence_level() const noexcept {
        return tid_ != database::invalid_target_id();
    }
    bool none() const noexcept {
        return !sequence_level() && !has_taxon();
    }

    target_id target() const noexcept {
        return tid_;
    }
    const taxon& tax() const noexcept { return *tax_; }

    taxon_rank rank() const noexcept {
        return sequence_level() ? taxon_rank::Sequence
                                : (tax_ ? tax_->rank : taxon_rank::none);
    }

private:
    target_id tid_;
    const taxon* tax_;
};



/*****************************************************************************
 *
 * @brief how taxon formatting will be done
 *
 *****************************************************************************/
enum class taxon_print_mode : unsigned char {
    name_only = 1, id_only = 2, id_name = 3
};



/*****************************************************************************
 *
 * @brief returns query taxon (ground truth for precision tests)
 *
 *****************************************************************************/
classification
ground_truth(const database&, const std::string& sequenceHeader);



/*****************************************************************************
 *
 * @brief returns lowest common rank of 2 classifications
 *
 *****************************************************************************/
taxon_rank
lowest_common_rank(const database&,
                   const classification&, const classification&);



/*****************************************************************************
 *
 * @brief print ranked lineage
 *
 *****************************************************************************/
void show_ranks(std::ostream&,
                const database&,
                const ranked_lineage&,
                taxon_print_mode = taxon_print_mode::name_only,
                taxon_rank lowest  = taxon_rank::Species,
                taxon_rank highest = taxon_rank::Domain);



/*****************************************************************************
 *
 * @brief print ranked lineage of target
 *
 *****************************************************************************/
void show_ranks_of_target(std::ostream&,
                          const database&,
                          target_id,
                          taxon_print_mode = taxon_print_mode::name_only,
                          taxon_rank lowest  = taxon_rank::Sequence,
                          taxon_rank highest = taxon_rank::Domain);



/*****************************************************************************
 *
 * @brief print per-rank classification statistics
 *
 *****************************************************************************/
void show_per_rank_statistics(std::ostream&,
                              const rank_statistics&,
                              const std::string& prefix = "");



/*****************************************************************************
 *
 * @brief add difference between result and truth to statistics
 *
 *****************************************************************************/
void update_coverage_statistics(
    const database&,
    const classification& result, const classification& truth,
    rank_statistics&);



/*****************************************************************************
 *
 * @brief  print target.window hit matches
 *
 *****************************************************************************/
void show_matches(std::ostream&,
                  const database&,
                  const match_result&,
                  taxon_rank lowest = taxon_rank::Sequence);


//-------------------------------------------------------------------
template<int n>
void show_matches(std::ostream& os,
    const database& db,
    const matches_in_contiguous_window_range_top<n,target_id>& cand,
    taxon_rank lowest = taxon_rank::Sequence)
{
    if(lowest == taxon_rank::Sequence) {
        if(cand.hits(0) > 0) {
            os  << db.sequence_id_of_target(cand.target_id(0))
                << ':' << cand.hits(0);
        }

        for(int i = 1; i < n && cand.hits(i) > 0; ++i) {
            os  << ',' << db.sequence_id_of_target(cand.target_id(i))
                << ':' << cand.hits(i);
            ++i;
        }
    }
    else {
        if(cand.hits(0) > 0) {
            auto taxid = db.ranks_of_target(
                         cand.target_id(0))[int(lowest)];

            os << taxid << ':' << cand.hits(0);
        }

        for(int i = 1; i < n && cand.hits(i) > 0; ++i) {
            auto taxid = db.ranks_of_target(
                         cand.target_id(i))[int(lowest)];

            os << ',' << taxid << ':' << cand.hits(i);
            ++i;
        }
    }
}



} // namespace  m

#endif

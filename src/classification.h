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

#ifndef MC_CLASSIFICATION_H_
#define MC_CLASSIFICATION_H_


#include "config.h"


namespace mc {


/*****************************************************************************
 *
 * @brief needed because taxon ids are not defined at sequence level
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



} // namespace  m

#endif

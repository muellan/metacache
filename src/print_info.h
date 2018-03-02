/*****************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2018 André Müller (muellan@uni-mainz.de)
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

#ifndef MC_PRINT_INFO_H_
#define MC_PRINT_INFO_H_

#include "config.h"


namespace mc {


/*****************************************************************************
 *
 * @brief how taxon formatting will be done
 *
 *****************************************************************************/
enum class taxon_print_mode : unsigned char {
    rank_id, rank_name, rank_name_id
};


/*****************************************************************************
 *
 * @brief print database information
 *
 *****************************************************************************/
void show_info(std::ostream&, const database&, const taxon&);



/*****************************************************************************
 *
 * @brief print ranked lineage
 *
 *****************************************************************************/
void show_ranks(std::ostream&,
                const ranked_lineage&,
                taxon_print_mode = taxon_print_mode::rank_name,
                taxon_rank lowest  = taxon_rank::Sequence,
                taxon_rank highest = taxon_rank::Domain,
                const std::string& columnSeparator = "");


} // namespace mc


#endif

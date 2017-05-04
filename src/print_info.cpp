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

#include "print_info.h"


namespace mc {


//-------------------------------------------------------------------
void show_ranks(
    std::ostream& os,
    const database& db,
    const ranked_lineage& lineage,
    taxon_print_mode mode, taxon_rank lowest, taxon_rank highest)
{
    //one rank only
    if(lowest == highest) {
        auto taxid = lineage[int(lowest)];
        os << taxonomy::rank_name(lowest) << ':';
        if(mode != taxon_print_mode::id_only) {
            if(taxid > 1)
                os << db.taxon_with_id(taxid).name();
            else
                os << "n/a";
            if(mode != taxon_print_mode::name_only)
                os << "(" << taxid << ")";
        }
        else {
            os << taxid;
        }
    }
    //range of ranks
    else {
        for(auto r = lowest; r <= highest; ++r) {
            auto taxid = lineage[int(r)];
            if(taxid > 1) {
                auto&& taxon = db.taxon_with_id(taxid);
                if(taxon.rank() >= lowest && taxon.rank() <= highest) {
                    os << taxon.rank_name() << ':';
                    if(mode != taxon_print_mode::id_only) {
                        os << taxon.name();
                        if(mode != taxon_print_mode::name_only) {
                            os << "(" << taxon.id() << ")";
                        }
                    }
                    else {
                        os << taxon.id();
                    }
                }
                if(r < highest) os << ',';
            }
            else if (r == lowest && taxid == 0) {
                os << taxonomy::rank_name(lowest) << ':';
                if(mode != taxon_print_mode::id_only) {
                    os << "n/a";
                    if(mode != taxon_print_mode::name_only)
                        os << "(" << taxid << ")";
                }
                else {
                    os << taxid;
                }
            }
        }
    }
}



//-------------------------------------------------------------------
void show_ranks_of_target(
    std::ostream& os,
    const database& db,
    target_id tid,
    taxon_print_mode mode, taxon_rank lowest, taxon_rank highest)
{
    //since targets don't have their own taxonId, print their sequence id
    if(lowest == taxon_rank::Sequence) {
        if(mode != taxon_print_mode::id_only) {
            os << "sequence:" << db.sequence_id_of_target(tid);
        } else {
            os << db.sequence_id_of_target(tid);
        }
    }
    if(highest == taxon_rank::Sequence) return;

    if(lowest == taxon_rank::Sequence) os << ',';

    show_ranks(os, db, db.ranks_of_target(tid),
                        mode, lowest, highest);
}


} // namespace mc

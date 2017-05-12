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
void show_info(std::ostream& os, const database& db, const taxon& tax)
{
    os  << "Target " << tax.name() << "):\n"
        << "    source:     "
        << tax.source().filename << " / " << tax.source().index;

    for(const taxon* t : db.ranks(tax)) {
        if(t) {
            auto rn = std::string(t->rank_name()) + ":";
            rn.resize(12, ' ');
            os << "\n    " << rn << "(" << t->id() << ") " << t->name();
        }
    }
    os << '\n';
}



//-------------------------------------------------------------------
void show_ranks(
    std::ostream& os,
    const ranked_lineage& lineage,
    taxon_print_mode mode, taxon_rank lowest, taxon_rank highest)
{
    if(lowest == taxon_rank::none) return;
    if(highest == taxon_rank::none) highest = taxon_rank::root;

    //one rank only
    if(lowest == highest) {
        const taxon* tax = lineage[int(lowest)];
        if(tax) {
            os << taxonomy::rank_name(lowest) << ':';
            if(mode != taxon_print_mode::id_only) {
                if(tax) {
                    os << tax->name();
                    if(mode != taxon_print_mode::name_only)
                        os << "(" << tax->id() << ")";
                } else {
                    os << "n/a";
                    if(mode != taxon_print_mode::name_only)
                        os << "(" << taxonomy::none_id() << ")";
                }
            } else {
                os << (tax ? tax->id() : taxonomy::none_id());
            }
        }
    }
    //range of ranks
    else {
        for(auto r = lowest; r <= highest; ++r) {
            const taxon* tax = lineage[int(r)];
            if(tax) {
                if(tax->rank() >= lowest && tax->rank() <= highest) {
                    os << tax->rank_name() << ':';
                    if(mode != taxon_print_mode::id_only) {
                        os << tax->name();
                        if(mode != taxon_print_mode::name_only) {
                            os << "(" << tax->id() << ")";
                        }
                    } else {
                        os << tax->id();
                    }
                }
                if(r < highest) os << ',';
            }
            else if(r == lowest) {
                os << taxonomy::rank_name(lowest) << ':';
                if(mode != taxon_print_mode::id_only) {
                    os << "n/a";
                    if(mode != taxon_print_mode::name_only)
                        os << "(" << taxonomy::none_id() << ")";
                } else {
                    os << taxonomy::none_id();
                }
            }
        }
    }
}



} // namespace mc

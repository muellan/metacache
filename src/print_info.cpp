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
        << tax.source().filename << " / " << tax.source().index
        << "\n    length:     " << tax.source().windows << " windows";

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
void show_taxon(std::ostream& os,
                taxon_print_mode mode,
                const taxon* tax,
                const std::string& separator)
{
    if(!tax) return;
    if(separator.empty()) {
        switch(mode) {
            default:
            case taxon_print_mode::rank_name:
                os << taxonomy::rank_name(tax->rank()) << ':' << tax->name();
                break;
            case taxon_print_mode::rank_id:
                os << taxonomy::rank_name(tax->rank()) << ':' << tax->id();
                break;
            case taxon_print_mode::rank_name_id:
                os << taxonomy::rank_name(tax->rank()) << ':' << tax->name()
                   << '(' << tax->id() << ')';
                break;
        }
    } else {
        switch(mode) {
            default:
            case taxon_print_mode::rank_name:
                os << taxonomy::rank_name(tax->rank()) << separator << tax->name();
                break;
            case taxon_print_mode::rank_id:
                os << taxonomy::rank_name(tax->rank()) << separator << tax->id();
                break;
            case taxon_print_mode::rank_name_id:
                os << taxonomy::rank_name(tax->rank())
                   << separator << tax->name() << separator << tax->id();
                break;
        }
    }
}



//-------------------------------------------------------------------
void show_no_taxon(std::ostream& os,
                   taxon_print_mode mode,
                   taxon_rank rank,
                   const std::string& separator)
{
    if(separator.empty()) {
        switch(mode) {
            default:
            case taxon_print_mode::rank_name:
                os << taxonomy::rank_name(rank) << ":n/a";
                break;
            case taxon_print_mode::rank_id:
                os << taxonomy::rank_name(rank)
                << '(' << taxonomy::none_id() << ')';
                break;
            case taxon_print_mode::rank_name_id:
                os << taxonomy::rank_name(rank) << ":n/a("
                << taxonomy::none_id() << ')';
                break;
        }
    } else {
        switch(mode) {
            default:
            case taxon_print_mode::rank_name:
                os << taxonomy::rank_name(rank) << separator << "n/a";
                break;
            case taxon_print_mode::rank_id:
                os << taxonomy::rank_name(rank) << separator
                << taxonomy::none_id();
                break;
            case taxon_print_mode::rank_name_id:
                os << taxonomy::rank_name(rank) << separator << "n/a"
                << separator << taxonomy::none_id();
                break;
        }
    }
}



//-------------------------------------------------------------------
void show_ranks(
    std::ostream& os,
    const ranked_lineage& lineage,
    taxon_print_mode mode, taxon_rank lowest, taxon_rank highest,
    const std::string& separator)
{
    if(lowest == taxon_rank::none) return;
    if(highest == taxon_rank::none) highest = taxon_rank::root;

    //one rank only
    if(lowest == highest) {
        show_taxon(os, mode, lineage[int(lowest)], separator);
    }
    //range of ranks
    else {
        for(auto r = lowest; r <= highest; ++r) {
            const taxon* tax = lineage[int(r)];

            if(tax) {
                show_taxon(os, mode, lineage[int(r)], separator);
            }
            else {
                show_no_taxon(os, mode, lowest, separator);
            }
            if(r < highest) {
                if(separator.empty()) {
                    os << ',';
                } else {
                    os << separator;
                }
            }

        }
    }

}



} // namespace mc

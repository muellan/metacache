/*****************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2018 André Müller (muellan@uni-mainz.de)
 *                       & Robin Kobus  (rkobus@uni-mainz.de)
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

#include "printing.h"


namespace mc {



//-------------------------------------------------------------------
void show_target_info(std::ostream& os, const database& db, const taxon& tax)
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
                const taxon* tax, taxon_print_mode mode,
                const std::string& separator = "")
{
    if(!tax) return;
    if(separator.empty()) {
        switch(mode) {
            default:
            case taxon_print_mode::rank_name:
                os << taxonomy::rank_name(tax->rank()) << ':';
                [[fallthrough]] // fall through
            case taxon_print_mode::name:
                os << tax->name();
                break;
            case taxon_print_mode::rank_id:
                os << taxonomy::rank_name(tax->rank()) << ':';
                [[fallthrough]] // fall through
            case taxon_print_mode::id:
                os << tax->id();
                break;
            case taxon_print_mode::rank_name_id:
                os << taxonomy::rank_name(tax->rank()) << ':';
                [[fallthrough]] // fall through
            case taxon_print_mode::name_id:
                os << tax->name() << '(' << tax->id() << ')';
                break;
        }
    } else {
        switch(mode) {
            default:
            case taxon_print_mode::rank_name:
                os << taxonomy::rank_name(tax->rank()) << separator;
                [[fallthrough]] // fall through
            case taxon_print_mode::name:
                os << tax->name();
                break;
            case taxon_print_mode::rank_id:
                os << taxonomy::rank_name(tax->rank()) << separator;
                [[fallthrough]] // fall through
            case taxon_print_mode::id:
                os << tax->id();
                break;
            case taxon_print_mode::rank_name_id:
                os << taxonomy::rank_name(tax->rank()) << separator;
                [[fallthrough]] // fall through
            case taxon_print_mode::name_id:
                os << tax->name() << separator << tax->id();
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
                os << taxonomy::rank_name(rank) << ':';
                [[fallthrough]] // fall through
            case taxon_print_mode::name:
                os << "--";
                break;
            case taxon_print_mode::rank_id:
                os << taxonomy::rank_name(rank) << ':';
                [[fallthrough]] // fall through
            case taxon_print_mode::id:
                os << taxonomy::none_id();
                break;
            case taxon_print_mode::rank_name_id:
                os << taxonomy::rank_name(rank) << ':';
                [[fallthrough]] // fall through
            case taxon_print_mode::name_id:
                os << "--(" << taxonomy::none_id() << ')';
                break;
        }
    } else {
        switch(mode) {
            default:
            case taxon_print_mode::rank_name:
                os << taxonomy::rank_name(rank) << separator;
                [[fallthrough]] // fall through
            case taxon_print_mode::name:
                os << "--";
                break;
            case taxon_print_mode::rank_id:
                os << taxonomy::rank_name(rank) << separator;
                [[fallthrough]] // fall through
            case taxon_print_mode::id:
                os << taxonomy::none_id();
                break;
            case taxon_print_mode::rank_name_id:
                os << taxonomy::rank_name(rank) << separator;
                [[fallthrough]] // fall through
            case taxon_print_mode::name_id:
                os << "--" << separator << taxonomy::none_id();
                break;
        }
    }
}



//-------------------------------------------------------------------
void show_lineage(
    std::ostream& os,
    const ranked_lineage& lineage,
    taxon_print_mode mode, taxon_rank lowest, taxon_rank highest,
    const std::string& separator)
{
    if(lowest == taxon_rank::none) return;
    if(highest == taxon_rank::none) highest = taxon_rank::root;

    //one rank only
    if(lowest == highest) {
        show_taxon(os, lineage[int(lowest)], mode, separator);
    }
    //range of ranks
    else {
        for(auto r = lowest; r <= highest; ++r) {
            const taxon* tax = lineage[int(r)];
            if(tax) {
                show_taxon(os, tax, mode, separator);
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

//-------------------------------------------------------------------
void show_blank_lineage(std::ostream& os,
                        taxon_print_mode mode,
                        taxon_rank lowest, taxon_rank highest,
                        const std::string& separator)
{
    for(auto r = lowest; r <= highest; ++r) {
        switch(mode) {
            default:
            case taxon_print_mode::rank_name:
                os << "--" << separator;
                [[fallthrough]] // fall through
            case taxon_print_mode::name:
                os << "--";
                break;
            case taxon_print_mode::rank_id:
                os << "--" << separator;
                [[fallthrough]] // fall through
            case taxon_print_mode::id:
                os << taxonomy::none_id();
                break;
            case taxon_print_mode::rank_name_id:
                os << "--" << separator;
                [[fallthrough]] // fall through
            case taxon_print_mode::name_id:
                os << "--" << separator << taxonomy::none_id();
                break;
        }
        if(r < highest) os << separator;
    }
}



//-------------------------------------------------------------------
void show_taxon(std::ostream& os,
                const database& db,
                const classification_output_options& out,
                const taxon* tax)
{
    if(!tax || tax->rank() > out.highestRank) {
        if(out.separateTaxaInfo) {
            const auto rmax = out.showLineage ? out.highestRank : out.lowestRank;

            show_blank_lineage(os, out.showTaxaAs,
                               out.lowestRank, rmax, out.separator);
        }
        else {
            os << "--";
        }
    }
    else {
        const auto rmin = out.lowestRank < tax->rank()
                          ? tax->rank() : out.lowestRank;

        const auto rmax = out.showLineage ? out.highestRank : rmin;

        if(out.separateTaxaInfo) {
            show_lineage(os, db.ranks(tax), out.showTaxaAs,
                         rmin, rmax, out.separator);
        }
        else {
            show_lineage(os, db.ranks(tax), out.showTaxaAs,
                         rmin, rmax, "");
        }
    }
}



//-------------------------------------------------------------------
void show_matches(std::ostream& os,
                  const database& db,
                  const classification_candidates& cand,
                  taxon_rank lowest)
{
    if(lowest == taxon_rank::Sequence) {
        for(int i = 0; i < cand.size() && cand[i].hits > 0; ++i) {
            if(i > 0) os << ',';
            if(cand[i].tax) os << cand[i].tax->name() << ':' << cand[i].hits;
        }
    }
    else {
        for(int i = 0; i < cand.size() && cand[i].hits > 0; ++i) {
            if(i > 0) os << ',';
            const taxon* tax = db.ancestor(cand[i].tax,lowest);
            if(tax) {
                os << tax->id() << ':' << cand[i].hits;
            } else {
                os << cand[i].tax->name() << ':' << int(cand[i].hits) << ',';
            }
        }
    }
}



//-------------------------------------------------------------------
void show_matches(std::ostream& os,
                  const database& db,
                  const matches_per_location& matches,
                  taxon_rank lowest)
{
    if(matches.empty()) return;

    if(lowest == taxon_rank::Sequence) {
        for(const auto& r : matches) {
            const taxon* tax = r.loc.tax;
            if(tax) os << tax->name()
                       << '/' << int(r.loc.win)
                       << ':' << int(r.hits) << ',';
        }
    }
    else {
        for(const auto& r : matches) {
            const taxon* tax = db.ancestor(r.loc.tax, lowest);
            if(tax) {
                os << tax->name() << ':' << int(r.hits) << ',';
            } else {
                os << r.loc.tax->name() << ':' << int(r.hits) << ',';
            }
        }
    }
}




//-------------------------------------------------------------------
void show_candidate_ranges(std::ostream& os,
                           const database& db,
                           const classification_candidates& cand)
{
    const auto w = db.target_window_stride();

    for(int i = 0; i < cand.size(); ++i) {
        os << '[' << (w * cand[i].pos.beg)
           << ',' << (w * cand[i].pos.end + db.target_window_size()) << "] ";
    }
}



//-------------------------------------------------------------------
void show_matches_per_targets(std::ostream& os,
                              const matches_per_target& tgtMatches,
                              const classification_output_options& out)
{
    os  << out.comment << " genome " << out.separator
        << "queryid/window:hits/window:hits/...,queryid/...\n";

    for(const auto& mapping : tgtMatches) {
        show_taxon(os, mapping.first, out.showTaxaAs);
        os << out.separator;

        bool first = true;
        for(const auto& candidate : mapping.second) {
            if(first) first = false; else os << ',';
            os << candidate.qeryid;
            // windows with hit count
            for(const auto& match : candidate.matches) {
                os << '/' << match.win << ':' << match.hits;
            }
        }
    }
}



//-------------------------------------------------------------------
void show_taxon_statistics(std::ostream& os,
                                    const classification_statistics& stats,
                                    const std::string& prefix)
{
    constexpr taxon_rank ranks[] {
          taxon_rank::Sequence,
          taxon_rank::subSpecies,
          taxon_rank::Species,    taxon_rank::Genus,
          taxon_rank::Family,     taxon_rank::Order,
          taxon_rank::Class,      taxon_rank::Phylum,
          taxon_rank::Kingdom,    taxon_rank::Domain
    };

    if(stats.assigned() < 1) {
        os << "None of the input sequences could be classified.\n";
        return;
    }

    if(stats.unassigned() > 0) {
        os << prefix << "unclassified: "
           << (100 * stats.unclassified_rate())
           << "% (" << stats.unassigned() << ")\n";
    }
    os << prefix << "classified:\n";
    for(auto r : ranks) {
        if(stats.assigned(r) > 0) {
            std::string rn = taxonomy::rank_name(r);
            rn.resize(11, ' ');
            os  << prefix <<"  "<< rn
                << (100 * stats.classification_rate(r))
                << "% (" << stats.assigned(r) << ")\n";
        }
    }

    if(stats.known() > 0) {
        if(stats.unknown() > 0) {
            os << prefix << "ground truth unknown: "
               << (100 * stats.unknown_rate())
               << "% (" << stats.unknown() << ")\n";
        }
        os << prefix << "ground truth known:\n";
        for(auto r : ranks) {
            if(stats.assigned(r) > 0) {
                std::string rn = taxonomy::rank_name(r);
                rn.resize(11, ' ');
                os  << prefix <<"  "<< rn
                    << (100 * stats.known_rate(r))
                    << "% (" << stats.known(r) << ")\n";
            }
        }

        os << prefix << "correctly classified:\n";
        for(auto r : ranks) {
            if(stats.assigned(r) > 0) {
                std::string rn = taxonomy::rank_name(r);
                rn.resize(11, ' ');
                os << prefix <<"  "<< rn << stats.correct(r) << '\n';
            }
        }

        os << prefix << "precision (correctly classified / classified) if ground truth known:\n";
        for(auto r : ranks) {
            if(stats.assigned(r) > 0) {
                std::string rn = taxonomy::rank_name(r);
                rn.resize(11, ' ');
                os << prefix <<"  "<< rn << (100 * stats.precision(r)) << "%\n";
            }
        }

        os << prefix << "sensitivity (correctly classified / all) if ground truth known:\n";
        for(auto r : ranks) {
            if(stats.assigned(r) > 0) {
                std::string rn = taxonomy::rank_name(r);
                rn.resize(11, ' ');
                os << prefix <<"  "<< rn << (100 * stats.sensitivity(r)) << "%\n";
            }
        }

        if(stats.coverage(taxon_rank::Domain).total() > 0) {
            os << prefix << "false positives (hit on taxa not covered in DB):\n";
            for(auto r : ranks) {
                if(stats.assigned(r) > 0) {
                    std::string rn = taxonomy::rank_name(r);
                    rn.resize(11, ' ');
                    os << prefix <<"  "<< rn
                       << stats.coverage(r).false_pos() << "\n";
                }
            }
        }
    }

}


} // namespace mc

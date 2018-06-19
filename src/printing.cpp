/******************************************************************************
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

#include <set>
#include <ostream>
#include <utility>

#include "candidates.h"
#include "classification_statistics.h"
#include "matches_per_target.h"
#include "sketch_database.h"
#include "stat_confusion.h"
#include "taxonomy.h"

#include "printing.h"


namespace mc {

//-------------------------------------------------------------------
void show_query_parameters(std::ostream& os, const query_options& opt)
{
    const auto& comment = opt.output.format.comment;
    if(opt.output.mapViewMode != map_view_mode::none) {
        os << comment << "Reporting per-read mappings (non-mapping lines start with '"
           << comment << "').\n";

        if(opt.output.showLineage) {
            os << comment
               << "The complete lineage will be reported "
               << "starting with the lowest match.\n";
        }
        else {
            os << comment
               << "Only the lowest matching rank will be reported.\n";
        }
    }
    else {
        os << comment << "Per-Read mappings will not be shown.\n";
    }

    os << comment
       << "Classification will be constrained to ranks from '"
       << taxonomy::rank_name(opt.classify.lowestRank) << "' to '"
       << taxonomy::rank_name(opt.classify.highestRank) << "'.\n";

    os << comment
       << "Classification hit threshold is "
       << opt.classify.hitsMin << " per query\n";

    os << comment
       << "At maximum "
       << opt.classify.maxNumCandidatesPerQuery
       << " classification candidates will be considered per query.\n";

    if(opt.evaluate.excludeRank != taxon_rank::none) {
        os << comment << "Clade Exclusion on Rank: "
           << taxonomy::rank_name(opt.evaluate.excludeRank);
    }

    if(opt.process.pairing == pairing_mode::files) {
        os << comment << "File based paired-end mode:\n"
           << comment << "  Reads from two consecutive files will be interleaved.\n"
           << comment << "  Max insert size considered " << opt.classify.insertSizeMax << ".\n";
    }
    else if(opt.process.pairing == pairing_mode::sequences) {
        os << comment << "Per file paired-end mode:\n"
           << comment << "  Reads from two consecutive sequences in each file will be paired up.\n"
           << comment << "  Max insert size considered " << opt.classify.insertSizeMax << ".\n";
    }

    if(opt.output.showAlignment) {
        os << comment << "Query sequences will be aligned to best candidate target => SLOW!\n";
    }

    if(opt.output.showHitsPerTargetList) {
        os << comment << "A list of hits per reference sequence "
           << "will be generated after the read mapping.\n";
    }

    os << comment << "Using " << opt.process.numThreads << " threads\n";
}



//-------------------------------------------------------------------
void show_taxon(std::ostream& os,
                const taxon* tax, taxon_print_mode mode,
                const formatting_strings& fmt)
{
    if(!tax) return;
    switch(mode) {
        default:
        case taxon_print_mode::rank_name:
            os << taxonomy::rank_name(tax->rank()) << fmt.rankSuffix;
            // fall through
        case taxon_print_mode::name:
            os << tax->name();
            break;
        case taxon_print_mode::rank_id:
            os << taxonomy::rank_name(tax->rank()) << fmt.rankSuffix;
            // fall through
        case taxon_print_mode::id:
            os << tax->id();
            break;
        case taxon_print_mode::rank_name_id:
            os << taxonomy::rank_name(tax->rank()) << fmt.rankSuffix;
            // fall through
        case taxon_print_mode::name_id:
            os << tax->name() << fmt.taxidPrefix << tax->id() << fmt.taxidSuffix;
            break;
    }
}



//-------------------------------------------------------------------
void show_no_taxon(std::ostream& os,
                   taxon_print_mode mode,
                   taxon_rank rank,
                   const formatting_strings& fmt)
{
    switch(mode) {
        default:
        case taxon_print_mode::rank_name:
            os << taxonomy::rank_name(rank) << fmt.rankSuffix;
            // fall through
        case taxon_print_mode::name:
            os << fmt.none;
            break;
        case taxon_print_mode::rank_id:
            os << taxonomy::rank_name(rank) << fmt.rankSuffix;
            // fall through
        case taxon_print_mode::id:
            os << taxonomy::none_id();
            break;
        case taxon_print_mode::rank_name_id:
            os << taxonomy::rank_name(rank) << fmt.rankSuffix;
            // fall through
        case taxon_print_mode::name_id:
            os << fmt.none << fmt.taxidPrefix << taxonomy::none_id()
               << fmt.taxidSuffix;
            break;
    }
}



//-------------------------------------------------------------------
void show_lineage(std::ostream& os,
                  const ranked_lineage& lineage,
                  taxon_print_mode mode, taxon_rank lowest, taxon_rank highest,
                  const formatting_strings& fmt)
{
    if(lowest == taxon_rank::none) return;
    if(highest == taxon_rank::none) highest = taxon_rank::root;

    for(auto r = lowest; r <= highest; ++r) {
        const taxon* tax = lineage[int(r)];
        if(tax) {
            show_taxon(os, tax, mode, fmt);
        }
        else {
            show_no_taxon(os, mode, r, fmt);
        }
        if(r < highest) {
            os << fmt.taxSeparator;
        }
    }

}

//-------------------------------------------------------------------
void show_blank_lineage(std::ostream& os,
                        taxon_print_mode mode,
                        taxon_rank lowest, taxon_rank highest,
                        const formatting_strings& fmt)
{
    for(auto r = lowest; r <= highest; ++r) {
        switch(mode) {
            default:
            case taxon_print_mode::rank_name:
                os << fmt.none << fmt.rankSuffix;
                // fall through
            case taxon_print_mode::name:
                os << fmt.none;
                break;
            case taxon_print_mode::rank_id:
                os << fmt.none << fmt.rankSuffix;
                // fall through
            case taxon_print_mode::id:
                os << taxonomy::none_id();
                break;
            case taxon_print_mode::rank_name_id:
                os << fmt.none << fmt.rankSuffix;
                // fall through
            case taxon_print_mode::name_id:
                os << fmt.none
                   << fmt.taxidPrefix << taxonomy::none_id()
                   << fmt.taxidSuffix;
                break;
        }
        if(r < highest) os << fmt.taxSeparator;
    }
}



//-------------------------------------------------------------------
void show_taxon_header(std::ostream& os,
                       const classification_output_options& opt,
                       const std::string& prefix)
{
    const auto rmax = opt.showLineage ? opt.highestRank : opt.lowestRank;

    if(opt.lowestRank == rmax) {
        switch(opt.showTaxaAs) {
            default:
            case taxon_print_mode::rank_name:
                os << prefix << "rank" << opt.format.rankSuffix;
                // fall through
            case taxon_print_mode::name:
                os << prefix << "taxname";
                break;
            case taxon_print_mode::rank_id:
                os << prefix << "rank" << opt.format.rankSuffix;
                // fall through
            case taxon_print_mode::id:
                os << prefix << "taxid";
                break;
            case taxon_print_mode::rank_name_id:
                os << prefix << "rank" << opt.format.rankSuffix;
                // fall through
            case taxon_print_mode::name_id:
                os << prefix << "taxname" << opt.format.taxidPrefix
                   << prefix << "taxid" << opt.format.taxidSuffix;
                break;
        }
    }
    else {
        for(auto r = opt.lowestRank; r <= rmax; ++r) {
            switch(opt.showTaxaAs) {
                default:
                case taxon_print_mode::rank_name:
                    os << prefix << taxonomy::rank_name(r)
                       << opt.format.rankSuffix;
                    // fall through
                case taxon_print_mode::name:
                    os << prefix << "taxname";
                    break;
                case taxon_print_mode::rank_id:
                    os << prefix << taxonomy::rank_name(r)
                       << opt.format.rankSuffix;
                    // fall through
                case taxon_print_mode::id:
                    os << prefix << "taxid";
                    break;
                case taxon_print_mode::rank_name_id:
                    os << prefix << taxonomy::rank_name(r)
                       << opt.format.rankSuffix;
                    // fall through
                case taxon_print_mode::name_id:
                    os << prefix << "taxname" << opt.format.taxidPrefix
                       << prefix << "taxid" << opt.format.taxidSuffix;
                    break;
            }
            if(r < rmax) os << opt.format.taxSeparator;
        }
    }
}



//-------------------------------------------------------------------
void show_taxon(std::ostream& os,
                const database& db,
                const classification_output_options& opt,
                const taxon* tax)
{
    if(!tax || tax->rank() > opt.highestRank) {
        if(opt.collapseUnclassified) {
            if(opt.showTaxaAs == taxon_print_mode::id) {
                os << taxonomy::none_id();
            } else {
                os << opt.format.none;
            }
        } else {
            const auto rmax = opt.showLineage ? opt.highestRank : opt.lowestRank;
            show_blank_lineage(os, opt.showTaxaAs, opt.lowestRank, rmax, opt.format);
        }
    }
    else {
        const auto rmin = opt.lowestRank < tax->rank() ? tax->rank() : opt.lowestRank;
        const auto rmax = opt.showLineage ? opt.highestRank : rmin;

        show_lineage(os, db.ranks(tax), opt.showTaxaAs, rmin, rmax, opt.format);
    }
}



//-------------------------------------------------------------------
void show_matches(std::ostream& os,
                  const database& db,
                  const classification_candidates& cand,
                  taxon_rank lowest)
{
    using size_t = classification_candidates::size_type;

    if(lowest == taxon_rank::Sequence) {
        for(size_t i = 0; i < cand.size() && cand[i].hits > 0; ++i) {
            if(i > 0) os << ',';
            if(cand[i].tax) os << cand[i].tax->name() << ':' << cand[i].hits;
        }
    }
    else {
        for(size_t i = 0; i < cand.size() && cand[i].hits > 0; ++i) {
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

    for(const auto& c : cand) {
        os << '[' << (w * c.pos.beg)
           << ',' << (w * c.pos.end + db.target_window_size()) << "] ";
    }
}



//-------------------------------------------------------------------
void show_matches_per_targets(std::ostream& os,
                              const database& db,
                              const matches_per_target& tgtMatches,
                              const classification_output_options& opt)
{
    os << opt.format.comment << "TABLE_LAYOUT: "
        << " sequence " << opt.format.column
        << " windows_in_sequence " << opt.format.column
        << " queryid/window_index:hits/window_index:hits/...,queryid/...\n";

    for(const auto& mapping : tgtMatches) {
        show_taxon(os, db, opt, mapping.first);
        os << opt.format.column << mapping.first->source().windows
           << opt.format.column;

        bool first = true;
        for(const auto& candidate : mapping.second) {
            if(first) first = false; else os << ',';
            os << candidate.qeryid;
            // windows with hit count
            for(const auto& match : candidate.matches) {
                os << '/' << match.win << ':' << match.hits;
            }
        }
        os << '\n';
    }
}


//-------------------------------------------------------------------
void show_num_matches_per_targets(std::ostream& os,
                                  const database& db,
                                  const matches_per_target& tgtMatches,
                                  const classification_output_options& opt)
{
    os << opt.format.comment << "TABLE_LAYOUT: "
        << " sequence " << opt.format.column
        << " windows_in_sequence " << opt.format.column
        << " number of matched windows (including duplicates)\n";

    for(const auto& mapping : tgtMatches) {
        show_taxon(os, db, opt, mapping.first);
        os << opt.format.column << mapping.first->source().windows
           << opt.format.column << mapping.second.size()
           << '\n';
    }
}

//-------------------------------------------------------------------
void show_features_of_targets(std::ostream& os,
                                  const database& db,
                                  const matches_per_target& tgtMatches,
                                  const classification_output_options& opt)
{
    os << opt.format.comment << "TABLE_LAYOUT: "
    << " sequence " << opt.format.column
    << " windows_in_sequence " << opt.format.column
    << " coverage percentage " << opt.format.column
    << " longest gap\n";

    //calculate coverage percentages
    for(const auto& mapping : tgtMatches) {
        const taxon* target = mapping.first;
        std::set<window_id> hitWindows;
        for(const auto& candidate : mapping.second) {
            for(const auto& windowMatch : candidate.matches) {
                hitWindows.emplace(windowMatch.win);
            }
        }
        const window_id targetSize = target->source().windows;
        const float covP = float(hitWindows.size()) / targetSize;

        window_id longestGapSize = 0;
        window_id prevWin = 0;
        for(const auto& win : hitWindows) {
            const window_id gapSize = win - prevWin;
            if(gapSize > longestGapSize)
                longestGapSize = gapSize;
            prevWin = win;
        }
        const window_id gapSize = (targetSize-1) - prevWin;
        if(gapSize > longestGapSize)
            longestGapSize = gapSize;

        show_taxon(os, db, opt, mapping.first);
        os << opt.format.column << targetSize
           << opt.format.column << covP
           << opt.format.column << longestGapSize
           << '\n';
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

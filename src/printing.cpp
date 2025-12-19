/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2026 André Müller (github.com/muellan)
 *                       & Robin Kobus  (github.com/funatiq)
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


#include "printing.hpp"

#include "classification.hpp"
#include "classification_statistics.hpp"
#include "database.hpp"
#include "matches_per_target.hpp"
#include "options.hpp"
#include "stat_confusion.hpp"
#include "taxonomy.hpp"
#include "typename.hpp"
#include "version.hpp"

#include <ostream>
#include <utility>
#include <iomanip>
#include <cmath>


namespace mc {

using std::cout;
using std::endl;


//-----------------------------------------------------------------------------
void show_query_parameters (std::ostream& os, const query_options& opt)
{
    const auto& analysis = opt.output.analysis;
    const auto& fmt = opt.output.format;
    const auto& comment = fmt.tokens.comment;

    if (fmt.mapViewMode != map_view_mode::none) {
        os << comment << "Reporting per-read mappings (non-mapping lines start with '"
           << comment << "').\n";

        if (fmt.showLineage) {
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

    if (opt.minReadLength > 0) {
        os << comment << "Only reads with a minimum length of "
           << opt.minReadLength << " bp will be mapped.\n";
    }
    if (opt.maxReadLength < std::numeric_limits<std::size_t>::max()) {
        os << comment << "Only reads with a maximum length of "
           << opt.maxReadLength << " bp will be mapped.\n";
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

    if (opt.pairing == pairing_mode::files) {
        os << comment << "File based paired-end mode:\n"
           << comment << "  Reads from two consecutive files will be interleaved.\n"
           << comment << "  Max insert size considered " << opt.classify.insertSizeMax << ".\n";
    }
    else if (opt.pairing == pairing_mode::sequences) {
        os << comment << "Per file paired-end mode:\n"
           << comment << "  Reads from two consecutive sequences in each file will be paired up.\n"
           << comment << "  Max insert size considered " << opt.classify.insertSizeMax << ".\n";
    }

    if (analysis.showAlignment) {
        os << comment << "Query sequences will be aligned to best candidate target => SLOW!\n";
    }

    if (analysis.showHitsPerTargetList) {
        os << comment << "A list of hits per reference sequence "
           << "will be generated after the read mapping.\n";
    }

    if (analysis.showTaxAbundances) {
        os << comment << "A list of absolute and relative abundances per taxon "
           << "will be generated after the read mapping.\n";
    }

    if (analysis.showAbundanceEstimatesOnRank != taxon_rank::none) {
        os << comment << "A list of absolute and relative abundances for each '"
           << taxonomy::rank_name(analysis.showAbundanceEstimatesOnRank)
           << "' will be generated after the read mapping.\n";
    }

    os << comment << "Using " << opt.performance.numThreads << " threads\n";
}



//-------------------------------------------------------------------
void show_taxon_header (std::ostream& os,
                        const classification_output_formatting& opt,
                        const std::string& prefix)
{
    const auto rmax = opt.showLineage ? opt.highestRank : opt.lowestRank;
    const auto& style = opt.taxonStyle;
    const auto& fmt = opt.tokens;

    if (opt.lowestRank == rmax) {
        if (style.showRankName) {
            os << prefix << "rank" << fmt.rankSuffix;
        }

        if (style.showName) {
            os << prefix << "taxname";
            if (style.showId) {
                os << fmt.taxidPrefix << prefix << "taxid" << fmt.taxidSuffix;
            }
        }
        else if (style.showId) {
            os << prefix << "taxid";
        }
    }
    else {
        for (auto r = opt.lowestRank; r <= rmax; ++r) {

            if (style.showRankName) {
                os << prefix << taxonomy::rank_name(r) << fmt.rankSuffix;
            }

            if (style.showName) {
                os << prefix << "taxname";
                if (style.showId) {
                    os << fmt.taxidPrefix << prefix << "taxid" << fmt.taxidSuffix;
                }
            }
            else if (style.showId) {
                os << prefix << "taxid";
            }

            if (r < rmax) os << opt.tokens.taxSeparator;
        }
    }
}



//-------------------------------------------------------------------
void print_taxon (std::ostream& os,
                  const std::string& taxName,
                  taxon_id id,
                  taxon_rank rank,
                  taxon_print_style style,
                  const formatting_tokens& fmt)
{
    if (style.showRankName) {
        if (rank == taxon_rank::none) {
            os << fmt.none << fmt.rankSuffix;
        } else {
            os << taxonomy::rank_name(rank) << fmt.rankSuffix;
        }
    }

    if (style.showName) {
        os << taxName;
        if (style.showId) {
            os << fmt.taxidPrefix << id << fmt.taxidSuffix;
        }
    }
    else if (style.showId) {
        os << id;
    }
}



//-------------------------------------------------------------------
void show_lineage (std::ostream& os,
                   const ranked_lineage& lineage,
                   taxon_print_style style, taxon_rank lowest, taxon_rank highest,
                   const formatting_tokens& fmt)
{
    if (lowest == taxon_rank::none) return;
    if (highest == taxon_rank::none) highest = taxon_rank::root;

    for (auto r = lowest; r <= highest; ++r) {
        const taxon* tax = lineage[int(r)];
        if (tax) {
            print_taxon(os, tax->name(), tax->id(), tax->rank(), style, fmt);
        }
        else {
            print_taxon(os, fmt.none, taxonomy::none_id(), r, style, fmt);
        }
        if (r < highest) {
            os << fmt.taxSeparator;
        }
    }

}



//-------------------------------------------------------------------
void show_blank_lineage (std::ostream& os,
                         taxon_print_style style,
                         taxon_rank lowest, taxon_rank highest,
                         const formatting_tokens& fmt)
{
    for (auto r = lowest; r <= highest; ++r) {
        print_taxon(os, fmt.none, taxonomy::none_id(), taxon_rank::none, style, fmt);
        if (r < highest) os << fmt.taxSeparator;
    }
}



//-------------------------------------------------------------------
void show_taxon (std::ostream& os,
                 const taxonomy_cache& taxonomy,
                 const classification_output_formatting& opt,
                 const taxon* tax)
{
    if (not tax || tax->rank() > opt.highestRank) {
        if (opt.collapseUnclassifiedLineages) {
            if (opt.taxonStyle.showId
               && not opt.taxonStyle.showName
               && not opt.taxonStyle.showRankName)
            {
                os << taxonomy::none_id();
            }
            else {
                os << opt.tokens.none;
            }
        }
        else {
            const auto rmax = opt.showLineage ? opt.highestRank : opt.lowestRank;
            show_blank_lineage(os, opt.taxonStyle, opt.lowestRank, rmax, opt.tokens);
        }
    }
    else {
        const auto rmin = opt.lowestRank < tax->rank() ? tax->rank() : opt.lowestRank;
        const auto rmax = opt.showLineage ? opt.highestRank : rmin;

        show_lineage(os, taxonomy.cached_ranks(tax), opt.taxonStyle, rmin, rmax, opt.tokens);
    }
}



//-------------------------------------------------------------------
void show_candidates (std::ostream& os,
                      const taxonomy_cache& taxonomy,
                      const span<const match_candidate> cand,
                      taxon_rank lowest)
{
    using size_t = span<match_candidate>::size_type;

    if (lowest == taxon_rank::Sequence) {
        for (size_t i = 0; i < cand.size() && cand[i].hits > 0; ++i) {
            if (i > 0) os << ',';
            if (cand[i].tax) os << cand[i].tax->name() << ':' << cand[i].hits;
        }
    }
    else {
        for (size_t i = 0; i < cand.size() && cand[i].hits > 0; ++i) {
            if (i > 0) os << ',';
            const taxon* tax = (cand[i].tax->rank() < lowest) ?
                               taxonomy.cached_ancestor(cand[i].tgt,lowest) :
                               cand[i].tax;
            if (tax) {
                os << tax->id();
            } else {
                os << cand[i].tax->name();
            }
            os << ':' << cand[i].hits;
        }
    }
}



//-------------------------------------------------------------------
void show_matches (std::ostream& os,
                   const taxonomy_cache& taxonomy,
                   const span<const location> matches,
                   taxon_rank lowest)
{
    if (matches.empty()) return;

    if (lowest == taxon_rank::Sequence) {
        auto cur = matches.begin();
        int count = 1;
        for (auto it = matches.begin()+1; it != matches.end(); ++it) {
            if (*cur == *it)
                ++count;
            else {
                const taxon* tax = taxonomy.cached_taxon_of_target(cur->tgt);
                if (tax) os << tax->name()
                           << '/' << int(cur->win)
                           << ':' << count << ',';
                cur = it;
                count = 1;
            }
        }
        const taxon* tax = taxonomy.cached_taxon_of_target(cur->tgt);
        if (tax) os << tax->name()
                   << '/' << int(cur->win)
                   << ':' << count << ',';
    }
    else {
        auto cur = matches.begin();
        int count = 1;
        for (auto it = matches.begin()+1; it != matches.end(); ++it) {
            if (*cur == *it)
                ++count;
            else {
                const taxon* tax = taxonomy.cached_ancestor(cur->tgt, lowest);
                if (not tax) {
                    tax = taxonomy.cached_taxon_of_target(cur->tgt);
                }
                os << tax->name() << ':' << count << ',';

                cur = it;
                count = 1;
            }
        }
        const taxon* tax = taxonomy.cached_ancestor(cur->tgt, lowest);
        if (not tax) {
            tax = taxonomy.cached_taxon_of_target(cur->tgt);
        }
        os << tax->name() << ':' << count << ',';
    }
}



//-------------------------------------------------------------------
void show_candidate_ranges (std::ostream& os,
                            const sketching_opt& targetSketching,
                            const span<const match_candidate> cand)
{
    const auto w = targetSketching.winstride;

    for (const auto& c : cand) {
        os << '[' << (w * c.pos.beg)
           << ',' << (w * c.pos.end + targetSketching.winlen) << "] ";
    }
}



//-------------------------------------------------------------------
void show_matches_per_targets (std::ostream& os,
                               const database& db,
                               const matches_per_target& tgtMatches,
                               const classification_output_formatting& opt)
{
    os << opt.tokens.comment
       << "--- list of hits for each reference sequence ---\n"
       << opt.tokens.comment
       << "window start position within sequence = window_index * window_stride(="
       << db.target_sketching().winstride << ")\n";

    os << opt.tokens.comment << "TABLE_LAYOUT: "
       << " sequence " << opt.tokens.column
       << " windows_in_sequence " << opt.tokens.column
       << "queryid/first_window_index+additional_windows:hits,queryid/...\n";

    const auto rmin = taxon_rank::Sequence;
    const auto rmax = opt.showLineage ? opt.highestRank : rmin;

    for (const auto& mapping : tgtMatches) {
        show_lineage(os, db.taxa().cached_ranks(mapping.first), opt.taxonStyle, rmin, rmax, opt.tokens);
        os << opt.tokens.column << db.taxa().cached_taxon_of_target(mapping.first)->source().windows
           << opt.tokens.column;

        bool first = true;
        for (const auto& candidate : mapping.second) {
            if (first) first = false; else os << ',';
            os << candidate.qid << '/'
               << candidate.pos.beg << '+'
               << (candidate.pos.end-candidate.pos.beg) << ':'
               << candidate.hits;
        }
        os << '\n';
    }
}



//-------------------------------------------------------------------
void show_abundance_table (std::ostream& os,
                           const taxon_count_map& allTaxCounts,
                           const classification_statistics& statistics,
                           const classification_output_formatting&  opt)
{
    os << opt.tokens.comment
       << "rank" << opt.tokens.rankSuffix
       << "name" << opt.tokens.column
       << "taxid" << opt.tokens.column
       << "number of reads" << opt.tokens.column
       << "abundance\n";

    double ipart = 0.0;

    for (const auto& tc : allTaxCounts) {
        if (tc.first) {
            os << tc.first->rank_name() << opt.tokens.rankSuffix
               << tc.first->name() << opt.tokens.column;
            if (tc.first->rank() == taxon_rank::Sequence)
               os << tc.first->parent_id();
            else
               os << tc.first->id();
        } else {
            os << "none";
        }

        os << opt.tokens.column;
        // ensure that integer values are printed as such
        // and not in scientific notation
        if (std::modf(tc.second, &ipart) == 0.0) {
            os << ipart;
        } else {
            // set higher precision to prevent scientific notation
            // and set back to default (6) after that
            os << std::setprecision(15) << tc.second << std::setprecision(6);
        }
        os << opt.tokens.column
           << (tc.second / double(statistics.total()) * 100) << "%\n";
    }
    os << "unclassified" << opt.tokens.column
       << "--" << opt.tokens.column
       << '0' << opt.tokens.column
       << statistics.unassigned() << opt.tokens.column
       << statistics.unclassified_rate() * 100 << "%\n";
}



//-------------------------------------------------------------------
void show_abundances (std::ostream& os,
                      const taxon_count_map& allTaxCounts,
                      const classification_statistics& statistics,
                      const classification_output_formatting& opt)
{
    os << opt.tokens.comment
       << "query summary: number of queries mapped per taxon\n";
    show_abundance_table(os, allTaxCounts, statistics, opt);
}



//-------------------------------------------------------------------
void show_abundance_estimates (std::ostream& os,
                               taxon_rank onRank,
                               const taxon_count_map& allTaxCounts,
                               const classification_statistics& statistics,
                               const classification_output_formatting& opt)
{
    os << opt.tokens.comment
       << "estimated abundance (number of queries) per "
       << taxonomy::rank_name(onRank) << "\n";

    show_abundance_table(os, allTaxCounts, statistics, opt);
}



//-------------------------------------------------------------------
void show_taxon_statistics (std::ostream& os,
                            const classification_statistics& stats,
                            const std::string& prefix)
{
    constexpr taxon_rank ranks[] {
          taxon_rank::Sequence,
          taxon_rank::subSpecies,
          taxon_rank::Species,    taxon_rank::Genus,
          taxon_rank::Family,     taxon_rank::Order,
          taxon_rank::Class,      taxon_rank::Phylum,
          taxon_rank::Kingdom,    taxon_rank::Domain,
          taxon_rank::root
    };

    if (stats.assigned() < 1) {
        os << "None of the input sequences could be classified.\n";
        return;
    }

    if (stats.unassigned() > 0) {
        os << prefix << "unclassified: "
           << (100 * stats.unclassified_rate())
           << "% (" << stats.unassigned() << ")\n";
    }
    os << prefix << "classified:\n";
    for (auto r : ranks) {
        if (stats.assigned(r) > 0) {
            std::string rn = taxonomy::rank_name(r);
            rn.resize(11, ' ');
            os  << prefix <<"  "<< rn
                << (100 * stats.classification_rate(r))
                << "% (" << stats.assigned(r) << ")\n";
        }
    }

    if (stats.known() > 0) {
        if (stats.unknown() > 0) {
            os << prefix << "ground truth unknown: "
               << (100 * stats.unknown_rate())
               << "% (" << stats.unknown() << ")\n";
        }
        os << prefix << "ground truth known:\n";
        for (auto r : ranks) {
            if (stats.assigned(r) > 0) {
                std::string rn = taxonomy::rank_name(r);
                rn.resize(11, ' ');
                os  << prefix <<"  "<< rn
                    << (100 * stats.known_rate(r))
                    << "% (" << stats.known(r) << ")\n";
            }
        }

        os << prefix << "correctly classified:\n";
        for (auto r : ranks) {
            if (stats.assigned(r) > 0) {
                std::string rn = taxonomy::rank_name(r);
                rn.resize(11, ' ');
                os << prefix <<"  "<< rn << stats.correct(r) << '\n';
            }
        }

        os << prefix << "precision (correctly classified / classified) if ground truth known:\n";
        for (auto r : ranks) {
            if (stats.assigned(r) > 0) {
                std::string rn = taxonomy::rank_name(r);
                rn.resize(11, ' ');
                os << prefix <<"  "<< rn << (100 * stats.precision(r)) << "%\n";
            }
        }

        os << prefix << "sensitivity (correctly classified / all) if ground truth known:\n";
        for (auto r : ranks) {
            if (stats.assigned(r) > 0) {
                std::string rn = taxonomy::rank_name(r);
                rn.resize(11, ' ');
                os << prefix <<"  "<< rn << (100 * stats.sensitivity(r)) << "%\n";
            }
        }

        if (stats.coverage(taxon_rank::Domain).total() > 0) {
            os << prefix << "false positives (hit on taxa not covered in DB):\n";
            for (auto r : ranks) {
                if (stats.assigned(r) > 0) {
                    std::string rn = taxonomy::rank_name(r);
                    rn.resize(11, ' ');
                    os << prefix <<"  "<< rn
                       << stats.coverage(r).false_pos() << "\n";
                }
            }
        }
    }

}


//-----------------------------------------------------------------------------
/**
 * @brief show summary and statistics of classification
 */
void show_summary (const query_options& opt,
                   const classification_results& results)
{
    const auto& statistics = results.statistics;
    const auto numQueries = (opt.pairing == pairing_mode::none)
                            ? statistics.total() : 2 * statistics.total();

    const auto speed = numQueries / results.time.minutes();
    const auto& comment = opt.output.format.tokens.comment;
    results.perReadOut
        << comment << "queries: " << numQueries << '\n'
        << comment << "time:    " << results.time.milliseconds() << " ms\n"
        << comment << "speed:   " << speed << " queries/min\n";

    if (statistics.total() > 0) {
        show_taxon_statistics(results.perReadOut, statistics, comment);
    } else {
        results.status << comment << "No valid query sequences found.\n";
    }
}



//-------------------------------------------------------------------
void print_static_properties (const database& db)
{
    using target_id = database::target_id;
    using window_id = database::window_id;
    using feature_t = database::feature;
    using bkt_sz_t  = database::bucket_size_type;

    cout<< "------------------------------------------------\n"
        << "MetaCache version  " << MC_VERSION_STRING << " (" << MC_VERSION << ")\n"
        << "database version   " << MC_DB_VERSION << '\n'
        << "------------------------------------------------\n"
        << "sequence type      " << type_name<database::sequence>() << '\n'
        << "target id type     " << type_name<target_id>() << " " << (sizeof(target_id)*CHAR_BIT) << " bits\n"
        << "target limit       " << std::uint64_t(db.max_target_count()) << '\n'
        << "------------------------------------------------\n"
        << "window id type     " << type_name<window_id>() << " " << (sizeof(window_id)*CHAR_BIT) << " bits\n"
        << "window limit       " << std::uint64_t(db.max_windows_per_target()) << '\n'
        << "window length      " << db.target_sketching().winlen << '\n'
        << "window stride      " << db.target_sketching().winstride << '\n'
        << "------------------------------------------------\n"
        << "sketcher type      " << type_name<database::sketcher>() << '\n'
        << "feature type       " << type_name<feature_t>() << " " << (sizeof(feature_t)*CHAR_BIT) << " bits\n"
        << "feature hash       " << type_name<database::feature_hash>() << '\n'
        << "kmer size          " << std::uint64_t(db.target_sketching().kmerlen) << '\n'
        << "kmer limit         " << std::uint64_t(db.target_sketching().max_kmer_size()) << '\n'
        << "sketch size        " << db.target_sketching().sketchlen << '\n'
        << "------------------------------------------------\n"
        << "bucket size type   " << type_name<bkt_sz_t>() << " " << (sizeof(bkt_sz_t)*CHAR_BIT) << " bits\n"
        << "max. locations     " << std::uint64_t(db.max_locations_per_feature()) << '\n'
        << "location limit     " << std::uint64_t(db.max_supported_locations_per_feature()) << '\n'
        << "------------------------------------------------"
        << endl;
}



//-----------------------------------------------------------------------------
void print_content_properties (const database& db)
{
    cout    << "------------------------------------------------\n"
            << "database parts     " << db.part_count() << '\n';

    if (db.target_count() > 0) {
        std::uint64_t numRankedTargets = 0;
        for (const auto& t : db.taxa().target_taxa()) {
            if (t.has_parent()) ++numRankedTargets;
        }

        cout<< "targets            " << db.target_count() << '\n'
            << "ranked targets     " << numRankedTargets << '\n'
            << "taxa in tree       " << db.taxa().non_target_taxon_count() << '\n';
    }

    if (db.feature_count() > 0) {
        auto lss = db.location_list_size_statistics();

        if (db.part_count() > 1) {
            cout
            << "------------------------------------------------\n"
            << "complete database (all parts):\n";
        }
        cout<< "buckets            " << db.bucket_count() << '\n'
            << "bucket size        " << "max: " << lss.max()
                                     << " mean: " << lss.mean()
                                     << " +/- " << lss.stddev()
                                     << " <> " << lss.skewness() << '\n'
            << "features           " << std::uint64_t(lss.size()) << '\n'
            << "dead features      " << db.dead_feature_count() << '\n'
            << "locations          " << std::uint64_t(lss.sum()) << '\n';
    }
    cout    << "------------------------------------------------\n";
}


} // namespace mc

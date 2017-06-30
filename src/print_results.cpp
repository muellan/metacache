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

#include "print_results.h"


namespace mc {


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
           << ',' << (w * cand[i].pos.end) << "] ";
    }
}



//-------------------------------------------------------------------
void show_classification_statistics(std::ostream& os,
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

    auto alscore = stats.alignment_scores();
    if(!alscore.empty()) {
        os << prefix << "semi-global alignment of query to best candidates:\n"
           << prefix << "  score: " << alscore.mean()
                     << " +/- " << alscore.stddev()
                     << " <> " << alscore.skewness() << '\n';
    }

}


} // namespace mc

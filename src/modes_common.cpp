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

#include "modes_common.h"
#include "sequence_io.h"


namespace mc {


//-------------------------------------------------------------------
classification
ground_truth(const database& db, const std::string& header)
{

    //try to extract query id and find the corresponding target in database
    auto tid = db.target_id_of_sequence(extract_ncbi_accession_version_number(header));
    //not found yet
    if(!db.is_valid(tid)) {
        tid = db.target_id_of_sequence(extract_ncbi_accession_number(header));
        //not found yet
//        if(!db.is_valid(tid)) {
//            tid = db.target_id_of_sequence(extract_genbank_identifier(header));
//        }
    }

    auto taxid = extract_taxon_id(header);
    //make sure to get the lowest taxon that has a rank assigned
    if(taxid > 0) {
        for(auto id : db.ranks(db.taxon_with_id(taxid))) {
            if(id > 0) {
                taxid = id;
                break;
            }
        }
    }

    //if target known and in db
    if(db.is_valid(tid)) {
        //if taxid could not be determined solely from header, use the one
        //from the target in the db
        if(taxid < 1) {
            //use lowest taxon that has a rank assigned
            for(auto id : db.ranks_of_target(tid)) {
                if(id > 0) {
                    taxid = id;
                    break;
                }
            }
        }

        return classification { tid, taxid > 0 ? &db.taxon_with_id(taxid) : nullptr };
    }

    return classification { taxid > 0 ? &db.taxon_with_id(taxid) : nullptr };
}



//-------------------------------------------------------------------
taxon_rank
lowest_common_rank(const database& db,
                   const classification& a,
                   const classification& b)
{
    if(a.sequence_level() && b.sequence_level() &&
        a.target() == b.target()) return taxon_rank::Sequence;

    if(a.has_taxon() && b.has_taxon()) {
        return db.ranked_lca(a.tax(), b.tax()).rank;
    }

    return taxon_rank::none;
}



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
                os << db.taxon_with_id(taxid).name;
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
                if(taxon.rank >= lowest && taxon.rank <= highest) {
                    os << taxon.rank_name() << ':';
                    if(mode != taxon_print_mode::id_only) {
                        os << taxon.name;
                        if(mode != taxon_print_mode::name_only) {
                            os << "(" << taxon.id << ")";
                        }
                    }
                    else {
                        os << taxon.id;
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



//-------------------------------------------------------------------
void show_matches(std::ostream& os,
                  const database& db,
                  const match_result& matches,
                  taxon_rank lowest)
{
    if(matches.empty()) return;

    if(lowest == taxon_rank::Sequence) {
        for(const auto& r : matches) {
            os << db.sequence_id_of_target(r.first.tgt)
               << '/' << int(r.first.win)
               << ':' << int(r.second) << ',';
        }
    }
    else {
        for(const auto& r : matches) {
            auto taxid = db.ranks_of_target(r.first.tgt)[int(lowest)];
            os << taxid << ':' << int(r.second) << ',';
        }
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
            auto rn = taxonomy::rank_name(r);
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
            auto rn = taxonomy::rank_name(r);
            rn.resize(11, ' ');
            os  << prefix <<"  "<< rn
                << (100 * stats.known_rate(r))
                << "% (" << stats.known(r) << ")\n";
        }

        os << prefix << "correctly classified:\n";
        for(auto r : ranks) {
            auto rn = taxonomy::rank_name(r);
            rn.resize(11, ' ');
            os << prefix <<"  "<< rn << stats.correct(r) << '\n';
        }

        os << prefix << "precision (correctly classified / classified) if ground truth known:\n";
        for(auto r : ranks) {
            auto rn = taxonomy::rank_name(r);
            rn.resize(11, ' ');
            os << prefix <<"  "<< rn << (100 * stats.precision(r)) << "%\n";
        }

        os << prefix << "sensitivity (correctly classified / all) if ground truth known:\n";
        for(auto r : ranks) {
            if(stats.assigned(r) > 0) {
                auto rn = taxonomy::rank_name(r);
                rn.resize(11, ' ');
                os << prefix <<"  "<< rn << (100 * stats.sensitivity(r)) << "%\n";
            }
        }

        if(stats.coverage(taxon_rank::Domain).total() > 0) {
            os << prefix << "false positives (hit on taxa not covered in DB):\n";
            for(auto r : ranks) {
                if(stats.assigned(r) > 0) {
                    auto rn = taxonomy::rank_name(r);
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



//-------------------------------------------------------------------
void update_coverage_statistics(
    const database& db,
    const classification& result, const classification& truth,
    classification_statistics& stats)
{
    const auto& lin = db.ranks(truth.tax());
    //check if taxa are covered in DB
    for(auto taxid : lin) {
        auto r = db.taxon_with_id(taxid).rank;
        if(db.covers_taxon(taxid)) {
            if(r < result.rank()) { //unclassified on rank
                stats.count_coverage_false_neg(r);
            } else { //classified on rank
                stats.count_coverage_true_pos(r);
            }
        }
        else {
            if(r < result.rank()) { //unclassified on rank
                stats.count_coverage_true_neg(r);
            } else { //classified on rank
                stats.count_coverage_false_pos(r);
            }
        }
    }
}


} // namespace mc

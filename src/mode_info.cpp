/*****************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * version 1.0
 *
 * Copyright (C) 2016 André Müller (muellan@uni-mainz.de)
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
#include <map>

#include "args_handling.h"
#include "sequence_io.h"

#include "modes_common.h"


namespace mc {

/*****************************************************************************
 *
 *
 *****************************************************************************/
void show_database_statistics(const args_parser& args)
{
    auto dbfilename = database_name(args);
    auto db = make_database<database>(dbfilename);
    print_statistics(db);
}



/*****************************************************************************
 *
 *
 *****************************************************************************/
void show_feature_map(const args_parser& args)
{
    auto dbfilename = database_name(args);
    auto db = make_database<database>(dbfilename);
    print_statistics(db);
    std::cout << "===================================================\n";
    db.print_feature_map(std::cout);
    std::cout << "===================================================\n";
}



/*****************************************************************************
 *
 *
 *****************************************************************************/
void show_ranks_of_genome(const database& db, database::genome_id gid)
{
    //if genomes don't have their own taxonId, print their sequence id
    std::cout << "    sequence:   " << db.sequence_id_of_genome(gid);

    for(auto taxid : db.ranks_of_genome(gid)) {
        if(taxid > 1) {
            auto&& taxon = db.taxon_with_id(taxid);
            auto rn = taxon.rank_name() + ":";
            rn.resize(12, ' ');
            std::cout << "\n    " << rn << "(" << taxon.id << ") " << taxon.name;
        }
    }

    std::cout << '\n';
}



/*****************************************************************************
 *
 * @brief
 *
 *****************************************************************************/
void show_sequence_info(const database& db, database::genome_id gid)
{
    std::cout
        << "Reference sequence " << gid << " ("
        << db.sequence_id_of_genome(gid) << "):\n"
        << "    origin:     " << db.origin_of_genome(gid).filename << " / "
        << db.origin_of_genome(gid).index << '\n';

    show_ranks_of_genome(db, gid);
}



/*****************************************************************************
 *
 * @brief
 *
 *****************************************************************************/
void show_sequence_info(const args_parser& args)
{
    auto dbfilename = database_name(args);

    auto sids = std::vector<std::string>{};
    for(int i = 2; i < args.non_prefixed_count(); ++i) {
        sids.push_back(args.non_prefixed(i));
    }

    if(sids.empty()) return;

    auto db = make_database_metadata_only<database>(dbfilename);

    for(const auto& sid : sids) {
        auto gid = db.genome_id_of_sequence(sid);
        if(gid < db.genome_count()) {
            show_sequence_info(db, gid);
        }
        else {
            std::cout << "Reference sequence " << sid
                      << " not found in database.\n";
        }
    }
}



/*****************************************************************************
 *
 * @brief
 *
 *****************************************************************************/
void show_lineage_table(const args_parser& args)
{
    using rank = taxonomy::taxon_rank;

    auto dbfilename = database_name(args);

    auto db = make_database_metadata_only<database>(dbfilename);
    if(db.genome_count() < 1) return;

    //table header
    std::cout << taxonomy::rank_name(rank::Sequence);
    for(auto r = rank::Sequence; r <= rank::Domain; ++r) {
        std::cout << '\t' << taxonomy::rank_name(r);
    }
    std::cout << '\n';

    //rows
    for(genome_id gid = 0; gid < db.genome_count(); ++gid) {
        std::cout << db.sequence_id_of_genome(gid);
        auto ranks = db.ranks_of_genome(gid);
        for(auto r = rank::Sequence; r <= rank::Domain; ++r) {
            std::cout << '\t' << ranks[int(r)];
        }
        std::cout << '\n';
    }
}



/*****************************************************************************
 *
 * @brief
 *
 *****************************************************************************/
void show_all_meta_info(const args_parser& args)
{
    auto dbfilename = database_name(args);

    auto db = make_database_metadata_only<database>(dbfilename);
    if(db.genome_count() < 1) return;

    std::cout << "Properties of database " << dbfilename << ":\n";
    print_statistics(db);

    std::cout << "Reference sequences in database:\n";
    for(genome_id gid = 0; gid < db.genome_count(); ++gid) {
        show_sequence_info(db, gid);
    }
}



/*****************************************************************************
 *
 * @brief
 *
 *****************************************************************************/
void show_rank_statistics(const args_parser& args)
{
    auto rankName = args.non_prefixed(3);
    auto rank = taxonomy::rank_from_name(rankName);
    if(rank == taxonomy::taxon_rank::none) {
        std::cout << "rank not recognized" << std::endl;
        return;
    }

    auto dbfilename = database_name(args);

    auto db = make_database_metadata_only<database>(dbfilename);

    std::map<taxonomy::taxon_id, std::size_t> stat;

    for(genome_id i = 0; i < db.genome_count(); ++i) {
        auto tax = db.ranks_of_genome(i)[int(rank)];
        auto it = stat.find(tax);
        if(it != stat.end()) {
            ++(it->second);
        } else {
            stat.insert({tax, 1});
        }
    }

    std::cout << "Sequence distribution for rank " << rankName << ":" << std::endl;
    for(const auto& s : stat) {
        std::cout << db.taxon_with_id(s.first).name << " \t " << s.second << '\n';
    }
}



/*****************************************************************************
 *
 * @brief shows database properties
 *
 *****************************************************************************/
void main_mode_info(const args_parser& args)
{
    if(args.non_prefixed_count() > 2) {
        if(args.non_prefixed(2) == "lineages") {
            show_lineage_table(args);
        }
        else if(args.non_prefixed(2) == "rank" && args.non_prefixed_count() > 3) {
            show_rank_statistics(args);
        }
        else if(args.non_prefixed(2) == "statistics") {
            show_database_statistics(args);
        }
        else if(args.non_prefixed(2) == "featuremap") {
            show_feature_map(args);
        }
        else {
            show_sequence_info(args);
        }
    }
    else {
        show_all_meta_info(args);
    }
}


} // namespace mc

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
#include <map>

#include "args_handling.h"
#include "sequence_io.h"
#include "print_info.h"


namespace mc {

/*****************************************************************************
 *
 *
 *****************************************************************************/
void show_database_config(const args_parser& args)
{
    auto dbfilename = database_name(args);
    auto db = make_database<database>(dbfilename, database::scope::metadata_only);
    print_properties(db);
    std::cout << std::endl;
}



/*****************************************************************************
 *
 *
 *****************************************************************************/
void show_database_statistics(const args_parser& args)
{
    auto dbfilename = database_name(args);
    auto db = make_database<database>(dbfilename);
    print_properties(db);
    std::cout << std::endl;
}



/*****************************************************************************
 *
 *
 *****************************************************************************/
void show_feature_map(const args_parser& args)
{
    auto dbfilename = database_name(args);
    auto db = make_database<database>(dbfilename);
    print_properties(db);
    std::cout << "\n===================================================\n";
    db.print_feature_map(std::cout);
    std::cout << "===================================================\n";
}




/*****************************************************************************
 *
 * @brief
 *
 *****************************************************************************/
void show_info(const args_parser& args)
{
    auto dbfilename = database_name(args);

    auto sids = std::vector<std::string>{};
    for(std::size_t i = 3; i < args.non_prefixed_count(); ++i) {
        sids.push_back(args.non_prefixed(i));
    }

    auto db = make_database_metadata_only<database>(dbfilename);

    if(!sids.empty()) {
        for(const auto& sid : sids) {
            const auto& tax = db.taxon_of_sequence(sid);
            if(!tax.is_none()) {
                show_info(std::cout, db, tax);
            }
            else {
                std::cout << "Target (reference sequence) " << sid
                          << " not found in database.\n";
            }
        }
    }
    else {
        std::cout << "Targets (reference sequences) in database:\n";
        for(const auto& tax : db.target_taxa()) {
            show_info(std::cout, db, tax);
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
    using rank = taxonomy::rank;

    auto dbfilename = database_name(args);

    auto db = make_database_metadata_only<database>(dbfilename);
    if(db.target_count() < 1) return;

    //table header
    std::cout << taxonomy::rank_name(rank::Sequence);
    for(auto r = rank::Sequence; r <= rank::Domain; ++r) {
        std::cout << '\t' << taxonomy::rank_name(r);
    }
    std::cout << '\n';

    //rows
    for(const auto& tax : db.taxa()) {
        std::cout << tax.name();
        auto ranks = db.ranks(tax);
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
void show_rank_statistics(const args_parser& args)
{
    auto rankName = args.non_prefixed(3);
    auto rank = taxonomy::rank_from_name(rankName);
    if(rank == taxonomy::rank::none) {
        std::cerr << "rank not recognized" << std::endl;
        return;
    }

    auto dbfilename = database_name(args);

    auto db = make_database_metadata_only<database>(dbfilename);

    std::map<taxonomy::taxon_id, std::size_t> stat;

    for(const auto& tax : db.taxa()) {
        auto tid = db.ranks(tax)[int(rank)];
        auto it = stat.find(tid);
        if(it != stat.end()) {
            ++(it->second);
        } else {
            stat.insert({tid, 1});
        }
    }

    std::cout << "Sequence distribution for rank " << rankName << ":" << std::endl;
    for(const auto& s : stat) {
        std::cout << db.taxon_with_id(s.first).name() << " \t " << s.second << '\n';
    }
}



/*****************************************************************************
 *
 * @brief
 *
 *****************************************************************************/
void show_basic_exec_info()
{
    database db;
    print_properties(db);
    std::cout << std::endl;
}



/*****************************************************************************
 *
 * @brief shows database properties
 *
 *****************************************************************************/
void main_mode_info(const args_parser& args)
{
    const auto nargs = args.non_prefixed_count();

    if(nargs < 2) {
        show_basic_exec_info();
    }
    else if(nargs == 2) {
        try {
            show_database_config(args);
        }
        catch(std::exception& e) {
            std::cout << '\n' << e.what() << "!\n";
        }
    }
    else {
        const auto mode = args.non_prefixed(2);

        if(mode == "target" || mode == "targets"  ||
                mode == "sequ"   || mode == "sequence" || mode == "sequences")
        {
            show_info(args);
        }
        else if(mode == "lineages" || mode == "lineage" || mode == "lin")
        {
            show_lineage_table(args);
        }
        else if(mode == "rank" && args.non_prefixed_count() > 3) {
            show_rank_statistics(args);
        }
        else if(mode == "statistics" || mode == "stat") {
            show_database_statistics(args);
        }
        else if(mode == "features" || mode == "featuremap") {
            show_feature_map(args);
        }
        else {
            std::cerr << "Info mode '" << mode << "' not recognized.\n"
                      << "Properties of database " << args.non_prefixed(1)
                      << ":" << std::endl;
            show_database_config(args);
        }
    }
}


} // namespace mc

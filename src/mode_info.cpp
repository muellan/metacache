/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2018 André Müller (muellan@uni-mainz.de)
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
#include <typeinfo>

#include "args_handling.h"
#include "sequence_io.h"
#include "printing.h"
#include "typename.h"


namespace mc {

using std::cout;
using std::cerr;
using std::endl;


/*************************************************************************//**
 *
 *
 *****************************************************************************/
database
make_database(const args_parser& args,
              database::scope scope = database::scope::everything)
{
    auto dbfilename = database_name(args);
    if(dbfilename.empty()) {
        throw std::invalid_argument{"No database filename provided."};
    }
    return make_database<database>(dbfilename, scope);
}



/*************************************************************************//**
 *
 *
 *****************************************************************************/
void show_database_config(const args_parser& args)
{
    print_static_properties(make_database(args, database::scope::metadata_only));
    cout << endl;
}



/*************************************************************************//**
 *
 *
 *****************************************************************************/
void print_query_config()
{
    cout << "hit classifier       "
         << type_name<classification_candidates>() << '\n'
         << "------------------------------------------------\n"
         << endl;
}



/*************************************************************************//**
 *
 *
 *****************************************************************************/
void show_database_statistics(const args_parser& args)
{
    auto db = make_database(args);
    print_static_properties(db);
    print_content_properties(db);
}




/*************************************************************************//**
 *
 *
 *****************************************************************************/
void show_feature_map(const args_parser& args)
{
    auto db = make_database(args);
    print_static_properties(db);
    print_content_properties(db);
    cout << "===================================================\n";
    db.print_feature_map(cout);
    cout << "===================================================\n";
}



/*************************************************************************//**
 *
 *
 *****************************************************************************/
void show_feature_counts(const args_parser& args)
{
    auto db = make_database(args);
    print_static_properties(db);
    print_content_properties(db);
    cout << "===================================================\n";
    db.print_feature_counts(cout);
    cout << "===================================================\n";
}



/*************************************************************************//**
 *
 * @brief
 *
 *****************************************************************************/
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


/*************************************************************************//**
 *
 * @brief
 *
 *****************************************************************************/
void show_target_info(const args_parser& args)
{
    auto sids = std::vector<std::string>{};
    for(std::size_t i = 3; i < args.non_prefixed_count(); ++i) {
        sids.push_back(args.non_prefixed(i));
    }

    auto db = make_database(args, database::scope::metadata_only);

    if(!sids.empty()) {
        for(const auto& sid : sids) {
            const taxon* tax = db.taxon_with_name(sid);
            if(tax) {
                show_target_info(cout, db, *tax);
            }
            else {
                cout << "Target (reference sequence) " << sid
                     << " not found in database.\n";
            }
        }
    }
    else {
        cout << "Targets (reference sequences) in database:\n";
        for(const auto& tax : db.target_taxa()) {
            show_target_info(cout, db, tax);
        }
    }
}



/*************************************************************************//**
 *
 * @brief
 *
 *****************************************************************************/
void show_lineage_table(const args_parser& args)
{
    using rank = taxonomy::rank;

    auto db = make_database(args, database::scope::metadata_only);
    if(db.target_count() < 1) return;

    //table header
    cout << "name";
    for(auto r = rank::Sequence; r <= rank::Domain; ++r) {
        cout << '\t' << taxonomy::rank_name(r);
    }
    cout << '\n';

    //rows
    for(const auto& tax : db.target_taxa()) {
        cout << tax.name();
        auto ranks = db.ranks(tax);
        for(auto r = rank::Sequence; r <= rank::Domain; ++r) {
            cout << '\t'
                 << (ranks[int(r)] ? ranks[int(r)]->id() : taxonomy::none_id());
        }
        cout << '\n';
    }
}




/*************************************************************************//**
 *
 * @brief
 *
 *****************************************************************************/
void show_rank_statistics(const args_parser& args)
{
    auto rankName = args.non_prefixed(3);
    auto rank = taxonomy::rank_from_name(rankName);
    if(rank == taxonomy::rank::none) {
        cerr << "rank not recognized" << endl;
        return;
    }

    auto db = make_database(args, database::scope::metadata_only);

    std::map<const taxon*, std::size_t> stat;

    for(const auto& tax : db.target_taxa()) {
        const taxon* t = db.ranks(tax)[int(rank)];
        if(t) {
            auto it = stat.find(t);
            if(it != stat.end()) {
                ++(it->second);
            } else {
                stat.insert({t, 1});
            }
        }
    }

    cout << "Sequence distribution for rank " << rankName << ":\n"
         << "taxid \t taxon_name \t sequences" << endl;

    for(const auto& s : stat) {
        cout << s.first->id() << " \t "
             << s.first->name() << " \t " << s.second << '\n';
    }
}



/*************************************************************************//**
 *
 * @brief
 *
 *****************************************************************************/
void show_basic_exec_info()
{
    database db;
    print_static_properties(db);
    print_query_config();
    cout << endl;
}



/*************************************************************************//**
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
            print_query_config();
        }
        catch(std::exception& e) {
            cout << '\n' << e.what() << "!\n";
        }
    }
    else {
        const auto mode = args.non_prefixed(2);

        if(mode == "target" || mode == "targets"  ||
           mode == "sequ"   || mode == "sequence" || mode == "sequences")
        {
            show_target_info(args);
        }
        else if(mode == "lineages" || mode == "lineage" || mode == "lin")
        {
            show_lineage_table(args);
        }
        else if(mode == "rank") {
            if(args.non_prefixed_count() > 3) {
                show_rank_statistics(args);
            } else {
                cerr << "Please specify a taxonomic rank:\n";
                for(auto r = taxon_rank::Sequence; r <= taxon_rank::Domain; ++r) {
                    cerr << "    " << taxonomy::rank_name(r) << '\n';
                }
            }
        }
        else if(mode == "statistics" || mode == "stat") {
            show_database_statistics(args);
        }
        else if(mode == "features" || mode == "featuremap") {
            show_feature_map(args);
        }
        else if(mode == "featurecounts" || mode == "featuresizes") {
            show_feature_counts(args);
        }
        else {
            cerr << "Info mode '" << mode << "' not recognized.\n"
                      << "Properties of database " << args.non_prefixed(1)
                      << ":" << endl;
            show_database_config(args);
        }
    }
}


} // namespace mc

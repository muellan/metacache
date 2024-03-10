/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2024 André Müller (muellan@uni-mainz.de)
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


#include "candidate_generation.h"
#include "database.h"
#include "options.h"
#include "printing.h"
#include "typename.h"

#include <array>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>


namespace mc {

using std::cout;
using std::cerr;
using std::endl;
using std::string;


/*************************************************************************//**
 *
 *
 *****************************************************************************/
void show_database_config(const string& dbfile)
{
    int dbPart = -1;
    auto db = make_database(dbfile, dbPart, database::scope::metadata_only);
    print_static_properties(db);
}



/*************************************************************************//**
 *
 *
 *****************************************************************************/
void print_query_config()
{
    cout << "hit classifier       "
         << type_name<classification_candidates>() << '\n'
         << "------------------------------------------------\n";
}



/*************************************************************************//**
 *
 *
 *****************************************************************************/
void show_database_statistics(const string& dbfile, int dbPart)
{
    auto db = make_database(dbfile, dbPart);
    print_static_properties(db);
    print_content_properties(db);
}




/*************************************************************************//**
 *
 *
 *****************************************************************************/
void show_feature_map(const string& dbfile, int dbPart)
{
    auto db = make_database(dbfile, dbPart);
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
void show_feature_counts(const string& dbfile, int dbPart)
{
    auto db = make_database(dbfile, dbPart);
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
void show_target_info(std::ostream& os, const taxon& tax, const ranked_lineage& lineage)
{
    os  << "Target " << tax.name() << "):\n"
        << "    source:     "
        << tax.source().filename << " / " << tax.source().index
        << "\n    length:     " << tax.source().windows << " windows";

    for (const taxon* t : lineage) {
        if (t) {
            auto rn = string(t->rank_name()) + ":";
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
void show_target_info(const info_options& opt)
{
    int dbPart = -1;
    auto db = make_database(opt.dbfile, dbPart, database::scope::metadata_only);
    const auto& taxonomy = db.taxo_cache();

    if (!opt.targetNames.empty()) {
        for (const auto& target : opt.targetNames) {
            const taxon* tax = taxonomy.taxon_with_name(target);
            if (tax && tax->id() < 0) {
                target_id tgt = -tax->id()-1;
                show_target_info(cout, *tax, taxonomy.cached_ranks(tgt));
            }
            else {
                cout << "Target (reference sequence) '" << target
                     << "' not found in database.\n";
            }
        }
    }
    else {
        cout << "Targets (reference sequences) in database:\n";
        for (const auto& ranks : taxonomy.target_lineages()) {
            show_target_info(cout, *ranks[0], ranks);
        }
    }
}



/*************************************************************************//**
 *
 * @brief
 *
 *****************************************************************************/
void show_lineage_table(const info_options& opt)
{
    using rank = taxonomy::rank;

    int dbPart = -1;
    auto db = make_database(opt.dbfile, dbPart, database::scope::metadata_only);
    if (db.target_count() < 1) return;

    //table header
    cout << "name";
    for (auto r = rank::Sequence; r <= rank::Domain; ++r) {
        cout << '\t' << taxonomy::rank_name(r);
    }
    cout << '\n';

    const auto& taxonomy = db.taxo_cache();

    //rows
    for (const auto& ranks : taxonomy.target_lineages()) {
        cout << ranks[0]->name();
        for (auto r = rank::Sequence; r <= rank::Domain; ++r) {
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
void show_rank_statistics(const info_options& opt)
{
    if (opt.rank == taxon_rank::none) {
        cerr << "Please specify a taxonomic rank:\n";
        for (auto r = taxon_rank::Sequence; r <= taxon_rank::Domain; ++r) {
            cerr << "    " << taxonomy::rank_name(r) << '\n';
        }
        return;
    }

    int dbPart = -1;
    auto db = make_database(opt.dbfile, dbPart, database::scope::metadata_only);
    const auto& taxonomy = db.taxo_cache();

    std::map<const taxon*, std::size_t> stat;

    for (target_id tgt = 0; tgt < db.target_count(); ++tgt) {
        const taxon* t = taxonomy.cached_ancestor(tgt, opt.rank);
        if (t) {
            auto it = stat.find(t);
            if (it != stat.end()) {
                ++(it->second);
            } else {
                stat.insert({t, 1});
            }
        }
    }

    cout << "Sequence distribution for rank '"
         << taxonomy::rank_name(opt.rank) << "':\n"
         << "taxid \t taxon_name \t sequences" << endl;

    for (const auto& s : stat) {
        cout << s.first->id() << " \t "
             << s.first->name() << " \t "
             << s.second << '\n';
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
void main_mode_info(const cmdline_args& args)
{
    auto opt = get_info_options(args);

    switch (opt.mode) {
        default:
        case info_mode::basic:
            show_basic_exec_info();
            break;
        case info_mode::targets:
            show_target_info(opt);
            break;
        case info_mode::tax_lineages:
            show_lineage_table(opt);
            break;
        case info_mode::tax_ranks:
            show_rank_statistics(opt);
            break;
        case info_mode::db_config:
            show_database_config(opt.dbfile);
            print_query_config();
            break;
        case info_mode::db_statistics:
            show_database_statistics(opt.dbfile, opt.dbpart);
            break;
        case info_mode::db_feature_map:
            show_feature_map(opt.dbfile, opt.dbpart);
            break;
        case info_mode::db_feature_counts:
            show_feature_counts(opt.dbfile, opt.dbpart);
            break;
    }
}


} // namespace mc

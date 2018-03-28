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

#include <cstdint>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "args_parser.h"
#include "filesys_utility.h"
#include "cmdline_utility.h"
#include "sequence_io.h"
#include "modes.h"


namespace mc {

using std::cout;
using std::cerr;
using std::endl;
using std::flush;


/*************************************************************************//**
 *
 *
 *
 *****************************************************************************/
using taxid_t = std::uint64_t;

enum class annotation_mode
{
    none, taxid
};



/*************************************************************************//**
 *
 * @brief database configuration parameters
 *
 *****************************************************************************/
struct annotation_options
{
    annotation_mode mode = annotation_mode::none;

    std::string infile = "";
    std::string outfile = "";
    std::string fieldSeparator = " ";
    std::string valueSeparator = "|";

    sequence_id_type idtype = sequence_id_type::acc_ver;

    std::vector<std::string> mappingFiles;
};



/*************************************************************************//**
 *
 *
 *
 *****************************************************************************/
std::vector<std::string>
mapping_filenames(const args_parser& args)
{
    //first 2 non-prefixed args are:
    //metacache mode, mapping mode and file to annotate
    //all other non-prefixed args that are not
    //preceded by options ("-xyz") should be mapping file or folder names
    auto n = args.non_prefixed_count();
    if(n < 4) {
        throw std::invalid_argument{"No taxonomic mapping filenames provided"};
    }

    auto files = std::vector<std::string>{};
    files.reserve(n-3);
    for(std::size_t i = 3; i < n; ++i) {
        if(args.is_preceded_by_prefixed_arg(i)) break;
        auto name = args.non_prefixed(i);

        auto fnames = files_in_directory(name);

        if(!fnames.empty()) {
            files.insert(files.end(), std::make_move_iterator(fnames.begin()),
                                      std::make_move_iterator(fnames.end()));
        } else {
            files.push_back(name);
        }
    }

    return files;
}



/*************************************************************************//**
 *
 *
 *
 *****************************************************************************/
annotation_options
get_annotation_options(const args_parser& args)
{
    const annotation_options defaults;
    annotation_options opt;

    if(args.non_prefixed_count() < 2) {
        throw std::invalid_argument{"No annotation mode provided."};
    }

    auto modestr = args.non_prefixed(1);
    if(modestr.empty()) {
        throw std::invalid_argument{"No annotation mode provided."};
    }

    if(modestr == "taxid") {
        opt.mode = annotation_mode::taxid;
    }
    else {
        throw std::invalid_argument{
            "Annotation mode '" + modestr + "' not recognized."};
    }

    if(args.non_prefixed_count() < 3) {
        throw std::invalid_argument{"No input filename provided."};
    }

    opt.infile = args.non_prefixed(2);

    opt.outfile = args.get<std::string>("out", defaults.outfile);

    opt.mappingFiles = mapping_filenames(args);

    opt.fieldSeparator = args.get<std::string>({"field-sep", "field-separator"},
                                                 defaults.fieldSeparator);

    opt.valueSeparator = args.get<std::string>({"value-sep", "value-separator"},
                                                 defaults.valueSeparator);

    if(args.contains("id=accver"))
        opt.idtype = sequence_id_type::acc_ver;
    else if(args.contains("id=acc"))
        opt.idtype = sequence_id_type::acc;
    else if(args.contains("id=gi"))
        opt.idtype = sequence_id_type::gi;

    return opt;
}



/*************************************************************************//**
 *
 *
 *
 *****************************************************************************/
void read_mappings_from_file(sequence_id_type idtype,
                             const std::string& file,
                             std::map<std::string,taxid_t>& map)
{
    std::ifstream is {file};

    if(!is.good()) {
        cerr << "File '" << file << "' could not be read." << endl;
        return;
    }

    cout << "Processing mappings in file '" << file << "'" << endl;

    const auto fsize = file_size(file);
    auto nextStatStep = fsize / 1000;
    auto nextStat = nextStatStep;
    bool showProgress = fsize > 100000000;

    if(showProgress) show_progress_indicator(cout, 0);

    std::string acc;
    std::string acc_ver;
    std::string gi;
    taxid_t taxid;

    //skip first line
    getline(is, acc);
    acc.clear();

    while(is >> acc >> acc_ver >> taxid >> gi) {

        switch(idtype) {
            case sequence_id_type::acc:     map.insert({acc, taxid}); break;
            case sequence_id_type::acc_ver: map.insert({acc_ver, taxid}); break;
            case sequence_id_type::gi:      map.insert({gi, taxid}); break;
            default: break;
        }

        if(showProgress) {
            auto pos = is.tellg();
            if(pos >= nextStat) {
                show_progress_indicator(cout, pos / float(fsize));
                nextStat = pos + nextStatStep;
            }
        }
    }

    if(showProgress) clear_current_line(cout);
}



/*************************************************************************//**
 *
 *
 *
 *****************************************************************************/
void annotate_with_taxid(const annotation_options& opt,
                         const std::map<std::string,taxid_t>& map,
                         std::istream& is,
                         std::ostream& os)
{
   bool firstLine = true;
   std::string line;
   do {
       getline(is, line);

       if(!line.empty()) {
           //if line is header
           if(line[0] == '>' || line[0] == '@') {
               std::string id;
               switch(opt.idtype) {
                   case sequence_id_type::acc:
                       id = extract_ncbi_accession_number(line);
                       break;
                   case sequence_id_type::acc_ver:
                       id = extract_ncbi_accession_version_number(line);
                       break;
                   case sequence_id_type::gi:
                       id = extract_genbank_identifier(line);
                       break;
                   default: break;
               }

               //delete old taxid if present
               {
                   auto j = line.find("taxid");
                   if(j != std::string::npos) {
                       auto k = line.find(opt.valueSeparator, j+4);
                       if(k != std::string::npos) {
                           auto l = line.find(opt.fieldSeparator, k+1);
                           if(l == std::string::npos)
                               l = line.find(opt.valueSeparator, k+1);
                           line.erase(j, l-j+1);
                       }
                   }
               }

               if(!id.empty()) {
                   taxid_t taxid = 0;
                   auto it = map.find(id);
                   if(it != map.end()) {
                       taxid = it->second;
                   }

                   //insert taxid after first field separator / end of line
                   auto ipos = line.find(opt.fieldSeparator);
                   if(ipos == std::string::npos) {
                       line += opt.fieldSeparator + "taxid" +
                               opt.valueSeparator + std::to_string(taxid) +
                               opt.fieldSeparator;
                   }
                   else {
                       line.insert(ipos+1, "taxid" +
                                   opt.valueSeparator + std::to_string(taxid) +
                                   opt.fieldSeparator);
                   }
               }
           }
       }

       if(firstLine) {
           os << line;
           firstLine = false;
       }
       else if(!line.empty()) {
           os << '\n' << line;
       }
   } while(is.good() && os.good() && !line.empty());

}



/*************************************************************************//**
 *
 *
 *
 *****************************************************************************/
void annotate_with_taxid(const annotation_options& opt)
{

    std::ifstream is {opt.infile};
    if(!is.good()) {
        cerr << "Input file " << opt.infile << " could not be opened." << endl;
        return;
    }

    auto map = std::map<std::string,taxid_t>{};

    for(const auto& file : opt.mappingFiles) {
        read_mappings_from_file(opt.idtype, file, map);
    }

    if(opt.outfile.empty()) {
        annotate_with_taxid(opt, map, is, cout);
    }
    else {
        std::ofstream os {opt.outfile};
        if(os.good()) {
            cout << "Mapping ... " << flush;
            annotate_with_taxid(opt, map, is, os);
            cout << "complete." << endl;
        }
        else {
            cerr << "Output file " << opt.outfile
                 << " could not be opened."
                 << endl;
        }
    }

}



/*************************************************************************//**
 *
 *
 *
 *****************************************************************************/
void main_mode_annotate(const args_parser& args)
{
    try {
        auto opt = get_annotation_options(args);

        switch(opt.mode) {
            case annotation_mode::taxid: {
                cout << "Annotating sequences in '"
                     << opt.infile << "' with taxonomic ids.\n"
                     << "Output will be written to ";

                if(!opt.outfile.empty()) {
                    cout << "'" << opt.outfile << "'" << endl;
                } else {
                    cout << "'stdout'" << endl;
                }
                annotate_with_taxid(opt);
                break;
            }
            case annotation_mode::none:
                break;
        }
    }
    catch(std::exception& e) {
        cout << e.what() << "\n\n";
        show_help_for_topic("annotate");
    }
}


} // namespace mc

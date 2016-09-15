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
#include <cstdint>
#include <fstream>

#include "timer.h"
#include "args_handling.h"
#include "filesys_utility.h"
#include "cmdline_utility.h"
#include "sequence_io.h"

#include "modes_common.h"
#include "modes.h"


namespace mc {

/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
using taxid_t = std::uint64_t;

enum class annotation_mode
{
    none, taxid
};



/*****************************************************************************
 *
 * @brief database configuration parameters
 *
 *****************************************************************************/
struct annotation_param
{
    annotation_mode mode = annotation_mode::none;

    std::string infile = "";
    std::string outfile = "";
    std::string fieldSeparator = " ";
    std::string valueSeparator = "|";

    sequence_id_type idtype = sequence_id_type::acc;

    std::vector<std::string> mappingFiles;
};



/*****************************************************************************
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
    for(int i = 3; i < n; ++i) {
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



/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
annotation_param
get_annotation_param(const args_parser& args)
{
    const annotation_param defaults;
    annotation_param param;

    if(args.non_prefixed_count() < 2) {
        throw std::invalid_argument{"No annotation mode provided."};
    }

    auto modestr = args.non_prefixed(1);
    if(modestr.empty()) {
        throw std::invalid_argument{"No annotation mode provided."};
    }

    if(modestr == "taxid") {
        param.mode = annotation_mode::taxid;
    }
    else {
        throw std::invalid_argument{
            "Annotation mode '" + modestr + "' not recognized."};
    }

    if(args.non_prefixed_count() < 3) {
        throw std::invalid_argument{"No input filename provided."};
    }

    param.infile = args.non_prefixed(2);

    param.outfile = args.get<std::string>("out", defaults.outfile);

    param.mappingFiles = mapping_filenames(args);

    param.fieldSeparator = args.get<std::string>("field_separator",
                                                 defaults.fieldSeparator);

    param.valueSeparator = args.get<std::string>("value_separator",
                                                 defaults.valueSeparator);

    if(args.contains("id=accver"))
        param.idtype = sequence_id_type::acc_ver;
    else if(args.contains("id=acc"))
        param.idtype = sequence_id_type::acc;
    else if(args.contains("id=gi"))
        param.idtype = sequence_id_type::gi;

    return param;
}



/*****************************************************************************
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
        std::cout << "File '" << file << "' could not be read." << std::endl;
        return;
    }

    std::cout << "Processing mappings in file '" << file << "'" << std::endl;

    const auto fsize = file_size(file);
    auto nextStatStep = fsize / 1000;
    auto nextStat = nextStatStep;
    bool showProgress = fsize > 100000000;

    if(showProgress) show_progress_indicator(0);

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
                show_progress_indicator(pos / float(fsize));
                nextStat = pos + nextStatStep;
            }
        }
    }

    if(showProgress) clear_current_line();
}



/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
void annotate_with_taxid(const annotation_param& param,
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
               switch(param.idtype) {
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
                       auto k = line.find(param.valueSeparator, j+4);
                       if(k != std::string::npos) {
                           auto l = line.find(param.fieldSeparator, k+1);
                           if(l == std::string::npos)
                               l = line.find(param.valueSeparator, k+1);
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
                   auto ipos = line.find(param.fieldSeparator);
                   if(ipos == std::string::npos) {
                       line += param.fieldSeparator + "taxid" +
                               param.valueSeparator + std::to_string(taxid) +
                               param.fieldSeparator;
                   }
                   else {
                       line.insert(ipos+1, "taxid" +
                                   param.valueSeparator + std::to_string(taxid) +
                                   param.fieldSeparator);
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



/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
void annotate_with_taxid(const annotation_param& param)
{

    std::ifstream is {param.infile};
    if(!is.good()) {
        std::cout << "Input file " << param.infile << " could not be opened."
                  << std::endl;
    }

    auto map = std::map<std::string,taxid_t>{};

    for(const auto& file : param.mappingFiles) {
        read_mappings_from_file(param.idtype, file, map);
    }

    if(param.outfile.empty()) {
        annotate_with_taxid(param, map, is, std::cout);
    }
    else {
        std::ofstream os {param.outfile};
        if(os.good()) {
            std::cout << "Mapping ... " << std::flush;
            annotate_with_taxid(param, map, is, os);
            std::cout << "complete." << std::endl;
        }
        else {
            std::cout << "Output file " << param.outfile
                << " could not be opened."
                << std::endl;
        }
    }

}



/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
void main_mode_annotate(const args_parser& args)
{
    try {
        auto param = get_annotation_param(args);

        switch(param.mode) {
            case annotation_mode::taxid: {
                std::cout
                    << "Annotating sequences in '"
                    << param.infile << "' with taxonomic ids.\n"
                    << "Output will be written to ";

                if(!param.outfile.empty()) {
                    std::cout << "'" << param.outfile << "'" << std::endl;
                } else {
                    std::cout << "'stdout'" << std::endl;
                }
                annotate_with_taxid(param);
                break;
            }
            case annotation_mode::none:
                break;
        }
    }
    catch(std::exception& e) {
        std::cout << e.what() << "\n\n";
        show_help_for_topic("annotate");
    }
}


} // namespace mc

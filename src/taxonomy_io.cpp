/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2020 André Müller (muellan@uni-mainz.de)
 *                       & Robin Kobus  (kobus@uni-mainz.de)
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

#include <fstream>
#include <limits>
#include <sstream>

#include "taxonomy_io.h"
#include "cmdline_utility.h"
#include "filesys_utility.h"


namespace mc {


using std::string;
using std::cout;
using std::cerr;

using taxon_id = taxonomy::taxon_id;


/*************************************************************************//**
 *
 * @brief forward stream past the next occurence of a specific character
 *
 *****************************************************************************/
inline void
forward(std::istream& is, char c)
{
    is.ignore(std::numeric_limits<std::streamsize>::max(), c);
}



//-------------------------------------------------------------------
taxonomy
make_taxonomic_hierarchy(const string& taxNodesFile,
                         const string& taxNamesFile,
                         const string& mergeTaxFile,
                         info_level infoLvl)
{
    using taxon_id = taxonomy::taxon_id;

    const bool showInfo = infoLvl != info_level::silent;

    //read scientific taxon names
    //failure to do so will not be fatal
    auto taxonNames = std::map<taxon_id,string>{};

    std::ifstream is{taxNamesFile};
    if(is.good()) {
        if(showInfo) cout << "Reading taxon names ... " << std::flush;
        taxon_id lastId = 0;
        taxon_id taxonId = 0;
        string category;

        while(is.good()) {
            is >> taxonId;
            if(taxonId != lastId) {
                forward(is, '|');
                string name;
                is >> name;
                string word;
                while(is.good()) {
                    is >> word;
                    if(word.find('|') != string::npos) break;
                    name += " " + word;
                };
                forward(is, '|');
                is >> category;
                if(category.find("scientific") != string::npos) {
                    lastId = taxonId;
                    taxonNames.insert({taxonId,std::move(name)});
                }
            }
            forward(is, '\n');
        }
        if(showInfo) cout << "done." << std::endl;
    }
    else {
        cerr << "Could not read taxon names file "
                  << taxNamesFile
                  << "; continuing with ids only." << std::endl;
    }

    taxonomy tax;

    //read merged taxa
    is.close();
    is.open(mergeTaxFile);
    auto mergedTaxa = std::map<taxon_id,taxon_id>{};
    if(is.good()) {
        if(showInfo) cout << "Reading taxonomic node mergers ... " << std::flush;
        taxon_id oldId;
        taxon_id newId;

        while(is.good()) {
            is >> oldId;
            forward(is, '|');
            is >> newId;
            forward(is, '\n');
            mergedTaxa.insert({oldId,newId});

            tax.emplace(oldId, newId, "", "");
        }
        if(showInfo) cout << "done." << std::endl;
    }

    //read taxonomic structure
    is.close();
    is.open(taxNodesFile);
    if(is.good()) {
        if(showInfo) cout << "Reading taxonomic tree ... " << std::flush;
        taxon_id taxonId;
        taxon_id parentId;
        string rankName;
        string rankNameExt;

        while(is.good()) {
            is >> taxonId;
            forward(is, '|');
            is >> parentId;
            forward(is, '|');
            is >> rankName;
            is >> rankNameExt;
            if(rankNameExt != "|")
                rankName += ' ' + rankNameExt;
            forward(is, '\n');

            //get taxon name
            auto it = taxonNames.find(taxonId);
            auto taxonName = (it != taxonNames.end())
                             ? it->second : string("--");
            if(taxonName.empty()) {
                taxonName = "<" + std::to_string(taxonId) + ">";
            }

            //replace ids with new ids according to mergers
            //TODO this is stupid, handle mergers properly
            auto mi = mergedTaxa.find(taxonId);
            if(mi != mergedTaxa.end()) taxonId = mi->second;
            mi = mergedTaxa.find(parentId);
            if(mi != mergedTaxa.end()) parentId = mi->second;

            tax.emplace(taxonId, parentId, taxonName, rankName);
        }
        if(showInfo) cout << tax.size() << " taxa read." << std::endl;
    }
    else {
        cerr << "Could not read taxonomic nodes file "
             << taxNodesFile << std::endl;
        return tax;
    }

    //set rank of root
    tax.reset_rank(1, taxonomy::rank::root);

    //make sure every taxon has a rank designation
//    tax.rank_all_unranked();

    return tax;
}



//-------------------------------------------------------------------
void read_sequence_to_taxon_id_mapping(const string& mappingFile,
                                       std::map<string,taxon_id>& map,
                                       info_level infoLvl)
{
    const bool showInfo = infoLvl != info_level::silent;

    std::ifstream is{mappingFile};
    if(is.good()) {
        if(showInfo) {
            cout << "Reading sequence to taxon mappings from " << mappingFile
                 << std::endl;
        }

        const auto fsize = file_size(mappingFile);

        bool showProgress = showInfo && fsize > 100000000;
        //assembly_summary files have up to 550K lines
        //update progress indicator every 128K lines
        size_t step = 0;
        size_t statStep = 1UL << 17;
        if(showProgress) show_progress_indicator(cout, 0);


        //read first line(s) and determine the columns which hold
        //sequence ids (keys) and taxon ids
        int headerRow = 0;
        {
            string line;
            for(int i = 0; i < 10; ++i, ++headerRow) {
                getline(is, line);
                if(line[0] != '#') break;
            }
            if(headerRow > 0) --headerRow;
        }

        //reopen and forward to header row
        is.close();
        is.open(mappingFile);
        {
            string line;
            for(int i = 0; i < headerRow; ++i) getline(is, line);
        }

        if(is.good()) {
            //process header row
            int keycol = 0;
            int taxcol = 0;
            {
                int col = 0;
                string header;
                getline(is, header);
                std::istringstream hs(header);
                //get rid of comment chars
                hs >> header;
                while(hs >> header) {
                    if(header == "taxid") {
                        taxcol = col;
                    }
                    else if (header == "accession.version" ||
                            header == "assembly_accession")
                    {
                        keycol = col;
                    }
                    ++col;
                }
            }
            //taxid column assignment not found
            if(taxcol < 1) {
                //reopen file and use 1st column as key and 2nd column as taxid
                is.close();
                is.open(mappingFile);
                taxcol = 1;
            }

            string key;
            taxon_id taxonId;
            while(is.good()) {
                //forward to column with key
                for(int i = 0; i < keycol; ++i) forward(is, '\t');
                is >> key;
                //forward to column with taxid
                for(int i = 0; i < taxcol; ++i) forward(is, '\t');
                is >> taxonId;
                forward(is, '\n');

                map.insert({key, taxonId});

                if(showProgress && !(++step % statStep)) {
                    auto pos = is.tellg();
                    show_progress_indicator(cout, pos / float(fsize));
                }
            }
            if(showProgress) clear_current_line(cout);
        }
    }
}




//-------------------------------------------------------------------
std::map<string,taxon_id>
make_sequence_to_taxon_id_map(const std::vector<string>& localMappingFilenames,
                              const std::vector<string>& globalMappingFilenames,
                              const std::vector<string>& infilenames,
                              info_level infoLvl)
{
    //gather all taxonomic mapping files that can be found in any
    //of the input directories
    auto indirs = unique_directories(infilenames);

    auto map = std::map<string,taxon_id>{};

    for(const auto& dir : indirs) {
        for(const auto& file : localMappingFilenames) {
            read_sequence_to_taxon_id_mapping(dir + "/" + file, map, infoLvl);
        }
    }

    for(const auto& file : globalMappingFilenames) {
        read_sequence_to_taxon_id_mapping(file, map, infoLvl);
    }

    return map;
}



} // namespace mc

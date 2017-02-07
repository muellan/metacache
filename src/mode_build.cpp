/*****************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
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

#include <set>

#include "timer.h"
#include "args_handling.h"
#include "filesys_utility.h"
#include "cmdline_utility.h"
#include "sequence_io.h"

#include "modes_common.h"


namespace mc {


/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
enum class build_info : unsigned char {
    silent, moderate, verbose
};



/*****************************************************************************
 *
 * @brief database creation parameters
 *
 *****************************************************************************/
struct build_options
{
    int kmerlen = 16;
    int sketchlen = 16;
    int winlen = 128;
    int winstride = 113;

    float maxLoadFactor = -1;           //< 0 : use database default
    int maxLocationsPerFeature = database::max_supported_locations_per_feature();

    taxon_rank removeAmbigFeaturesOnRank = taxon_rank::none;
    int maxTaxaPerFeature = 1;

    taxonomy_options taxonomy;

    build_info infoMode = build_info::moderate;

    std::string dbfile;
    std::vector<std::string> infiles;
};



/*****************************************************************************
 *
 * @brief command line args -> database creation parameters
 *
 *****************************************************************************/
build_options
get_build_options(const args_parser& args)
{
    const build_options defaults;

    build_options param;

    param.dbfile = database_name(args);

    param.infiles = sequence_filenames(args);

    if(args.contains("silent")) {
        param.infoMode = build_info::silent;
    } else if(args.contains("verbose")) {
        param.infoMode = build_info::verbose;
    }

    param.kmerlen   = args.get<int>({"kmerlen"}, defaults.kmerlen);
    param.sketchlen = args.get<int>({"sketchlen"}, defaults.sketchlen);
    param.winlen    = args.get<int>({"winlen"}, defaults.winlen);
    param.winstride = args.get<int>({"winstride"}, param.winlen - param.kmerlen + 1);

    param.maxLoadFactor = args.get<float>({"max-load-fac", "max_load_fac"},
                                          defaults.maxLoadFactor);

    param.maxLocationsPerFeature = args.get<int>({"max-locations-per-feature",
                                                  "max_locations_per_feature" },
                                                 defaults.maxLocationsPerFeature);

    param.removeAmbigFeaturesOnRank = taxonomy::rank_from_name(
        args.get<std::string>({"remove-ambig-features",
                               "remove_ambig_features"}, "domain"));

    param.maxTaxaPerFeature = args.get<int>({"max-ambig-per-feature",
                                             "max_ambig_per_feature"},
                                            defaults.maxTaxaPerFeature);

    param.taxonomy = get_taxonomy_options(args);

    return param;
}



/*****************************************************************************
 *
 * @brief reads taxonomic tree + names from NCBI's taxnonomy files
 *
 *****************************************************************************/
taxonomy
make_taxonomic_hierarchy(const std::string& taxNodesFile,
                         const std::string& taxNamesFile = "",
                         const std::string& mergeTaxFile = "")
{
    using taxon_id = taxonomy::taxon_id;

    //read scientific taxon names
    //failure to do so will not be fatal
    auto taxonNames = std::map<taxon_id,std::string>{};

    std::ifstream is{taxNamesFile};
    if(is.good()) {
        std::cout << "Reading taxon names ... " << std::flush;
        taxon_id lastId = 0;
        taxon_id taxonId = 0;
        std::string category;

        while(is.good()) {
            is >> taxonId;
            if(taxonId != lastId) {
                is.ignore(std::numeric_limits<std::streamsize>::max(), '|');
                std::string name;
                is >> name;
                std::string word;
                while(is.good()) {
                    is >> word;
                    if(word.find('|') != std::string::npos) break;
                    name += " " + word;
                };
                is.ignore(std::numeric_limits<std::streamsize>::max(), '|');
                is >> category;
                if(category.find("scientific") != std::string::npos) {
                    lastId = taxonId;
                    taxonNames.insert({taxonId,std::move(name)});
                }
            }
            is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        std::cout << "done." << std::endl;
    }
    else {
        std::cerr << "Could not read taxon names file "
                  << taxNamesFile
                  << "; continuing with ids only." << std::endl;
    }

    //read merged taxa
    is.close();
    is.open(mergeTaxFile);
    auto mergedTaxa = std::map<taxon_id,taxon_id>{};
    if(is.good()) {
        std::cout << "Reading taxonomic node mergers ... " << std::flush;
        taxon_id oldId;
        taxon_id newId;

        while(is.good()) {
            is >> oldId;
            is.ignore(std::numeric_limits<std::streamsize>::max(), '|');
            is >> newId;
            is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            mergedTaxa.insert({oldId,newId});
        }
        std::cout << "done." << std::endl;
    }

    //read taxonomic structure
    taxonomy tax;

    is.close();
    is.open(taxNodesFile);
    if(is.good()) {
        std::cout << "Reading taxonomic tree ... " << std::flush;
        taxon_id taxonId;
        taxon_id parentId;
        std::string rankName;

        while(is.good()) {
            is >> taxonId;
            is.ignore(std::numeric_limits<std::streamsize>::max(), '|');
            is >> parentId;
            is.ignore(std::numeric_limits<std::streamsize>::max(), '|');
            is >> rankName;
            is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            //get taxon name
            auto it = taxonNames.find(taxonId);
            auto taxonName = (it != taxonNames.end())
                             ? it->second : std::string("--");
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
        std::cout << tax.taxon_count() << " taxa read." << std::endl;
    }
    else {
        std::cerr << "Could not read taxonomic nodes file "
                  << taxNodesFile << std::endl;
        return tax;
    }

    //make sure every taxon has a rank designation
//    tax.rank_all_unranked();

    return tax;
}




/*****************************************************************************
 *
 * @brief adds taxonomic information to database
 *
 *****************************************************************************/
void load_taxonomy_into_database(database& db,
                                 const build_options& param)
{
    db.apply_taxonomy( make_taxonomic_hierarchy(param.taxonomy.nodesFile,
                                                param.taxonomy.namesFile,
                                                param.taxonomy.mergeFile) );

    std::cout << "Taxonomy applied to database." << std::endl;
}




/*****************************************************************************
 *
 * @brief Reads a text file that maps sequences (ids or filenames) to taxon ids.
 *     Originally intended for use with the NCBI RefSeq "assembly_summary.txt"
 *     files but can be used with any text file provided that:
 *       - All lines in the file have the identical number of
 *         tab separated columns.
 *       - Either the first column contains the key (sequence id or file name).
 *       - Either the first lines contain the token "taxid" in the
 *         column which holds the taxids or the file has only two columns.
 *
 *****************************************************************************/
void read_sequence_to_taxon_id_mapping(
    const std::string& mappingFile,
    std::map<std::string,database::taxon_id>& map)
{
    using taxon_id = database::taxon_id;

    std::ifstream is{mappingFile};
    if(is.good()) {
        std::cout << "Reading sequence to taxon mappings from " << mappingFile
                  << std::endl;

        const auto fsize = file_size(mappingFile);
        auto nextStatStep = fsize / 1000;
        auto nextStat = nextStatStep;
        bool showProgress = fsize > 100000000;
        if(showProgress) {
            show_progress_indicator(0);
        }

        //read first line(s) and determine the columns which hold
        //sequence ids (keys) and taxon ids
        int headerRow = 0;
        {
            std::string line;
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
            std::string line;
            for(int i = 0; i < headerRow; ++i) {
                getline(is, line);
            }
        }

        if(is.good()) {
            //process header row
            int keycol = 0;
            int taxcol = 0;
            {
                int col = 0;
                std::string header;
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

            std::string key;
            taxon_id taxonId;
            while(is.good()) {
                for(int i = 0; i < keycol; ++i) { //forward to column with key
                    is.ignore(std::numeric_limits<std::streamsize>::max(), '\t');
                }
                is >> key;
                for(int i = 0; i < taxcol; ++i) { //forward to column with taxid
                    is.ignore(std::numeric_limits<std::streamsize>::max(), '\t');
                }
                is >> taxonId;
                is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    //            std::cout << key << " -> " << taxonId << std::endl;

                map.insert({key, taxonId});

                if(showProgress) {
                    auto pos = is.tellg();
                    if(pos >= nextStat) {
                        show_progress_indicator(pos / float(fsize));
                        nextStat = pos + nextStatStep;
                    }
                }
            }
            if(showProgress) {
                clear_current_line();
            }
        }
    }
}




/*****************************************************************************
 *
 * @brief reads file -> taxon id mappings files from the same directories
 *        as the input files
 *
 *****************************************************************************/
std::map<std::string,database::taxon_id>
make_sequence_to_taxon_id_map(const std::vector<std::string>& mappingFilenames,
                              const std::vector<std::string>& infilenames)
{
    //gather all taxonomic mapping files that can be found in any
    //of the input directories
    auto indirs = unique_directories(infilenames);

    auto map = std::map<std::string,database::taxon_id>{};

    for(const auto& dir : indirs) {
        for(const auto& file : mappingFilenames) {
            read_sequence_to_taxon_id_mapping(dir + "/" + file, map);
        }
    }

    return map;
}




/*****************************************************************************
 *
 * @brief Alternative way of providing taxonomic mappings.
 *        The input file must be a text file with each line in the format:
 *        accession  accession.version taxid gi
 *        (like the NCBI's *.accession2version files)
 *
 *****************************************************************************/
void rank_targets_post_process(database& db,
                               std::set<target_id>& gids,
                               const std::string& mappingFile)
{
    if(gids.empty()) return;

    std::ifstream is {mappingFile};
    if(!is.good()) return;

    const auto fsize = file_size(mappingFile);
    auto nextStatStep = fsize / 1000;
    auto nextStat = nextStatStep;
    bool showProgress = fsize > 100000000;

    std::cout << "Try to map sequences to taxa using '" << mappingFile
              << "' (" << std::max(std::streamoff(1),
                                   fsize/(1024*1024)) << " MB)" << std::endl;

    if(showProgress) show_progress_indicator(0);

    std::string acc;
    std::string accver;
    std::uint64_t taxid;
    std::string gi;

    //skip header
    getline(is, acc);
    acc.clear();

    while(is >> acc >> accver >> taxid >> gi) {
        //target in database?
        //accession.version is the default
        auto tid = db.target_id_of_sequence(accver);
        if(!db.is_valid(tid)) {
            tid = db.target_id_of_sequence(acc);
            if(!db.is_valid(tid)) {
                tid = db.target_id_of_sequence(gi);
            }
        }

        //if in database then map to taxon
        if(db.is_valid(tid)) {
            auto it = gids.find(tid);
            if(it != gids.end()) {
                db.rank_target(tid, taxid);
                gids.erase(it);
                if(gids.empty()) break;
            }
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
std::set<target_id>
unranked_targets(const database& db)
{
    auto res = std::set<target_id>{};

    for(target_id i = 0; i < db.target_count(); ++i) {
        if(db.taxon_of_target(i).none()) {
             res.insert(i);
        }
    }

    return res;
}




/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
void try_to_rank_unranked_targets(database& db, const build_options& param)
{
    auto unranked = unranked_targets(db);

    if(!unranked.empty()) {
        std::cout << unranked.size()
                  << " targets could not be ranked." << std::endl;

        for(const auto& file : param.taxonomy.mappingPostFiles) {
            rank_targets_post_process(db, unranked, file);
        }
    }

    unranked = unranked_targets(db);
    if(unranked.empty()) {
        std::cout << "All targets are ranked." << std::endl;
    }
    else {
        std::cout << unranked.size()
                  << " targets remain unranked." << std::endl;
    }

}




/*****************************************************************************
 *
 * @brief
 *
 *****************************************************************************/
std::string
extract_sequence_id(const std::string& text)
{
    auto sequid = extract_ncbi_accession_version_number(text);
    if(!sequid.empty()) return sequid;

    sequid = extract_genbank_identifier(text);
    if(!sequid.empty()) return sequid;

    return extract_ncbi_accession_number(text);
}




/*****************************************************************************
 *
 * @brief adds reference sequences to database
 *
 *****************************************************************************/
void add_targets_to_database(database& db,
    const std::vector<std::string>& infiles,
    const std::map<std::string,database::taxon_id>& sequ2taxid,
    build_info infoMode = build_info::moderate)
{
    int n = infiles.size();
    int i = 0;

    //read sequences
    for(const auto& filename : infiles) {
        if(infoMode == build_info::verbose) {
            std::cout << "  " << filename << " ... " << std::flush;
        } else {
            show_progress_indicator(i/float(n));
        }

        try {
            auto reader = make_sequence_reader(filename);
            while(reader->has_next()) {
                database::sequence_origin origin;
                origin.filename = filename;
                origin.index = 0;

                auto sequ = reader->next();
                if(!sequ.data.empty()) {
                    auto seqId = extract_sequence_id(sequ.header);
                    auto fileId = extract_sequence_id(filename);

                    //make sure sequence id is not empty,
                    //use entire header if neccessary
                    if(seqId.empty()) {
                        if(!fileId.empty()) {
                            seqId = fileId;
                        } else {
                            seqId = sequ.header;
                        }
                    }

                    //targets need to have a sequence id
                    //look up taxon id
                    taxon_id taxid = 0;
                    if(!sequ2taxid.empty()) {
                        auto it = sequ2taxid.find(seqId);
                        if(it != sequ2taxid.end()) {
                            taxid = it->second;
                        } else {
                            it = sequ2taxid.find(fileId);
                            if(it != sequ2taxid.end()) taxid = it->second;
                        }
                    }
                    //no valid taxid assigned -> try to find one in annotation
                    if(taxid > 0) {
                        if(infoMode == build_info::verbose)
                            std::cout << "[" << seqId << ":" << taxid << "] ";
                    }
                    else {
                        taxid = extract_taxon_id(sequ.header);
                        if(infoMode == build_info::verbose)
                            std::cout << "[" << seqId << "] ";
                    }

                    //try to add to database
                    bool added = db.add_target(sequ.data, seqId, taxid, origin);
                    if(infoMode == build_info::verbose && !added) {
                        std::cout << seqId << " not added to database" << std::endl;
                    }
                }
                ++origin.index; //track sequence index in file
            }
            if(infoMode == build_info::verbose) {
                std::cout << "done." << std::endl;
            }
        }
        catch(database::target_limit_exceeded_error&) {
            std::cout << std::endl;
            std::cerr
                << "Reached maximum number of targets per database ("
                << db.max_target_count() << ").\n"
                << "See 'README.md' on how to compile MetaCache with "
                << "support for databases with more reference targets.\n"
                << std::endl;
            break;
        }
        catch(std::exception& e) {
            if(infoMode == build_info::verbose) {
                std::cout << "FAIL: " << e.what() << std::endl;
            }
        }

        ++i;
    }
}




/*****************************************************************************
 *
 * @brief prepares datbase for build, adds targets and writes database to disk
 *
 *****************************************************************************/
void add_to_database(database& db, const build_options& param)
{
    if(param.maxLocationsPerFeature > 0) {
        db.max_locations_per_feature(param.maxLocationsPerFeature);
    }

    if(param.maxLoadFactor > 0) {
        db.max_load_factor(param.maxLoadFactor);
    }

    if(!param.taxonomy.path.empty()) {
        load_taxonomy_into_database(db, param);
    }

    if(db.taxon_count() < 1) {
        std::cout << "The datbase doesn't contain a taxonomic hierarchy yet.\n"
                  << "You can add one or update later via:\n"
                  << "   ./metacache add <database> -taxonomy <directory>"
                  << std::endl;
    }

    if(param.removeAmbigFeaturesOnRank != taxon_rank::none) {
        if(db.taxon_count() > 1) {
            std::cout << "Ambiguous features on rank "
                      << taxonomy::rank_name(param.removeAmbigFeaturesOnRank)
                      << " will be removed afterwards.\n";
        } else {
            std::cout << "Could not determine amiguous features "
                      << "due to missing taxonomic information.\n";
        }
    }
    
    print_properties(db);

    timer time;
    time.start();

    if(!param.infiles.empty()) {
        std::cout << "\nProcessing reference sequences." << std::endl;
        const auto initNumTargets = db.target_count();

        add_targets_to_database(db,
            param.infiles,
            make_sequence_to_taxon_id_map(param.taxonomy.mappingPreFiles,
                                          param.infiles),
            param.infoMode);


        if(param.infoMode == build_info::moderate) {
            clear_current_line();
        }
        std::cout << "Added "
                  << (db.target_count() - initNumTargets) << " reference sequences "
                  << "in " << time.seconds() << " s" << std::endl;

    }

    try_to_rank_unranked_targets(db, param);

    if(param.removeAmbigFeaturesOnRank != taxon_rank::none &&
        db.taxon_count() > 1)
    {
        std::cout << "\nRemoving ambiguous features on rank "
                  << taxonomy::rank_name(param.removeAmbigFeaturesOnRank)
                  << "..." << std::endl;

        db.remove_ambiguous_features(param.removeAmbigFeaturesOnRank,
                                     param.maxTaxaPerFeature);

        print_properties(db);
        std::cout << '\n';
    }

    write_database(db, param.dbfile);

    time.stop();

    std::cout << "Total build time: " << time.seconds() << " s" << std::endl;

    //prevents slow deallocation
    db.clear_without_deallocation();
}




/*****************************************************************************
 *
 * @brief builds a database from reference input sequences
 *
 *****************************************************************************/
void main_mode_build(const args_parser& args)
{
    std::cout << "Building new database from reference sequences." << std::endl;

    auto param = get_build_options(args);

    if(param.infiles.empty()) {
        throw std::invalid_argument{
            "Nothing to do - no reference sequences provided."};
    }

    //configure sketching scheme
    auto sketcher = database::sketcher{};
    sketcher.kmer_size(param.kmerlen);
    sketcher.sketch_size(param.sketchlen);

    auto db = database{sketcher};
    db.target_window_size(param.winlen);
    db.target_window_stride(param.winstride);

    add_to_database(db, param);
}




/*****************************************************************************
 *
 * @brief adds reference sequences to an existing database
 *
 *****************************************************************************/
void main_mode_build_modify(const args_parser& args)
{
    auto param = get_build_options(args);

    std::cout << "Modify database " << param.dbfile << std::endl;

    auto db = make_database<database>(param.dbfile);

    if(!param.infiles.empty()) {
        std::cout << "Adding reference sequences to database..." << std::endl;
    }

    add_to_database(db, param);
}


} // namespace mc

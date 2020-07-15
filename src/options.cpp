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

/*************************************************************************//**
 *
 * @file Facilities that generate option objects from command line arguments.
 * Argument parsing and doc generation uses the 'CLIPP' library.
 * If an option has multiple associated flag strings like, e.g.,
 * "-no-map" or "-nomap" the first one will appear in the documentation
 * and is considered to be the canonical flag while following flags are
 * there to preserve backwards compatibility with older versions of MetaCach.
 *
 *****************************************************************************/

#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <regex>

#include "options.h"
#include "filesys_utility.h"
#include "database.h"

#include "../dep/clipp.h"


namespace mc {

using std::size_t;
using std::vector;
using std::string;
using std::to_string;

using namespace std::string_literals;



/*************************************************************************//**
 *
 * @brief collects all command line interface error messages
 *
 *****************************************************************************/
class error_messages {
public:
    error_messages& operator += (const string& message) {
        messages_.push_back(message);
        return *this;
    }
    error_messages& operator += (string&& message) {
        messages_.push_back(std::move(message));
        return *this;
    }

    bool any() const noexcept   { return !messages_.empty(); }

    string str() const {
        string s;
        for(const auto msg : messages_) {
            if(!msg.empty()) s += msg + '\n';
        }
        return s;
    }

//    auto begin() const noexcept { return messages_.begin(); }
//    auto end()   const noexcept { return messages_.end();   }

private:
    vector<string> messages_;
};




/*****************************************************************************
 *
 *  H E L P E R S
 *
 *****************************************************************************/
string taxon_rank_names(const string& separator = ", ")
{
    string s;
    for(auto r = taxon_rank::Sequence; r < taxon_rank::Domain; ++r) {
        s += taxonomy::rank_name(r);
        s += separator;
    }
    s += taxonomy::rank_name(taxon_rank::Domain);
    return s;
}



//-------------------------------------------------------------------
/// @return database filename with extension
string sanitize_database_name(string name)
{
    if(name.find(".db") == string::npos) {
        name += ".db";
    }
    return name;
}



//-------------------------------------------------------------------
/// @return replaces '\t' with tab char and remove other special chars
string sanitize_special_chars(const string& text)
{
    return std::regex_replace( std::regex_replace(text,
        // no newlines, vertical tabs, etc. allowed
        std::regex(R"((\\n)|(\\v)|(\\r)|(\\a))"), ""),
        // substitute literat "\t" with tab character
        std::regex(R"(\\t)"), "\t");
}


//-------------------------------------------------------------------
void replace_directories_with_contained_files(vector<string>& names)
{
    vector<string> result;
    result.reserve(names.size());

    for(const auto& name : names) {

        auto fnames = files_in_directory(name);

        if(!fnames.empty()) {
            result.insert(result.end(), std::make_move_iterator(fnames.begin()),
                                        std::make_move_iterator(fnames.end()));
        } else {
            result.push_back(name);
        }
    }

    std::swap(result, names);
}



//-------------------------------------------------------------------
auto cli_doc_formatting()
{
    return clipp::doc_formatting{}
        .first_column(0)
        .doc_column(22)
        .last_column(80)
        .indent_size(4)
        .line_spacing(1)
        .alternatives_min_split_size(2)
        .paragraph_spacing(2)
        .max_flags_per_param_in_usage(1)
        .max_flags_per_param_in_doc(1)
        ;
}

auto cli_usage_formatting()
{
    return cli_doc_formatting().first_column(4).line_spacing(0);
}



//-------------------------------------------------------------------
/// @brief adds 'parameter' that catches unknown args with '-' prefix
clipp::parameter
catch_unknown(error_messages& err) {
    return clipp::any(clipp::match::prefix{"-"},
        [&](const string& arg) { err += "unknown argument: "s + arg; });
}



//-------------------------------------------------------------------
taxon_rank rank_from_name(const string& name, error_messages& err)
{
    auto r = taxonomy::rank_from_name(name);

    if(r == taxon_rank::none) {
        err += "Unknown taxonomic rank '"s + name + "'!\n";
        err += "Valid rank names are:\n    " + taxon_rank_names("\n    ") + "\n";
    }

    return r;
}



//-------------------------------------------------------------------
void raise_default_error(const error_messages& err,
                         const string& mode = "",
                         const string& usage = "",
                         const string& examples = "")
{
    auto msg = err.str();

    if(!msg.empty())      msg += "\n";

    if(!usage.empty())    msg += "USAGE:\n" + usage + "\n\n";
    if(!examples.empty()) msg += "EXAMPLES:\n" + examples + "\n\n";

    if(!mode.empty()) {
        msg += "\nYou can view the full interface documentation of mode '"s
            + mode + "' with:\n    metacache help " + mode + " | less";
    }

    throw std::invalid_argument{std::move(msg)};
}




/*****************************************************************************
 *
 *
 *  S H A R E D
 *
 *
 *****************************************************************************/

/// @brief
clipp::parameter
database_parameter(string& filename, error_messages& err)
{
    using namespace clipp;

    return value(match::prefix_not{"-"}, "database")
        .call([&](const string& arg){ filename = sanitize_database_name(arg); })
        .if_missing([&]{ err += "Database filename is missing!"; })
        .doc("database file name;\n"
             "A MetaCache database contains taxonomic information and "
             "min-hash signatures of reference sequences "
             "(complete genomes, scaffolds, contigs, ...).\n");
}



//-------------------------------------------------------------------
clipp::group
info_level_cli(info_level& lvl, error_messages& err)
{
    using namespace clipp;

    return one_of (
        option("-silent").set(lvl, info_level::silent),
        option("-verbose").set(lvl, info_level::verbose)
        .if_conflicted([&]{
            err += "Info level must be either '-silent' or '-verbose'!";
        })
    )
        % "information level during build:\n"
          "silent => none / verbose => most detailed\n"
          "default: neither => only errors/important info";
}



//-------------------------------------------------------------------
/// @brief shared command-line options for taxonomy
clipp::group
taxonomy_cli(taxonomy_options& opt, error_messages& err)
{
    using namespace clipp;

    return (
    (   option("-taxonomy") &
        value("path", opt.path)
            .if_missing([&]{ err += "Taxonomy path is missing after '-taxonomy'!"; })
    )
        % "directory with taxonomic hierarchy data (see NCBI's taxonomic data files)\n"
    ,
    (   option("-taxpostmap") &
        values("file", opt.mappingPostFiles)
            .if_missing([&]{ err += "Taxonomy mapping files are missing after '-taxpostmap'!"; })
    )
        % "Files with sequence to taxon id mappings that are used "
          "as alternative source in a post processing step.\n"
          "default: 'nucl_(gb|wgs|est|gss).accession2taxid'"
    );

}



//-------------------------------------------------------------------
/// @brief shared command-line options for sequence sketching
clipp::group
sketching_options_cli(sketching_options& opt, error_messages& err)
{
    using namespace clipp;
    return (
    (   option("-kmerlen") &
        integer("k", opt.kmerlen)
            .if_missing([&]{ err += "Number missing after '-kmerlen'!"; })
    )
        %("number of nucleotides/characters in a k-mer\n"
          "default: "s + (opt.kmerlen > 0 ? to_string(opt.kmerlen)
                                          : "determined by database"s))
    ,
    (   option("-sketchlen") &
        integer("s", opt.sketchlen)
            .if_missing([&]{ err += "Number missing after '-sketchlen'!"; })
    )
        %("number of features (k-mer hashes) per sampling window\n"
          "default: "s + (opt.sketchlen > 0 ? to_string(opt.sketchlen)
                                            : "determined by database"s))
    ,
    (   option("-winlen") &
        integer("w", opt.winlen)
            .if_missing([&]{ err += "Number missing after '-winlen'!"; })
    )
        %("number of letters in each sampling window\n"
          "default: "s + (opt.winlen > 0 ? to_string(opt.winlen)
                                         : "determined by database"s))
    ,
    (   option("-winstride") &
        integer("l", opt.winstride)
            .if_missing([&]{ err += "Number missing after '-winstride'!"; })
    )
        %("distance between window starting positions\n"
          "default: "s +
          (opt.winlen > 0 && opt.kmerlen > 0
              ? (to_string(opt.winlen - opt.kmerlen + 1) + " (w-k+1)")
              : "determined by database"s)
    )
    );
}



//-------------------------------------------------------------------
/// @brief shared command-line options for sequence sketching
clipp::group
database_storage_options_cli(database_storage_options& opt, error_messages& err)
{
    using namespace clipp;

    const database defaultDb;

    return (
#ifndef GPU_MODE
    (   option("-parts") &
        integer("#", opt.numParts)
            .if_missing([&]{ err += "Number missing after '-parts'!"; })
    )
        %("Sets the number of database parts to use."
          "default: 1"s)
    ,
#else
    (   option("-gpus") &
        integer("#", opt.numParts)
            .if_missing([&]{ err += "Number missing after '-gpus'!"; })
    )
        %("Sets the maximum number of GPUs to use."
          "default: all available GPUs"s)
    ,
#endif
    (   option("-max-locations-per-feature") &
        integer("#", opt.maxLocationsPerFeature)
            .if_missing([&]{ err += "Number missing after '-max-locations-per-feature'!"; })
    )
        %("maximum number of reference sequence locations to be stored per feature;\n"
          "If the value is too high it will significantly impact querying speed. "
          "Note that an upper hard limit is always imposed by the data type "
          "used for the hash table bucket size (set with "
          "compilation macro '-DMC_LOCATION_LIST_SIZE_TYPE').\n"
          "default: "s + to_string(defaultDb.max_locations_per_feature()))
    ,
    (
        option("-remove-overpopulated-features")
            .set(opt.removeOverpopulatedFeatures)
        %("Removes all features that have reached the maximum allowed "
          "amount of locations per feature. This can improve querying "
          "speed and can be used to remove non-discriminative features.\n"
          "default: "s + (opt.removeOverpopulatedFeatures ? "on" : "off"))
    )
    ,
    (   option("-remove-ambig-features") &
        value("rank", [&](const string& name) {
                opt.removeAmbigFeaturesOnRank = rank_from_name(name, err);
            })
            .if_missing([&]{ err += "Taxonomic rank missing after '-remove-ambig-features'!"; })
    )
        %("Removes all features that have more distinct reference sequence "
          "on the given taxonomic rank than set by '-max-ambig-per-feature'. "
          "This can decrease the database size significantly at the expense "
          "of sensitivity. Note that the lower the given taxonomic rank is, "
          "the more pronounced the effect will be.\n"
          "Valid values: "s + taxon_rank_names() + "\n"s +
          "default: "s + (opt.removeAmbigFeaturesOnRank != taxon_rank::none ? "on" : "off"))
    ,
    (   option("-max-ambig-per-feature") &
        integer("#", opt.maxTaxaPerFeature)
            .if_missing([&]{ err += "Number missing after '-max-ambig-per-feature'!"; })
    )
        % "Maximum number of allowed different reference sequence taxa per feature "
          "if option '-remove-ambig-features' is used.\n"
    ,
    (   option("-max-load-fac", "-max-load-factor") &
        number("factor", opt.maxLoadFactor)
            .if_missing([&]{ err += "Number missing after '-max-load-fac'!"; })
    )
        %("maximum hash table load factor;\n"
          "This can be used to trade off larger memory consumption for "
          "speed and vice versa. A lower load factor will improve speed, "
          "a larger one will improve memory efficiency.\n"
          "default: "s + to_string(defaultDb.max_load_factor())
    )
    );
}



//-------------------------------------------------------------------
void augment_taxonomy_options(taxonomy_options& opt)
{
    if(!opt.path.empty() && opt.path.back() != '/') opt.path += '/';

    opt.nodesFile = opt.path + "nodes.dmp";
    opt.namesFile = opt.path + "names.dmp";
    opt.mergeFile = opt.path + "merged.dmp";

    opt.mappingPreFilesLocal.push_back("assembly_summary.txt");
    opt.mappingPreFilesGlobal.push_back(opt.path + "assembly_summary_refseq.txt");
    opt.mappingPreFilesGlobal.push_back(opt.path + "assembly_summary_refseq_historical.txt");
    opt.mappingPreFilesGlobal.push_back(opt.path + "assembly_summary_genbank.txt");
    opt.mappingPreFilesGlobal.push_back(opt.path + "assembly_summary_genbank_historical.txt");

    //default NCBI accession to taxon map file names
    opt.mappingPostFiles.push_back(opt.path + "nucl_gb.accession2taxid");
    opt.mappingPostFiles.push_back(opt.path + "nucl_wgs.accession2taxid");
    opt.mappingPostFiles.push_back(opt.path + "nucl_est.accession2taxid");
    opt.mappingPostFiles.push_back(opt.path + "nucl_gss.accession2taxid");

    //find additional maps by file extension ".accession2taxid"
    for(const auto f : files_in_directory(opt.path)) {
        if(f.find(".accession2taxid") != string::npos) {
            if(std::find(opt.mappingPostFiles.begin(),
                         opt.mappingPostFiles.end(), f)
               == opt.mappingPostFiles.end())
            {
                opt.mappingPostFiles.push_back(f);
            }
        }
    }
}





/*****************************************************************************
 *
 *
 *  B U I L D   M O D E
 *
 *
 *****************************************************************************/
/// @brief build mode command-line options
clipp::group
build_mode_cli(build_options& opt, error_messages& err)
{
    using namespace clipp;

    return (
    "REQUIRED PARAMETERS" %
    (
        database_parameter(opt.dbfile, err)
        ,
        values(match::prefix_not{"-"}, "sequence file/directory", opt.infiles)
            .if_missing([&]{
                err += "No reference sequence files provided or found!";
            })
            % "FASTA or FASTQ files containing genomic sequences "
              "(complete genomes, scaffolds, contigs, ...) that shall be"
              "used as representatives of an organism/taxon.\n"
              "If directory names are given, they will be searched for "
              "sequence files (at most 10 levels deep).\n"
    ),
    "BASIC OPTIONS" %
    (
        taxonomy_cli(opt.taxonomy, err),
        info_level_cli(opt.infoLevel, err)
    ),
    "SKETCHING (SUBSAMPLING)" %
        sketching_options_cli(opt.sketching, err)
    ,
    "ADVANCED OPTIONS" %
    (
        option("-reset-taxa", "-reset-parents").set(opt.resetParents)
            %("Attempts to re-rank all sequences after the main build phase "
              "using '.accession2taxid' files. This will reset the taxon id "
              "of a reference sequence even if a taxon id could be obtained "
              "from other sources during the build phase.\n"
              "default: "s + (opt.resetParents ? "on" : "off"))
        ,
        database_storage_options_cli(opt.dbconfig, err)
    ),
    catch_unknown(err)
    );

}



//-------------------------------------------------------------------
build_options
get_build_options(const cmdline_args& args, build_options opt)
{
    error_messages err;

    auto cli = build_mode_cli(opt, err);

    auto result = clipp::parse(args, cli);

    if(!result || err.any()) {
        raise_default_error(err, "build", build_mode_usage());
    }

    augment_taxonomy_options(opt.taxonomy);
    replace_directories_with_contained_files(opt.infiles);

    if(opt.dbconfig.maxLocationsPerFeature < 0)
        opt.dbconfig.maxLocationsPerFeature = database::max_supported_locations_per_feature();

    auto& sk = opt.sketching;
    if(sk.winstride < 0) sk.winstride = sk.winlen - sk.kmerlen + 1;

    return opt;
}



//-------------------------------------------------------------------
string build_mode_usage() {
    return
    "    metacache build <database> <sequence file/directory>... [OPTION]...\n\n"
    "    metacache build <database> [OPTION]... <sequence file/directory>...";
}



//-------------------------------------------------------------------
string build_mode_examples() {
    return
    "    Build database 'mydb' from sequence file 'genomes.fna':\n"
    "        metacache build mydb genomes.fna\n"
    "\n"
    "    Build database with latest complete genomes from the NCBI RefSeq\n"
    "        download-ncbi-genomes refseq/bacteria myfolder\n"
    "        download-ncbi-genomes refseq/viruses myfolder\n"
    "        download-ncbi-taxonomy myfolder\n"
    "        metacache build myRefseq myfolder -taxonomy myfolder\n"
    "\n"
    "    Build database 'mydb' from two sequence files:\n"
    "        metacache build mydb mrsa.fna ecoli.fna\n"
    "\n"
    "    Build database 'myBacteria' from folder containing sequence files:\n"
    "        metacache build myBacteria all_bacteria\n";
}



//-------------------------------------------------------------------
string build_mode_docs() {

    build_options opt;
    error_messages err;

    auto cli = build_mode_cli(opt, err);

    string docs = "SYNOPSIS\n\n";

    docs += build_mode_usage();

    docs += "\n\n\n"
        "DESCRIPTION\n"
        "\n"
        "    Create a new database of reference sequences (usually genomic sequences).\n"
        "\n\n";

    docs += clipp::documentation(cli, cli_doc_formatting()).str();

    docs += "\n\nEXAMPLES\n\n";
    docs += build_mode_examples();

    return docs;
}





/*****************************************************************************
 *
 *
 *  M O D I F Y   M O D E
 *
 *
 *****************************************************************************/

/// @brief build mode command-line options
clipp::group
modify_mode_cli(build_options& opt, error_messages& err)
{
    using namespace clipp;

    return (
    "REQUIRED PARAMETERS" %
    (
        database_parameter(opt.dbfile, err)
        ,
        values(match::prefix_not{"-"}, "sequence file/directory", opt.infiles)
            .if_missing([&]{
                err += "No reference sequence files provided or found!";
            })
            % "FASTA or FASTQ files containing genomic sequences "
              "(complete genomes, scaffolds, contigs, ...) that shall be"
              "used as representatives of an organism/taxon.\n"
              "If directory names are given, they will be searched for "
              "sequence files (at most 10 levels deep).\n"
    ),
    "BASIC OPTIONS" %
    (
        taxonomy_cli(opt.taxonomy, err),
        info_level_cli(opt.infoLevel, err)
    ),
    "ADVANCED OPTIONS" %
    (
        option("-reset-taxa", "-reset-parents")
            .set(opt.resetParents)
            %("Attempts to re-rank all sequences after the main build phase "
              "using '.accession2taxid' files. This will reset the taxon id "
              "of a reference sequence even if a taxon id could be obtained "
              "from other sources during the build phase.\n"
              "default: "s + (opt.resetParents ? "on" : "off"))
        ,
        database_storage_options_cli(opt.dbconfig, err)
    ),
    catch_unknown(err)
    );

}



//-------------------------------------------------------------------
build_options
get_modify_options(const cmdline_args& args, modify_options opt)
{
    error_messages err;

    auto cli = modify_mode_cli(opt, err);

    auto result = clipp::parse(args, cli);

    if(!result || err.any()) {
        raise_default_error(err, "modify", modify_mode_usage());
    }

    // use settings from database file as defaults
    auto db = make_database(opt.dbfile, database::scope::metadata_only);

    const auto& ts = db.target_sketcher();
    auto& sk = opt.sketching;
    sk.kmerlen   = ts.kmer_size();
    sk.sketchlen = ts.sketch_size();
    sk.winlen    = ts.window_size();
    sk.winstride = ts.window_stride();

    opt.dbconfig.maxLoadFactor = db.max_load_factor();
    opt.dbconfig.maxLocationsPerFeature = db.max_locations_per_feature();

    // parse again
    clipp::parse(args, cli);

    augment_taxonomy_options(opt.taxonomy);
    replace_directories_with_contained_files(opt.infiles);

    if(opt.dbconfig.maxLocationsPerFeature < 0)
        opt.dbconfig.maxLocationsPerFeature = database::max_supported_locations_per_feature();

    return opt;
}



//-------------------------------------------------------------------
string modify_mode_usage() {
    return
    "    metacache modify <database> <sequence file/directory>... [OPTION]...\n\n"
    "    metacache modify <database> [OPTION]... <sequence file/directory>...";
}



//-------------------------------------------------------------------
string modify_mode_examples() {
    return
    "    Add reference sequence 'penicillium.fna' to database 'fungi'\n"
    "        metacache modify fungi penicillium.fna\n"
    "\n"
    "    Add taxonomic information from NCBI to database 'myBacteria'\n"
    "        download_ncbi_taxonomy myTaxo\n"
    "        metacache modify myBacteria -taxonomy myTaxo\n";
}



//-------------------------------------------------------------------
string modify_mode_docs() {

    build_options opt;
    error_messages err;

    auto cli = modify_mode_cli(opt, err);

    string docs = "SYNOPSIS\n\n";

    docs += modify_mode_usage();

    docs += "\n\n\n"
        "DESCRIPTION\n"
        "\n"
        "    Add reference sequence and/or taxonomic information to an existing database.\n"
        "\n\n";

    docs += clipp::documentation(cli, cli_doc_formatting()).str();

    docs += "\n\n\nEXAMPLES\n";
    docs += modify_mode_examples();

    return docs;
}





/*************************************************************************//**
 *
 *
 *  Q U E R Y   M O D E
 *
 *
 *****************************************************************************/
/// @brief command line interface for classification parameter tuning
clipp::group
classification_params_cli(classification_options& opt, error_messages& err)
{
    using namespace clipp;

    return (
    (   option("-lowest") &
        value("rank", [&](const string& name) {
                auto r = rank_from_name(name, err);
                if(opt.lowestRank < taxon_rank::root) opt.lowestRank = r;
            })
            .if_missing([&]{ err += "Taxonomic rank missing after '-lowest'!"; })
    )
        %("Do not classify on ranks below <rank>\n"
          "(Valid values: "s + taxon_rank_names() + ")\n"s +
          "default: "s + taxonomy::rank_name(opt.lowestRank))
    ,
    (   option("-highest") &
        value("rank", [&](const string& name) {
                auto r = rank_from_name(name, err);
                if(opt.highestRank <= taxon_rank::root) opt.highestRank = r;
            })
            .if_missing([&]{ err += "Taxonomic rank missing after '-highest'!"; })
    )
        %("Do not classify on ranks above <rank>\n"
          "(Valid values: "s + taxon_rank_names() + ")\n"s +
          "default: "s + taxonomy::rank_name(opt.highestRank))
    ,
    (   option("-hitmin", "-hit-min", "-hits-min", "-hitsmin") &
        integer("t", opt.hitsMin)
            .if_missing([&]{ err += "Number missing after '-hitmin'!"; })
    )
        %("Sets classification threshhold to <t>.\n"
          "A read will not be classified if less than t features "
          "from the database match. Higher values will increase "
          "precision at the expense of sensitivity.\n"
          "default: "s + to_string(opt.hitsMin))
    ,
    (   option("-hitdiff", "-hit-diff", "-hitsdiff", "-hits-diff") &
        number("t", opt.hitsDiffFraction)
            .if_missing([&]{ err += "Number missing after '-hitdiff'!"; })
    )
        %("Sets classification threshhold to <t>.\n"
          "A read will not be classified if less than t features "
          "from the database match. Higher values will increase "
          "precision at the expense of sensitivity.\n"
          "default: "s + to_string(opt.hitsMin))
    ,
    (   option("-maxcand", "-max-cand") &
        integer("#", opt.maxNumCandidatesPerQuery)
            .if_missing([&]{ err += "Number missing after '-maxcand'!"; })
    )
        %("maximum number of reference taxon candidates to "
          "consider for each query;\n"
          "A large value can significantly decrease the querying speed!.\n"
          "default: "s + to_string(opt.maxNumCandidatesPerQuery))
    ,
    (   option("-cov-percentile") &
        number("p", opt.covPercentile)
            .if_missing([&]{ err += "Number missing after '-cov-percentile'!"; })
    )
        %("Remove the p-th percentile of hit reference sequences "
          "with the lowest coverage. Classification is done using "
          "only the remaining reference sequences. "
          "This can help to reduce false positives, especially when"
          "your input data has a high sequencing coverage.\n"
          "This feature decreases the querying speed!\n"
          "default: "s + (opt.covPercentile > 1e-3 ? "on" : "off")
    )
    );
}



//-------------------------------------------------------------------
/// @brief build mode command-line options
clipp::group
classification_output_format_cli(classification_output_formatting& opt,
                                 error_messages& err)
{
    using namespace clipp;
    return (
        one_of(
            option("-no-map", "-nomap").set(opt.mapViewMode, map_view_mode::none)
            %("Don't report classification for each individual query "
              "sequence; show summaries only (useful for quick tests).\n"
              "default: "s + (opt.mapViewMode == map_view_mode::none ? "on" : "off"))
            ,
            option("-mapped-only", "-mappedonly").set(opt.mapViewMode, map_view_mode::mapped_only)
            %("Don't list unclassified reads/read pairs.\n"
              "default: "s + (opt.mapViewMode == map_view_mode::mapped_only ? "on" : "off"))
        )
        ,
        option("-taxids", "-taxid").set(opt.taxonStyle.showId)
            %("Print taxon ids in addition to taxon names.\n"
              "default: "s + (opt.taxonStyle.showId ? "on" : "off"))
        ,
        option("-taxids-only", "-taxidsonly")
            .set(opt.taxonStyle.showId).set(opt.taxonStyle.showName,false)
            %("Print taxon ids instead of taxon names.\n"
              "default: "s + (opt.taxonStyle.showId && !opt.taxonStyle.showName ? "on" : "off"))
        ,
        option("-omit-ranks", "-omitranks").set(opt.taxonStyle.showRankName,false)
            %("Do not print taxon rank names.\n"
              "default: "s + (!opt.taxonStyle.showRankName ? "on" : "off"))
        ,
        option("-separate-cols", "-separatecols").set(opt.useSeparateCols)
            %("Prints *all* mapping information (rank, taxon name, taxon ids) "
              "in separate columns (see option '-separator').\n"
              "default: "s + (opt.useSeparateCols ? "on" : "off"))
        ,
        (   option("-separator") &
            value("text", [&](const string& arg) {
                    opt.tokens.column = sanitize_special_chars(arg);
                })
                .if_missing([&]{ err += "Text missing after '-separator'!"; })
        )
            % "Sets string that separates output columns.\n"
              "default: '\\t|\\t'"
        ,
        (   option("-comment") &
            value("text", opt.tokens.comment)
                .if_missing([&]{ err += "Text missing after '-comment'!"; })
        )
            %("Sets string that precedes comment (non-mapping) lines.\n"
              "default: '"s + opt.tokens.comment + "'")
        ,
        option("-queryids", "-query-ids").set(opt.showQueryIds)
            %("Show a unique id for each query.\n"
              "Note that in paired-end mode a query is a pair of two "
              "read sequences. This option will always be activated if "
              "option '-hits-per-seq' is given.\n"
              "default: "s + (opt.showQueryIds ? "on" : "off"))
        ,
        option("-lineage", "-lineages").set(opt.showLineage)
            %("Report complete lineage for per-read classification "
              "starting with the lowest rank found/allowed and "
              "ending with the highest rank allowed. See also "
              "options '-lowest' and '-highest'.\n"
              "default: "s + (opt.showLineage ? "on" : "off"))
    );
}



//-------------------------------------------------------------------
/// @brief build mode command-line options
clipp::group
classification_analysis_cli(classification_analysis_options& opt, error_messages& err)
{
    using namespace clipp;
    return (
        "ANALYSIS: ABUNDANCES" %
        (
            (option("-abundances", "-abundance").set(opt.showTaxAbundances) &
             opt_value("file", opt.abundanceFile))
                %("Show absolute and relative abundance of each taxon.\n"
                  "If a valid filename is given, the list will be written to "
                  "this file.\n"
                  "default: "s + (opt.showTaxAbundances ? "on" : "off"))
            ,
            (   option("-abundance-per") &
                value("rank", [&](const string& name) {
                        auto r = rank_from_name(name, err);
                        if(r < taxon_rank::root) opt.showAbundanceEstimatesOnRank = r;
                    })
                    .if_missing([&]{ err += "Taxonomic rank missing after '-abundance-per'!"; })
            )
                %("Show absolute and relative abundances for each "
                  "taxon on one specific rank.\n"
                  "Classifications on higher ranks will be estimated by "
                  "distributing them down according to the relative "
                  "abundances of classifications on or below the given rank. "
                  "(Valid values: "s + taxon_rank_names() + ")\n"s +
                  "If '-abundances <file>' was given, this list will "
                  "be printed to the same file.\n"
                  "default: "s + (false ? "on" : "off"))
        )
        ,
        "ANALYSIS: RAW DATABASE HITS" %
        (
            option("-tophits", "-top-hits").set(opt.showTopHits)
                %("For each query, print top feature hits in database.\n"
                  "default: "s + (opt.showTopHits ? "on" : "off"))
            ,
            option("-allhits", "-all-hits").set(opt.showAllHits)
                %("For each query, print all feature hits in database.\n"
                  "default: "s + (opt.showAllHits ? "on" : "off"))
            ,
            option("-locations").set(opt.showLocations).set(opt.showTopHits)
                %("Show locations in candidate reference sequences.\n"
                  "Activates option '-tophits'.\n"
                  "default: "s + (opt.showLocations ? "on" : "off"))
            ,
            (   option("-hits-per-ref", "-hits-per-seq",
                       "-hits-per-tgt", "-hits-per-target")
                    .set(opt.showHitsPerTargetList) &
                opt_value("file", opt.targetMappingsFile)
            )
                %("Shows a list of all hits for each reference sequence.\n"
                  "If this condensed list is all you need, you should "
                  "deactive the per-read mapping output with '-no-map'.\n"
                  "If a valid filename is given after '-hits-per-seq', "
                  "the list will be written to a separate file.\n"
                  "Option '-queryids' will be activated and the lowest "
                  "classification rank will be set to 'sequence'.\n"
                  "default: "s + (opt.showHitsPerTargetList ? "on" : "off"))
        )
        ,
        "ANALYSIS: ALIGNMENTS" %
        group(
            option("-align", "-alignment").set(opt.showAlignment)
                %("Show semi-global alignment to best candidate reference sequence.\n"
                  "Original files of reference sequences must be available.\n"
                  "This feature decreases the querying speed!\n"
                  "default: "s + (opt.showAlignment ? "on" : "off"))
        )
    );
}



//-------------------------------------------------------------------
/// @brief build mode command-line options
clipp::group
classification_evaluation_cli(classification_evaluation_options& opt,
                              error_messages&)
{
    using namespace clipp;
    return (
        option("-ground-truth", "-groundtruth")
            .set(opt.determineGroundTruth).set(opt.showGroundTruth)
            %("Report correct query taxa if known.\n"
              "Queries need to have either a 'taxid|<number>' entry in "
              "their header or a sequence id that is also present in "
              "the database.\n"
              "This feature decreases the querying speed!\n"
              "default: "s + (opt.determineGroundTruth ? "on" : "off"))
        ,
        option("-precision").set(opt.precision).set(opt.determineGroundTruth)
            %("Report precision & sensitivity "
              "by comparing query taxa (ground truth) and mapped taxa.\n"
              "Queries need to have either a 'taxid|<number>' entry in "
              "their header or a sequence id that is also found in "
              "the database.\n"
              "This feature decreases the querying speed!\n"
              "default: "s + (opt.precision ? "on" : "off"))
        ,
        option("-taxon-coverage")
            .set(opt.taxonCoverage).set(opt.precision)
            .set(opt.determineGroundTruth)
            %("Report true/false positives and true/false negatives."
              "This option turns on '-precision', so ground truth data "
              "needs to be available.\n"
              "This feature decreases the querying speed!\n"
              "default: "s + (opt.taxonCoverage ? "on" : "off"))
    );
}



//-------------------------------------------------------------------
/// @brief build mode command-line options
clipp::group
performance_options_cli(performance_tuning_options& opt, error_messages& err)
{
    using namespace clipp;
    return (
    (   option("-threads") &
        integer("#", opt.numThreads)
            .if_missing([&]{ err += "Number missing after '-threads'!"; })
    )
        %("Sets the maximum number of parallel threads to use."
          "default (on this machine): "s + to_string(opt.numThreads))
    ,
    (   option("-batch-size", "-batchsize") &
        integer("#", opt.batchSize)
            .if_missing([&]{ err += "Number missing after '-batch-size'!"; })
    )
        %("Process <#> many queries (reads or read pairs) per thread at once.\n"
          "default (on this machine): "s + to_string(opt.batchSize))
    ,
    (   option("-query-limit", "-querylimit") &
        integer("#", opt.queryLimit)
            .if_missing([&]{ err += "Number missing after '-query-limit'!"; })
    )
        %("Classify at max. <#> queries (reads or read pairs) per input file.\n"
          "default: "s + (opt.queryLimit < 1 ? "none"s : to_string(opt.queryLimit))
    )
    );
}



//-------------------------------------------------------------------
/// @brief build mode command-line options
clipp::group
query_mode_cli(query_options& opt, error_messages& err)
{
    using namespace clipp;

    return (
    "BASIC PARAMETERS" %
    (
        database_parameter(opt.dbfile, err)
        ,
        opt_values(match::prefix_not{"-"}, "sequence file/directory", opt.infiles)
            % "FASTA or FASTQ files containing genomic sequences "
              "(short reads, long reads, contigs, complete genomes, ...) "
              "that shall be classified.\n"
              "* If directory names are given, they will be searched for "
              "sequence files (at most 10 levels deep).\n"
              "* If no input filenames or directories are given, MetaCache will "
              "run in interactive query mode. This can be used to load the database into "
              "memory only once and then query it multiple times with different "
              "query options. "
    ),
    "MAPPING RESULTS OUTPUT" %
    one_of(
        (   option("-out") &
            value("file", opt.queryMappingsFile)
                .if_missing([&]{ err += "Output filename missing after '-out'!"; })
        )
            % "Redirect output to file <file>.\n"
              "If not specified, output will be written to stdout. "
              "If more than one input file was given all output "
              "will be concatenated into one file."
        ,
        (   option("-split-out", "-splitout").set(opt.splitOutputPerInput) &
            value("file", opt.queryMappingsFile)
                .if_missing([&]{ err += "Output filename missing after '-split-out'!"; })
        )
            % "Generate output and statistics for each input file "
              "separately. For each input file <in> an output file "
              "with name <file>_<in> will be written."
    ),
    "PAIRED-END READ HANDLING" %
    (   one_of(
            option("-pairfiles", "-pair-files", "-paired-files")
            .set(opt.pairing, pairing_mode::files)
            % "Interleave paired-end reads from two consecutive files, "
              "so that the nth read from file m and the nth read "
              "from file m+1 will be treated as a pair. "
              "If more than two files are provided, their names "
              "will be sorted before processing. Thus, the order "
              "defined by the filenames determines the pairing not "
              "the order in which they were given in the command line."
            ,
            option("-pairseq", "-pair-seq", "-paired-seq")
            .set(opt.pairing, pairing_mode::sequences)
            % "Two consecutive sequences (1+2, 3+4, ...) from each file "
              "will be treated as paired-end reads."
        ),

        (   option("-insertsize", "-insert-size") &
            integer("#", opt.classify.insertSizeMax)
                .if_missing([&]{ err += "Number missing after '-insertsize'!"; })
        )
            % "Maximum insert size to consider.\n"
              "default: sum of lengths of the individual reads"
    )
    ,
    "CLASSIFICATION" %
        classification_params_cli(opt.classify, err)
    ,
    "GENERAL OUTPUT FORMATTING" % (
        option("-no-summary", "-nosummary").set(opt.output.showSummary,false)
            %("Dont't show result summary & mapping statistics at the "
              "end of the mapping output\n"
              "default: "s + (!opt.output.showSummary ? "on" : "off"))
        ,
        option("-no-query-params", "-no-queryparams", "-noqueryparams")
            .set(opt.output.showQueryParams,false)
            %("Don't show query settings at the beginning of the "
              "mapping output\n"
              "default: "s + (!opt.output.showQueryParams ? "on" : "off"))
        ,
        option("-no-err", "-noerr", "-no-errors").set(opt.output.showErrors,false)
            %("Suppress all error messages.\n"
              "default: "s + (!opt.output.showErrors ? "on" : "off"))
    )
    ,
    "CLASSIFICATION RESULT FORMATTING" %
        classification_output_format_cli(opt.output.format, err)
    ,
    classification_analysis_cli(opt.output.analysis, err)
    ,
    "ADVANCED: GROUND TRUTH BASED EVALUATION" %
        classification_evaluation_cli(opt.output.evaluate, err)
    ,
    "ADVANCED: CUSTOM QUERY SKETCHING (SUBSAMPLING)" %
        sketching_options_cli(opt.sketching, err)
    ,
    "ADVANCED: DATABASE MODIFICATION" %
        database_storage_options_cli(opt.dbconfig, err)
    ,
    "ADVANCED: PERFORMANCE TUNING / TESTING" %
        performance_options_cli(opt.performance, err)
    ,
    catch_unknown(err)
    );
}



//-------------------------------------------------------------------
query_options
get_query_options(const cmdline_args& args, query_options opt)
{
    error_messages err;

    auto cli = query_mode_cli(opt, err);

    auto result = clipp::parse(args, cli);

    if(!result || err.any()) {
        raise_default_error(err, "query", query_mode_usage());
    }

    replace_directories_with_contained_files(opt.infiles);

    if(opt.pairing == pairing_mode::files) {
        if(opt.infiles.size() > 1) {
            std::sort(opt.infiles.begin(), opt.infiles.end());
        } else {
            // TODO warning that pairing_mode::files requires at least 2 files
            opt.pairing = pairing_mode::none;
        }
    }

    // interprest numbers > 1 as percentage
    auto& cl = opt.classify;
    if(cl.hitsDiffFraction > 1) cl.hitsDiffFraction *= 0.01;
    if(cl.covPercentile    > 1) cl.covPercentile    *= 0.01;

    if(cl.maxNumCandidatesPerQuery < 1) {
        cl.maxNumCandidatesPerQuery = std::numeric_limits<size_t>::max();
    }


    // classification rank consistency checks
    if(cl.lowestRank  > cl.highestRank) cl.lowestRank  = cl.highestRank;
    if(cl.highestRank < cl.lowestRank)  cl.highestRank = cl.lowestRank;


    // processing option checks
    auto& perf = opt.performance;
    if(perf.numThreads < 1) perf.numThreads = 1;
    if(perf.batchSize  < 1) perf.batchSize  = 1;
    if(perf.queryLimit < 0) perf.queryLimit = 0;


    //output file consistency checks
    auto& ana = opt.output.analysis;
    if(ana.targetMappingsFile == opt.queryMappingsFile) ana.targetMappingsFile.clear();
    if(ana.abundanceFile == opt.queryMappingsFile) ana.abundanceFile.clear();

    // output option checks and consistency

    //always show query ids if hits per target list requested
    auto& fmt = opt.output.format;
    // output ranks are the same as classification ranks
    fmt.lowestRank = cl.lowestRank;
    fmt.highestRank = cl.highestRank;

    if(ana.showHitsPerTargetList) fmt.showQueryIds = true;

    // modify output tokens for separate column printig
    if(fmt.useSeparateCols) {
        fmt.collapseUnclassifiedLineages = false;
        fmt.tokens.taxSeparator = fmt.tokens.column;
        fmt.tokens.rankSuffix   = fmt.tokens.column;
        fmt.tokens.taxidPrefix  = fmt.tokens.column;
        fmt.tokens.taxidSuffix  = "";
    }

    // showing hits changes the mapping mode!
    if(fmt.mapViewMode == map_view_mode::none && ana.showTopHits) {
        fmt.mapViewMode = map_view_mode::mapped_only;
    }
    else if(ana.showAllHits) {
        fmt.mapViewMode = map_view_mode::all;
    }

    return opt;
}



//-------------------------------------------------------------------
string query_mode_usage() {
    return
    "    metacache query <database>\n\n"
    "    metacache query <database> <sequence file/directory>... [OPTION]...\n\n"
    "    metacache query <database> [OPTION]... <sequence file/directory>...";
}



//-------------------------------------------------------------------
string query_mode_examples() {
    return
    "    Query all sequences in 'myreads.fna' against pre-built database 'refseq':\n"
    "        metacache query refseq myreads.fna -out results.txt\n"
    "\n"
    "    Query all sequences in multiple files against database 'refseq':\n"
    "        metacache query refseq reads1.fna reads2.fna reads3.fna\n"
    "\n"
    "    Query all sequence files in folder 'test' againgst database 'refseq':\n"
    "        metacache query refseq test\n"
    "\n"
    "    Query multiple files and folder contents against database 'refseq':\n"
    "        metacache query refseq file1.fna folder1 file2.fna file3.fna folder2\n"
    "\n"
    "    Perform a precision test and show all ranks for each classification result:\n"
    "        metacache query refseq reads.fna -precision -allranks -out results.txt\n"
    "\n"
    "    Load database in interactive query mode, then query multiple read batches\n"
    "        metacache query refseq\n"
    "        reads1.fa reads2.fa -pairfiles -insertsize 400\n"
    "        reads3.fa -pairseq -insertsize 300\n"
    "        reads4.fa -lineage\n";
}



//-------------------------------------------------------------------
string query_mode_docs() {

    query_options opt;
    error_messages err;
    auto cli = query_mode_cli(opt, err);

    string docs = "SYNOPSIS\n\n";

    docs += query_mode_usage();

    docs += "\n\n\n"
        "DESCRIPTION\n"
        "\n"
        "    Map sequences (short reads, long reads, genome fragments, ...)\n"
        "    to their most likely taxon of origin.\n"
        "\n\n";

    docs += clipp::documentation(cli, cli_doc_formatting()).str();

    docs += "\n\n\nEXAMPLES\n\n";
    docs += query_mode_examples();

    docs += "\n\n"
        "OUTPUT FORMAT\n"
        "\n"
        "    MetaCache's default read mapping output format is:\n"
        "    read_header | rank:taxon_name\n"
        "\n"
        "    This will not be changed in the future to avoid breaking anyone's\n"
        "    pipelines. Command line options won't change in the near future for the\n"
        "    same reason. The following table shows some of the possible mapping\n"
        "    layouts with their associated command line arguments:\n"
        "\n"
        "    read mapping layout                          command line arguments\n"
        "    ---------------------------------------      ---------------------------------\n"
        "    read_header | taxon_id                       -taxids-only -omit-ranks\n"
        "    read_header | taxon_name                     -omit-ranks\n"
        "    read_header | taxon_name(taxon_id)           -taxids -omit-ranks\n"
        "    read_header | taxon_name | taxon_id          -taxids -omit-ranks -separate-cols\n"
        "    read_header | rank:taxon_id                  -taxids-only\n"
        "    read_header | rank:taxon_name\n"
        "    read_header | rank:taxon_name(taxon_id)      -taxids\n"
        "    read_header | rank | taxon_id                -taxids-only -separate-cols\n"
        "    read_header | rank | taxon_name              -separate-cols\n"
        "    read_header | rank | taxon_name | taxon_id   -taxids -separate-cols\n"
        "\n"
        "    Note that the separator '\\t|\\t' can be changed to something else with\n"
        "    the command line option '-separator <text>'.\n"
        "\n"
        "    Note that the default lowest taxon rank is 'sequence'. Sequence-level taxon\n"
        "    ids have negative numbers in order to not interfere with NCBI taxon ids.\n"
        "    Each reference sequence is added as its own taxon below the\n"
        "    lowest known NCBI taxon for that sequence. If you do not want to classify\n"
        "    at sequence-level, you can set a higher rank as lowest classification rank\n"
        "    with the '-lowest' command line option: '-lowest species' or\n"
        "    '-lowest subspecies' or '-lowest genus', etc.\n";

    return docs;
}





/*************************************************************************//**
 *
 *
 *  M E R G E   M O D E
 *
 *
 *****************************************************************************/

/// @brief build mode command-line options
clipp::group
merge_mode_cli(merge_options& opt, error_messages& err)
{
    using namespace clipp;

    auto& qry = opt.query;

    return (
    "REQUIRED PARAMETERS" %
    (
        values(match::prefix_not{"-"}, "result file/directory", opt.infiles)
            .if_missing([&]{ err += "No result filenames provided!"; })
            % "MetaCache result files.\n"
              "If directory names are given, they will be searched for "
              "sequence files (at most 10 levels deep).\n"
              "IMPORTANT: Result files must have been produced with:\n"
              "    -tophits -queryids -lowest species\n"
              "and must NOT be run with options that suppress or alter the "
              "default output like, e.g.: -no-map, -no-summary, -separator, etc.\n"
        ,
        (   required("-taxonomy")
                .if_missing([&]{ err += "Taxonomy path missing. Use '-taxonomy <path>'!"; })
            &
            value("path", opt.taxonomy.path)
                .if_missing([&]{ err += "Taxonomy path missing after '-taxonomy'!"; })
        )
            % "directory with taxonomic hierarchy data (see NCBI's taxonomic data files)"
    )
    ,
    "CLASSIFICATION" %
        classification_params_cli(qry.classify, err)
    ,
    "GENERAL OUTPUT" % (
        info_level_cli(opt.infoLevel, err)
        ,
        option("-no-summary", "-nosummary").set(qry.output.showSummary,false)
            %("Dont't show result summary & mapping statistics at the "
              "end of the mapping output\n"
              "default: "s + (!qry.output.showSummary ? "on" : "off"))
        ,
        option("-no-query-params", "-no-queryparams", "-noqueryparams")
            .set(qry.output.showQueryParams,false)
            %("Don't show query settings at the beginning of the "
              "mapping output\n"
              "default: "s + (!qry.output.showQueryParams ? "on" : "off"))
        ,
        option("-no-err", "-noerr", "-no-errors").set(qry.output.showErrors,false)
            %("Suppress all error messages.\n"
              "default: "s + (!qry.output.showErrors ? "on" : "off"))
    )
    ,
    "CLASSIFICATION RESULT FORMATTING" %
        classification_output_format_cli(qry.output.format, err)
    ,
    "ANALYSIS" %
        classification_analysis_cli(qry.output.analysis, err)
    ,
    "ADVANCED: GROUND TRUTH BASED EVALUATION" %
        classification_evaluation_cli(qry.output.evaluate, err)
    ,
    "ADVANCED: CUSTOM QUERY SKETCHING (SUBSAMPLING)" %
        sketching_options_cli(qry.sketching, err)
    ,
    "ADVANCED: DATABASE MODIFICATION" %
        database_storage_options_cli(qry.dbconfig, err)
    ,
    "ADVANCED: PERFORMANCE TUNING / TESTING" %
        performance_options_cli(qry.performance, err)
    ,
    catch_unknown(err)
    );
}



//-------------------------------------------------------------------
merge_options
get_merge_options(const cmdline_args& args, merge_options opt)
{
    error_messages err;

    auto cli = merge_mode_cli(opt, err);

    auto result = clipp::parse(args, cli);

    if(!result || err.any()) {
        raise_default_error(err, "merge", merge_mode_usage());
    }

    replace_directories_with_contained_files(opt.infiles);
    std::sort(opt.infiles.begin(), opt.infiles.end());

    auto& qo = opt.query;
    qo = get_query_options(args, {});

    if(qo.classify.hitsMin == 0) {
        qo.classify.hitsMin = 5;
    }
    if(qo.classify.lowestRank < taxon_rank::Species) {
        qo.classify.lowestRank = taxon_rank::Species;
    }
    if(qo.output.format.lowestRank < taxon_rank::Species) {
        qo.output.format.lowestRank = taxon_rank::Species;
    }
    if(qo.performance.numThreads > 1) {
        qo.performance.numThreads = 1;
    }

    return opt;
}



//-------------------------------------------------------------------
string merge_mode_usage() {
    return
    "    metacache merge <result file/directory>... -taxonomy <path> [OPTION]...\n\n"
    "    metacache merge -taxonomy <path> [OPTION]... <result file/directory>...";
}



//-------------------------------------------------------------------
string merge_mode_examples() {
    return "";
}



//-------------------------------------------------------------------
string merge_mode_docs() {

    merge_options opt;
    error_messages err;
    auto cli = merge_mode_cli(opt, err);

    string docs = "SYNOPSIS\n\n";

    docs += merge_mode_usage();

    docs += "\n\n\n"
        "DESCRIPTION\n"
        "\n"
        "    This mode classifies reads by merging the results of multiple, independent\n"
        "    queries. These might have been obtained by querying one database with\n"
        "    different parameters or by querying different databases with different\n"
        "    reference sequences or build options.\n"
        "\n"
        "    IMPORTANT: In order to be mergable, independent queries\n"
        "    need to be run with options:\n"
        "        -tophits -queryids -lowest species\n"
        "    and must NOT be run with options that suppress or alter default output\n"
        "    like, e.g.: -no-map, -no-summary, -separator, etc.\n"
        "\n"
        "    Possible Use Case:\n"
        "    If your system has not enough memory for one large database, you can\n"
        "    split up the set of reference genomes into several databases and query these\n"
        "    in succession. The results of these independent query runs can then be\n"
        "    merged to obtain a classification based on the whole set of genomes.\n"
        "\n\n";

    docs += clipp::documentation(cli, cli_doc_formatting()).str();

    return docs;
}




/*************************************************************************//**
 *
 *
 *  I N F O   M O D E
 *
 *
 *****************************************************************************/
/// @brief shared command-line options for taxonomy
clipp::group
info_mode_cli(info_options& opt, error_messages& err)
{
    using namespace clipp;

    return "PARAMETERS" % (
        database_parameter(opt.dbfile, err)
            .required(false).set(opt.mode, info_mode::db_config)
        ,
        one_of(
            command(""), // dummy
            (
                command("reference", "references", "ref",
                        "target", "targets", "tgt",
                        "sequence", "sequences", "seq")
                    .set(opt.mode, info_mode::targets),
                opt_values("sequence_id", opt.targetIds)
            ),
            (
                command("rank").set(opt.mode, info_mode::tax_ranks),
                value("rank_name", [&](const string& name) {
                        opt.rank = rank_from_name(name, err);
                    })
                    .if_missing([&]{ err += "Rank name missing!"; })
                    .doc("Valid values: "s + taxon_rank_names())
            ),
            command("lineages", "lineage", "lin")
                .set(opt.mode, info_mode::tax_lineages)
            ,
            command("statistics", "stat")
                .set(opt.mode, info_mode::db_statistics)
            ,
            command("locations", "loc", "featuremap", "features")
                .set(opt.mode, info_mode::db_feature_map)
            ,
            command("featurecounts")
                .set(opt.mode, info_mode::db_feature_counts)
        )
        ,
        catch_unknown(err)
    );

}



//-------------------------------------------------------------------
info_options
get_info_options(const cmdline_args& args)
{
    info_options opt;
    error_messages err;

    auto cli = info_mode_cli(opt, err);

    auto result = clipp::parse(args, cli);

    if(!result || err.any()) {
        raise_default_error(err, "info", info_mode_usage());
    }

    return opt;
}



//-------------------------------------------------------------------
string info_mode_usage()
{
    info_options opt;
    error_messages err;
    const auto cli = info_mode_cli(opt, err);
    return clipp::usage_lines(cli, "metacache info", cli_usage_formatting()).str();
}



//-------------------------------------------------------------------
string info_mode_examples() {
    return
    "    List metadata for all reference sequences in database 'refseq':\n"
    "        metacache info refseq.db ref\n"
    "\n"
    "    List metadata for the sequence with id NC_12345.6 in database 'refseq':\n"
    "        metacache info refseq.db ref NC_12345.6\n"
    "\n"
    "    List distribution of the number of sequences on rank 'phylum':\n"
    "        metacache info refseq.db rank phylum\n";
}



//-------------------------------------------------------------------
string info_mode_docs() {

    info_options opt;
    error_messages err;
    const auto cli = info_mode_cli(opt, err);

    string docs = "SYNOPSIS\n\n";

    docs += clipp::usage_lines(cli, "metacache info", cli_usage_formatting()).str();

    docs += "\n\n\n"
        "DESCRIPTION\n"
        "\n"
        "    Display (meta-)information stored in a database.\n"
        "\n\n"
        "SUB-MODES\n"
        "\n"
        "    metacache info\n"
        "        show basic properties of MetaCache executable (data type widths, etc.)\n"
        "\n"
        "    matacache info <database>\n"
        "        show basic properties of <database>\n"
        "\n"
        "    matacache info <database> ref[erence]\n"
        "       list meta information for all reference sequences in <database>\n"
        "\n"
        "    matacache info <database> ref[erence] <sequence_id>...\n"
        "       list meta information for specific reference sequences\n"
        "\n"
        "    matacache info <database> rank <rank_name>\n"
        "       list reference sequence distribution on rank <rank_name>\n"
        "\n"
        "    matacache info <database> lin[eages]\n"
        "       print table with ranked lineages for all reference sequences\n"
        "\n"
        "    matacache info <database> stat[istics]\n"
        "       print database statistics / hash table properties\n"
        "\n"
        "    matacache info <database> loc[ations]\n"
        "       print map (feature -> list of reference locations)\n"
        "\n"
        "    matacache info <database> featurecounts\n"
        "       print map (feature -> number of reference locations)\n"

        "\n\n";

    docs += clipp::documentation{cli, cli_doc_formatting()}.str();

    docs += "\n\n\nEXAMPLES\n\n";
    docs += info_mode_examples();

    return docs;
}



} // namespace mc

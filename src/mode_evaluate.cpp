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

#include "args_handling.h"
#include "timer.h"
#include "filesys_utility.h"
#include "cmdline_utility.h"
#include "sequence_io.h"
#include "sequence_view.h"
#include "alignment.h"

#include "modes_common.h"


namespace mc {

/*****************************************************************************
 *
 *
 *****************************************************************************/
using taxon_rank     = database::taxon_rank;
using taxon          = database::taxon;
using target_id      = database::target_id;
using taxon_id       = database::taxon_id;
using ranked_lineage = database::ranked_lineage;
using match_result   = database::match_result;




/*****************************************************************************
 *
 * @brief query options
 *
 *****************************************************************************/
struct eval_options
{
    //show all ranks that a sequence could be classified on
    bool showLineage = false;
    //what to show of a taxon
    taxon_print_mode showTaxaAs = taxon_print_mode::name_only;
    //prefix for each non-mapping line
    std::string comment = "# ";
    //separates individual mapping fields
    std::string outSeparator = "\t|\t";

    //-----------------------------------------------------
    bool testPrecision = false;
    bool testCoverage = false;

    //-------------------------------------------
    //filenames
    std::string dbfile;
    std::string groundTruthFile;
    std::vector<std::string> infiles;
    std::string outfile;
};




/*****************************************************************************
 *
 * @brief command line args -> query options
 *
 *****************************************************************************/
eval_options
get_eval_options(const args_parser& args)
{
    const eval_options defaults;

    eval_options param;

    //files
    param.dbfile = database_name(args);

    param.infiles = sequence_filenames(args);
    if(param.infiles.empty()) {
        throw std::invalid_argument{"No sequence filenames provided."};
    }

    param.testCoverage = args.contains("coverage");
    param.testPrecision = param.testCoverage || args.contains("precision");

    //output formatting
    param.showLineage = args.contains("lineage");

    param.outfile = args.get<std::string>("out", "");
    param.outSeparator = args.get<std::string>("separator", "\t|\t");

    if(args.contains({"taxidsonly","taxids-only","taxids_only",
                      "taxidonly", "taxid-only", "taxid_only"}))
    {
        param.showTaxaAs = taxon_print_mode::id_only;
    }
    else if(args.contains({"taxids", "taxid"})) {
        param.showTaxaAs = taxon_print_mode::id_name;
    }
    else {
        param.showTaxaAs = taxon_print_mode::name_only;
    }

    param.groundTruthFile = args.contains({"ground-truth", "ground_truth",
                                           "groundtruth"});

    return param;
}




/*****************************************************************************
 *
 * @brief
 *
 *****************************************************************************/
void main_mode_evaluate(const args_parser& )
{
//    std::cout << "Evaluating read mappings against ground truth." << std::endl;

//    auto param = get_eval_options(args);

//    auto db = make_database<database>(param.dbfile);

}


} // namespace mc

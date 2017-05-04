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

#include "classification.h"
#include "sequence_io.h"


namespace mc {


//-------------------------------------------------------------------
classification
ground_truth(const database& db, const std::string& header)
{

    //try to extract query id and find the corresponding target in database
    auto tid = db.target_id_of_sequence(extract_ncbi_accession_version_number(header));
    //not found yet
    if(!db.is_valid(tid)) {
        tid = db.target_id_of_sequence(extract_ncbi_accession_number(header));
        //not found yet
//        if(!db.is_valid(tid)) {
//            tid = db.target_id_of_sequence(extract_genbank_identifier(header));
//        }
    }

    auto taxid = extract_taxon_id(header);
    //make sure to get the lowest taxon that has a rank assigned
    if(taxid > 0) {
        for(auto id : db.ranks(db.taxon_with_id(taxid))) {
            if(id > 0) {
                taxid = id;
                break;
            }
        }
    }

    //if target known and in db
    if(db.is_valid(tid)) {
        //if taxid could not be determined solely from header, use the one
        //from the target in the db
        if(taxid < 1) {
            //use lowest taxon that has a rank assigned
            for(auto id : db.ranks_of_target(tid)) {
                if(id > 0) {
                    taxid = id;
                    break;
                }
            }
        }

        return classification { tid, taxid > 0 ? &db.taxon_with_id(taxid) : nullptr };
    }

    return classification { taxid > 0 ? &db.taxon_with_id(taxid) : nullptr };
}



//-------------------------------------------------------------------
taxon_rank
lowest_common_rank(const database& db,
                   const classification& a,
                   const classification& b)
{
    if(a.sequence_level() && b.sequence_level() &&
        a.target() == b.target()) return taxon_rank::Sequence;

    if(a.has_taxon() && b.has_taxon()) {
        return db.ranked_lca(a.tax(), b.tax()).rank();
    }

    return taxon_rank::none;
}


} // namespace mc

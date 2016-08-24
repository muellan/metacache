#ifndef MC_CONFIG_H_
#define MC_CONFIG_H_


#include "sketch_database.h"
#include "min_hasher.h"


namespace mc {

/*****************************************************************************
 *
 *
 *
 *****************************************************************************/
using sequence = std::string;

using sketcher = single_function_min_hasher; //default, for 0 <= k <= 32
// using sketcher = single_function_min_hasher64; //for 33 <= k <= 64
// using sketcher = multi_function_min_hasher; //different hash fun for each sketch value

using database   = sketch_database<sequence,sketcher>;
using taxon_rank = database::taxon_rank;
using genome_id  = database::genome_id;

using top_matches_in_contiguous_window_range
        = matches_in_contiguous_window_range_top<2>;
//        = matches_in_contiguous_window_range_top<8>;  //will use majority voting scheme


} // namespace mc


#endif

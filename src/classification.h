#ifndef MC_CLASSIFICATION_H_
#define MC_CLASSIFICATION_H_

#include "classification_statistics.h"
#include "config.h"


namespace mc {

void process_database_answer(
     const database& db, const query_options& opt,
     const std::string& header, const sequence& query1, const sequence& query2,
     matches_per_location&& hits,
     classification_statistics& stats, std::ostream& os);


} // namespace mc

#endif

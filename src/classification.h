#ifndef MC_CLASSIFICATION_H_
#define MC_CLASSIFICATION_H_

#include <vector>
#include <string>

#include "classification_statistics.h"
#include "config.h"


namespace mc {


/*************************************************************************//**
 *
 * @brief classification result target
 *
 *****************************************************************************/
struct classification_results
{
    explicit
    classification_results(std::ostream& outTarget = std::cout,
                           std::ostream& statusTarget = std::cerr)
    :
        out(outTarget), status(statusTarget)
    {}

    std::ostream& out;
    std::ostream& status;
    classification_statistics statistics;
};


//-------------------------------------------------------------------
void map_reads_to_targets(
    const std::vector<std::string>& inputFilenames,
    const database&, const query_options&,
    classification_results&);




} // namespace mc

#endif

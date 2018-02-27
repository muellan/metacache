#ifndef MC_QUERYING_H_
#define MC_QUERYING_H_

#include <vector>

#include "config.h"
#include "sequence_io.h"
#include "cmdline_utility.h"
#include "parallel_task_queue.h"
#include "query_options.h"


namespace mc {


/*****************************************************************************
 *
 * @brief query db with batches of reads from one sequence source
 *
 * @tparam ResultProcessor  function/function object taking query results
 *
 *****************************************************************************/
template<class ResultProcessor>
void batched_query(
    parallel_queue& queue,
    const database& db, const query_options& opt,
    sequence_reader& reader, std::ostream& os,
    ResultProcessor&& process)
{
    const auto load = 32 * queue.concurrency();
    const auto batchSize = opt.batchSize > 0 ? opt.batchSize
                                             : 4096 * queue.concurrency();
    std::mutex mtx;
    std::atomic<std::uint64_t> queryLimit(opt.queryLimit);

    while(reader.has_next() && queryLimit > 0) {
        if(queue.unsafe_waiting() < load) {
            queue.enqueue([&]{
                std::ostringstream obuf;

                for(std::size_t i = 0;
                    i < batchSize && queryLimit > 0 && reader.has_next();
                    ++i, --queryLimit)
                {
                    auto query = reader.next();
                    if(!query.header.empty()) {
                        process(query.header, query.data, sequence{},
                                db.matches(query.data), obuf);
                    }
                }
                //flush output buffer
                std::lock_guard<std::mutex> lock(mtx);
                os << obuf.str();
            }); //enqueue
        }
    }
    //wait for all enqueued tasks to finish
    queue.wait();
}



/*****************************************************************************
 *
 * @brief query db with batches of reads from two sequence sources
 *
 * @tparam ResultProcessor  function/function object taking query results
 *
 *****************************************************************************/
template<class ResultProcessor>
void batched_pair_query(
    parallel_queue& queue,
    const database& db, const query_options& opt,
    sequence_reader& reader1, sequence_reader& reader2,
    std::ostream& os,
    ResultProcessor&& process)

{
    const auto load = 32 * queue.concurrency();
    const auto batchSize = opt.batchSize > 0 ? opt.batchSize
                                             : 4096 * queue.concurrency();
    std::mutex mtx1;
    std::mutex mtx2;
    std::atomic<std::uint64_t> queryLimit(opt.queryLimit);

    while(reader1.has_next() && reader2.has_next() && queryLimit > 0) {
        if(queue.unsafe_waiting() < load) {
            queue.enqueue([&]{
                std::ostringstream obuf;

                for(std::size_t i = 0; i < batchSize && queryLimit > 0;
                    ++i, --queryLimit)
                {
                    //make sure that both queries are read in sequence
                    std::unique_lock<std::mutex> lock1(mtx1);
                    if(!reader1.has_next() || !reader2.has_next()) break;
                    auto query1 = reader1.next();
                    auto query2 = reader2.next();
                    lock1.unlock();

                    if(!query1.header.empty()) {
                        auto matches = db.matches(query1.data);
                        db.accumulate_matches(query2.data, matches);

                        process(query1.header, query1.data, query2.data,
                                std::move(matches), obuf);
                    }
                }
                //flush output buffer
                std::lock_guard<std::mutex> lock(mtx2);
                os << obuf.str();
            }); //enqueue
        }
    }
    //wait for all enqueued tasks to finish
    queue.wait();
}



/*****************************************************************************
 *
 * @brief classify sequences per input file
 *
 * @tparam ResultProcessor  function/function object taking query results
 *
 *****************************************************************************/
template<class ResultProcessor>
void query_per_file(
    parallel_queue& queue,
    const database& db, const query_options& opt,
    const std::vector<std::string>& infilenames,
    std::ostream& os, ResultProcessor&& process)
{
    for(size_t i = 0; i < infilenames.size(); ++i) {
        const auto& fname = infilenames[i];
        if(opt.mapViewMode != map_view_mode::none) {
            os << opt.comment << fname << '\n';
        }
        else if(!opt.outfile.empty() && infilenames.size() > 1) {
            show_progress_indicator(i / float(infilenames.size()));
        }
        try {
            auto reader = make_sequence_reader(fname);
            if(opt.pairing == pairing_mode::sequences) {
                batched_pair_query(queue, db, opt, *reader, *reader, os, process);
            } else {
                batched_query(queue, db, opt, *reader, os, process);
            }
        }
        catch(std::exception& e) {
            if(opt.mapViewMode != map_view_mode::none) {
                os << "FAIL: " << e.what() << std::endl;
            }
        }
    }
}



/*****************************************************************************
 *
 * @brief classifies paired-end sequences in pairs of files
 *
 * @tparam ResultProcessor  function/function object taking query results
 *
 *****************************************************************************/
template<class ResultProcessor>
void query_with_file_pairs(
    parallel_queue& queue,
    const database& db, const query_options& opt,
    const std::vector<std::string>& infilenames,
    std::ostream& os, ResultProcessor&& process)
{
    //pair up reads from two consecutive files in the list
    for(size_t i = 0; i < infilenames.size(); i += 2) {
        const auto& fname1 = infilenames[i];
        const auto& fname2 = infilenames[i+1];
        if(opt.mapViewMode != map_view_mode::none) {
            os << opt.comment << fname1 << " + " << fname2 << '\n';
        }
        else if(!opt.outfile.empty() && infilenames.size() > 1) {
            show_progress_indicator(i / float(infilenames.size()));
        }
        try {
            auto reader1 = make_sequence_reader(fname1);
            auto reader2 = make_sequence_reader(fname2);
            batched_pair_query(queue, db, opt, *reader1, *reader2, os, process);
        }
        catch(std::exception& e) {
            if(opt.mapViewMode != map_view_mode::none) {
                os << "FAIL: " << e.what() << std::endl;
            }
        }
    }
}


} // namespace mc


#endif

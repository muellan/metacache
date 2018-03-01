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
 * @brief queries database with batches of reads from one sequence source
 *        produces one match list per sequence
 *
 * @tparam BufferSource  returns a per-batch buffer object
 *
 * @tparam BufferUpdate  takes database matches of one query and a buffer
 *                       must be thread-safe (only const operations on DB!)
 *
 * @tparam BufferSink    recieves buffer after batch is finished
 *
 *****************************************************************************/
template<class BufferSource, class BufferUpdate, class BufferSink>
void query_batched(
    parallel_queue& queue, sequence_reader& reader,
    const database& db, const query_processing_options& opt,
    BufferSource&& getBuffer, BufferUpdate&& update, BufferSink&& finalize)
{
    const auto load = 32 * queue.concurrency();

    std::mutex mtx;
    std::atomic<std::uint64_t> queryLimit(opt.queryLimit);

    while(reader.has_next() && queryLimit > 0) {
        if(queue.unsafe_waiting() < load) {
            queue.enqueue([&] {
                auto buffer = getBuffer();

                for(std::size_t i = 0;
                    i < opt.batchSize && queryLimit > 0 && reader.has_next();
                    ++i, --queryLimit)
                {
                    auto query = reader.next();
                    if(!query.header.empty()) {
                        update(buffer, db.matches(query.data),
                               query.header, query.data, sequence{});
                    }
                }
                std::lock_guard<std::mutex> lock(mtx);
                finalize(std::move(buffer));
            }); //enqueue
        }
    }
    //wait for all enqueued tasks to finish
    queue.wait();
}



/*****************************************************************************
 *
 * @brief queries database with batches of reads from two sequence sources
 *        produces one merged match list per sequence pair
 *
 * @tparam BufferSource  returns a per-batch buffer object
 *
 * @tparam BufferUpdate  takes database matches of one query and a buffer
 *                       must be thread-safe (only const operations on DB!)
 *
 * @tparam BufferSink    recieves buffer after batch is finished
 *
 *****************************************************************************/
template<class BufferSource, class BufferUpdate, class BufferSink>
void query_batched_merged_pairs(
    parallel_queue& queue, sequence_reader& reader1, sequence_reader& reader2,
    const database& db, const query_processing_options& opt,
    BufferSource&& getBuffer, BufferUpdate&& update, BufferSink&& finalize)
{
    const auto load = 32 * queue.concurrency();

    std::mutex mtx1;
    std::mutex mtx2;
    std::atomic<std::uint64_t> queryLimit(opt.queryLimit);

    while(reader1.has_next() && reader2.has_next() && queryLimit > 0) {
        if(queue.unsafe_waiting() < load) {
            queue.enqueue([&]{
                auto buffer = getBuffer();

                for(std::size_t i = 0; i < opt.batchSize && queryLimit > 0;
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
                        //merge matches with hits of second sequence
                        db.accumulate_matches(query2.data, matches);

                        update(buffer, std::move(matches),
                               query1.header, query1.data, query2.data);
                    }
                }
                std::lock_guard<std::mutex> lock(mtx2);
                finalize(std::move(buffer));
            }); //enqueue
        }
    }
    //wait for all enqueued tasks to finish
    queue.wait();
}



/*****************************************************************************
 *
 * @brief queries database with sequences per input file
 *        produces one match list per sequence or sequence pair
 *
 * @tparam BufferSource  returns a per-batch buffer object
 *
 * @tparam BufferUpdate  takes database matches of one query and a buffer
 *                       must be thread-safe (only const operations on DB!)
 *
 * @tparam BufferSink    recieves buffer after batch is finished
 *
 *****************************************************************************/
template<
    class BufferSource, class BufferUpdate, class BufferSink,
    class StatusCallback
>
void query_per_file(
    parallel_queue& queue,
    const std::vector<std::string>& infilenames,
    const database& db,
    const query_processing_options& opt,
    BufferSource&& bufsrc, BufferUpdate&& bufupdate, BufferSink&& bufsink,
    StatusCallback&& displayStatus)
{
    for(size_t i = 0; i < infilenames.size(); ++i) {
        const auto& fname = infilenames[i];

        const auto progress = infilenames.size() > 1 ? i / float(infilenames.size()) : -1;
        displayStatus(fname, progress);

        try {
            auto reader = make_sequence_reader(fname);
            if(opt.pairing == pairing_mode::sequences) {
                query_batched_merged_pairs(queue, *reader, *reader, db, opt,
                                           std::forward<BufferSource>(bufsrc),
                                           std::forward<BufferUpdate>(bufupdate),
                                           std::forward<BufferSink>(bufsink));
            } else {
                query_batched(queue, *reader, db, opt,
                              std::forward<BufferSource>(bufsrc),
                              std::forward<BufferUpdate>(bufupdate),
                              std::forward<BufferSink>(bufsink));
            }
        }
        catch(std::exception& e) {
            displayStatus(std::string("FAIL: ") + e.what(), progress);
        }
    }
}



/*****************************************************************************
 *
 * @brief queries database with paired-end sequences in pairs of files
 *        produces one merged match list per sequence pair
 *
 * @tparam BufferSource  returns a per-batch buffer object
 *
 * @tparam BufferUpdate  takes database matches of one query and a buffer
 *                       must be thread-safe (only const operations on DB!)
 *
 * @tparam BufferSink    recieves buffer after batch is finished
 *
 *****************************************************************************/
template<
    class BufferSource, class BufferUpdate, class BufferSink,
    class StatusCallback
>
void query_with_file_pairs(
    parallel_queue& queue,
    const std::vector<std::string>& infilenames,
    const database& db,
    const query_processing_options& opt,
    BufferSource&& bufsrc, BufferUpdate&& bufupdate, BufferSink&& bufsink,
    StatusCallback&& displayStatus)
{
    for(size_t i = 0; i < infilenames.size(); i += 2) {
        //pair up reads from two consecutive files in the list
        const auto& fname1 = infilenames[i];
        const auto& fname2 = infilenames[i+1];

        const auto progress = infilenames.size() > 1 ? i / float(infilenames.size()) : -1;
        displayStatus(fname1 + " + " + fname2, progress);

        try {
            auto reader1 = make_sequence_reader(fname1);
            auto reader2 = make_sequence_reader(fname2);

            query_batched_merged_pairs(queue, *reader1, *reader2, db, opt,
                                       std::forward<BufferSource>(bufsrc),
                                       std::forward<BufferUpdate>(bufupdate),
                                       std::forward<BufferSink>(bufsink));
        }
        catch(std::exception& e) {
            displayStatus(std::string("FAIL: ") + e.what(), progress);
        }
    }
}



/*****************************************************************************
 *
 * @brief queries database
 *
 * @tparam BufferSource  returns a per-batch buffer object
 *
 * @tparam BufferUpdate  takes database matches of one query and a buffer
 *                       must be thread-safe (only const operations on DB!)
 *
 * @tparam BufferSink    recieves buffer after batch is finished
 *
 *****************************************************************************/
template<
    class BufferSource, class BufferUpdate, class BufferSink,
    class StatusCallback
>
void query_database(
    const std::vector<std::string>& infilenames,
    const database& db,
    const query_processing_options& opt,
    BufferSource&& bufsrc, BufferUpdate&& bufupdate, BufferSink&& bufsink,
    StatusCallback&& displayStatus)
{
    parallel_queue queue(opt.numThreads);

    if(opt.pairing == pairing_mode::files) {
        query_with_file_pairs(queue, infilenames, db, opt,
                                     std::forward<BufferSource>(bufsrc),
                                     std::forward<BufferUpdate>(bufupdate),
                                     std::forward<BufferSink>(bufsink),
                                     std::forward<StatusCallback>(displayStatus));
    }
    else {
        query_per_file(queue, infilenames, db, opt,
                       std::forward<BufferSource>(bufsrc),
                       std::forward<BufferUpdate>(bufupdate),
                       std::forward<BufferSink>(bufsink),
                       std::forward<StatusCallback>(displayStatus));
    }
}


} // namespace mc


#endif

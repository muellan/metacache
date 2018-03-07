/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2018 André Müller (muellan@uni-mainz.de)
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

#ifndef MC_QUERYING_H_
#define MC_QUERYING_H_

#include <vector>

#include "config.h"
#include "sequence_io.h"
#include "cmdline_utility.h"
#include "parallel_task_queue.h"
#include "query_options.h"


namespace mc {


/*************************************************************************//**
 *
 * @brief
 *
 *****************************************************************************/
struct sequence_query
{
    explicit
    sequence_query(query_id qid, std::string headerText,
                   sequence s1, sequence s2 = sequence{}) noexcept
    :
        id{qid}, header(std::move(headerText)),
        seq1(std::move(s1)), seq2(std::move(s2)),
        groundTruth{nullptr}
    {}

    bool empty() const noexcept { return header.empty() || seq1.empty(); }

    query_id id;
    std::string header;
    sequence seq1;
    sequence seq2;  //2nd part of paired-end read
    const taxon* groundTruth;
};



/*************************************************************************//**
 *
 * @brief queries database with batches of reads from one sequence source
 *        produces one match list per sequence
 *
 * @tparam BufferSource  returns a per-batch buffer object
 *
 * @tparam BufferUpdate  takes database matches of one query and a buffer;
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
    std::atomic<query_id> qid(0);

    while(reader.has_next() && queryLimit > 0) {
        if(queue.unsafe_waiting() < load) {
            queue.enqueue([&] {
                auto buffer = getBuffer();

                for(std::size_t i = 0;
                    i < opt.batchSize && queryLimit > 0 && reader.has_next();
                    ++i, --queryLimit)
                {
                    auto seq = reader.next();
                    if(!seq.header.empty()) {
                        update(buffer,
                               sequence_query{++qid,
                                              std::move(seq.header),
                                              std::move(seq.data)},
                               db.matches(seq.data) );
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



/*************************************************************************//**
 *
 * @brief queries database with batches of reads from two sequence sources
 *        produces one merged match list per sequence pair
 *
 * @tparam BufferSource  returns a per-batch buffer object
 *
 * @tparam BufferUpdate  takes database matches of one query and a buffer;
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
    std::atomic<query_id> qid(0);

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
                    auto seq1 = reader1.next();
                    auto seq2 = reader2.next();
                    lock1.unlock();

                    if(!seq1.header.empty()) {
                        auto matches = db.matches(seq1.data);
                        //merge matches with hits of second sequence
                        db.accumulate_matches(seq2.data, matches);

                        update(buffer,
                               sequence_query{++qid,
                                              std::move(seq1.header),
                                              std::move(seq1.data),
                                              std::move(seq2.data)},
                               std::move(matches) );
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



/*************************************************************************//**
 *
 * @brief queries database with sequences per input file
 *        produces one match list per sequence or sequence pair
 *
 * @tparam BufferSource  returns a per-batch buffer object
 *
 * @tparam BufferUpdate  takes database matches of one query and a buffer;
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
    StatusCallback&& showStatus)
{
    for(size_t i = 0; i < infilenames.size(); ++i) {
        const auto& fname = infilenames[i];

        const auto progress = infilenames.size() > 1 ? i / float(infilenames.size()) : -1;
        showStatus(fname, progress);

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
            showStatus(std::string("FAIL: ") + e.what(), progress);
        }
    }
}



/*************************************************************************//**
 *
 * @brief queries database with paired-end sequences in pairs of files
 *        produces one merged match list per sequence pair
 *
 * @tparam BufferSource  returns a per-batch buffer object
 *
 * @tparam BufferUpdate  takes database matches of one query and a buffer;
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
    StatusCallback&& showStatus)
{
    for(size_t i = 0; i < infilenames.size(); i += 2) {
        //pair up reads from two consecutive files in the list
        const auto& fname1 = infilenames[i];
        const auto& fname2 = infilenames[i+1];

        const auto progress = infilenames.size() > 1 ? i / float(infilenames.size()) : -1;
        showStatus(fname1 + " + " + fname2, progress);

        try {
            auto reader1 = make_sequence_reader(fname1);
            auto reader2 = make_sequence_reader(fname2);

            query_batched_merged_pairs(queue, *reader1, *reader2, db, opt,
                                       std::forward<BufferSource>(bufsrc),
                                       std::forward<BufferUpdate>(bufupdate),
                                       std::forward<BufferSink>(bufsink));
        }
        catch(std::exception& e) {
            showStatus(std::string("FAIL: ") + e.what(), progress);
        }
    }
}



/*************************************************************************//**
 *
 * @brief queries database
 *
 * @tparam BufferSource  returns a per-batch buffer object
 *
 * @tparam BufferUpdate  takes database matches of one query and a buffer;
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
    StatusCallback&& showStatus)
{
    parallel_queue queue(opt.numThreads);

    if(opt.pairing == pairing_mode::files) {
        query_with_file_pairs(queue, infilenames, db, opt,
                                     std::forward<BufferSource>(bufsrc),
                                     std::forward<BufferUpdate>(bufupdate),
                                     std::forward<BufferSink>(bufsink),
                                     std::forward<StatusCallback>(showStatus));
    }
    else {
        query_per_file(queue, infilenames, db, opt,
                       std::forward<BufferSource>(bufsrc),
                       std::forward<BufferUpdate>(bufupdate),
                       std::forward<BufferSink>(bufsink),
                       std::forward<StatusCallback>(showStatus));
    }
}


} // namespace mc


#endif

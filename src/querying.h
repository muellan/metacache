/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2019 André Müller (muellan@uni-mainz.de)
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
#include <iostream>
#include <future>
#include <chrono>

#include "config.h"
#include "sequence_io.h"
#include "cmdline_utility.h"
#include "query_options.h"
#include "cmdline_utility.h"

#include "../dep/queue/concurrentqueue.h"


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
 * @brief queries database with a batch of reads
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
    class BufferSource, class BufferUpdate, class BufferSink
>
void process_sequence_batch(
    std::vector<sequence_pair_reader::sequence_pair>& sequenceBatch,
    const database& db, std::mutex& finalizeMtx,
    BufferSource& getBuffer, BufferUpdate& update, BufferSink& finalize)
{
    auto batchBuffer = getBuffer();
    bool bufferEmpty = true;
    match_locations matches;
    match_target_locations targetMatches;

    for(auto& seq : sequenceBatch) {
        if(!seq.first.header.empty()) {
            bufferEmpty = false;
            targetMatches.clear();

            db.accumulate_matches(seq.first.data, targetMatches);
            db.accumulate_matches(seq.second.data, targetMatches);
            targetMatches.sort();

            matches.clear();
            for(auto& m : targetMatches)
                matches.emplace_back(db.taxon_of_target(m.tgt), m.win);

            update(batchBuffer,
                    sequence_query{seq.first.index,
                                   std::move(seq.first.header),
                                   std::move(seq.first.data),
                                   std::move(seq.second.data)},
                    matches);
        }
    }

    sequenceBatch.clear();

    if(!bufferEmpty) {
        std::lock_guard<std::mutex> lock(finalizeMtx);
        finalize(std::move(batchBuffer));
    }
}



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
template<
    class BufferSource, class BufferUpdate, class BufferSink,
    class LogCallback
>
query_id query_batched(
    const std::string& filename1, const std::string& filename2,
    const database& db, const query_processing_options& opt,
    const query_id& startId,
    BufferSource&& bufsrc, BufferUpdate&& bufupdate, BufferSink&& bufsink,
    LogCallback&& log)
{
    using sequence_batch = std::vector<sequence_pair_reader::sequence_pair>;

    std::mutex finalizeMtx;
    std::atomic_bool queryingDone{0};
    moodycamel::ConcurrentQueue<sequence_batch> queue(opt.numThreads);

    //spawn consumers
    std::vector<std::future<void>> threads;
    for(int i = 0; i < opt.numThreads-1; ++i) {
        threads.emplace_back(std::async(std::launch::async, [&] {
                sequence_batch sequenceBatch;

                while(!queryingDone.load() || queue.size_approx() > 0) {
                    if(queue.try_dequeue(sequenceBatch)) {
                        process_sequence_batch(
                            sequenceBatch, db, finalizeMtx,
                            bufsrc, bufupdate, bufsink);
                    }
                }
            }));
    }

    std::int_least64_t queryLimit = opt.queryLimit;
    query_id endId = startId;

    //read sequences from file
    try {
        sequence_pair_reader reader{filename1, filename2};
        reader.index_offset(startId);

        while(reader.has_next() && queryLimit > 0) {
            sequence_batch sequenceBatch;
            sequenceBatch.reserve(opt.batchSize);

            //fill batch with sequences
            for(std::size_t i = 0; i < opt.batchSize && queryLimit > 0; ++i, --queryLimit) {
                if(!reader.has_next()) break;
                sequenceBatch.emplace_back(reader.next());
            }

            //enqueue batch or process if single threaded
            if(opt.numThreads > 1) {
                while(!queue.try_enqueue(std::move(sequenceBatch))) {
                    std::this_thread::sleep_for(std::chrono::milliseconds{10});
                };
            }
            else {
                process_sequence_batch(
                    sequenceBatch, db, finalizeMtx,
                    bufsrc, bufupdate, bufsink);
            }
        }

        endId = reader.index();
    }
    catch(file_access_error& e) {
        log(std::string("FAIL: ") + e.what());
    }
    catch(std::exception& e) {
        log(std::string("FAIL: ") + e.what());
    }

    queryingDone.store(1);
    // wait for all threads to finish
    for(unsigned int threadId = 0; threadId < threads.size(); ++threadId) {
        if(threads[threadId].valid()) {
            threads[threadId].get();
        }
    }

    return endId;
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
    class InfoCallback, class ProgressHandler, class LogHandler
>
void query_database(
    const std::vector<std::string>& infilenames,
    const database& db,
    const query_processing_options& opt,
    BufferSource&& bufsrc, BufferUpdate&& bufupdate, BufferSink&& bufsink,
    InfoCallback&& showInfo, ProgressHandler&& showProgress, LogHandler&& log)
{
    const size_t stride = opt.pairing == pairing_mode::files ? 1 : 0;
    const std::string nofile;
    query_id readIdOffset = 0;

    for(size_t i = 0; i < infilenames.size(); i += stride+1) {
        //pair up reads from two consecutive files in the list
        const auto& fname1 = infilenames[i];

        const auto& fname2 = (opt.pairing == pairing_mode::none)
                             ? nofile : infilenames[i+stride];

        if(opt.pairing == pairing_mode::files) {
            showInfo(fname1 + " + " + fname2);
        } else {
            showInfo(fname1);
        }
        showProgress(infilenames.size() > 1 ? i/float(infilenames.size()) : -1);

        readIdOffset = query_batched(fname1, fname2, db, opt, readIdOffset,
                                     std::forward<BufferSource>(bufsrc),
                                     std::forward<BufferUpdate>(bufupdate),
                                     std::forward<BufferSink>(bufsink),
                                     std::forward<LogHandler>(log));
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
    class BufferSource, class BufferUpdate, class BufferSink, class InfoCallback
>
void query_database(
    const std::vector<std::string>& infilenames,
    const database& db,
    const query_processing_options& opt,
    BufferSource&& bufsrc, BufferUpdate&& bufupdate, BufferSink&& bufsink,
    InfoCallback&& showInfo)
{
    query_database(infilenames, db, opt,
                   std::forward<BufferSource>(bufsrc),
                   std::forward<BufferUpdate>(bufupdate),
                   std::forward<BufferSink>(bufsink),
                   std::forward<InfoCallback>(showInfo),
                   [] (float p) { show_progress_indicator(std::cerr, p); },
                   [] (const std::string& s) {std::cerr << s << '\n'; }
                   );
}


} // namespace mc


#endif

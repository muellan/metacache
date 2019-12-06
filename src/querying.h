/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2019 André Müller (muellan@uni-mainz.de)
 *                       & Robin Kobus  (rkobus@uni-mainz.de)
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
 * @brief single query = id + header + read(pair)
 *
 *****************************************************************************/
struct sequence_query
{
    sequence_query() = default;
    sequence_query(const sequence_query&) = default;
    sequence_query(sequence_query&&) = default;
    sequence_query& operator = (const sequence_query&) = default;
    sequence_query& operator = (sequence_query&&) = default;

    explicit
    sequence_query(query_id qid, std::string headerText,
                   sequence s1, sequence s2 = sequence{}) noexcept
    :
        id{qid}, header(std::move(headerText)),
        seq1(std::move(s1)), seq2(std::move(s2)),
        groundTruth{nullptr}
    {}

    bool empty() const noexcept { return header.empty() || seq1.empty(); }

    query_id id = 0;
    std::string header;
    sequence seq1;
    sequence seq2;  // 2nd part of paired-end read
    const taxon* groundTruth = nullptr;
};


/** @brief batch of sequence queries */
using sequence_batch = std::vector<sequence_query>;




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
    sequence_batch& batch, const database& db, std::mutex& finalizeMtx,
    BufferSource& getBuffer, BufferUpdate& update, BufferSink& finalize)
{
    auto resultsBuffer = getBuffer();
    match_locations matches;
    match_target_locations targetMatches;

    for(auto& seq : batch) {
        targetMatches.clear();

        db.accumulate_matches(seq.seq1, targetMatches);
        db.accumulate_matches(seq.seq2, targetMatches);
        targetMatches.sort();

        matches.clear();
        for(auto& m : targetMatches) {
            matches.emplace_back(db.taxon_of_target(m.tgt), m.win);
        }

        update(resultsBuffer, seq, matches);
    }

    std::lock_guard<std::mutex> lock(finalizeMtx);
    finalize(std::move(resultsBuffer));
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
    std::mutex finalizeMtx;
    std::atomic_bool queryingDone{0};

    using batch_queue = moodycamel::ConcurrentQueue<sequence_batch>;

    if(opt.queryLimit < 1) return startId;
    const auto queryLimit = size_t(opt.queryLimit > 0 ? opt.queryLimit : std::numeric_limits<size_t>::max());
    const auto queueSize = opt.numThreads > 1 ? opt.numThreads+ 4 : 0;

    // queue for batches that can be reused
    batch_queue storageQueue(queueSize);
    for(auto i = 0; i < queueSize; ++i) {
        sequence_batch batch;
        batch.resize(opt.batchSize);
        storageQueue.enqueue(std::move(batch));
    }

    // queue for batches awaiting processing
    batch_queue workQueue(queueSize);
    moodycamel::ProducerToken ptok(workQueue);

    // spawn consumers
    std::vector<std::future<void>> threads;
    for(int i = 0; i < opt.numThreads-1; ++i) {
        threads.emplace_back(std::async(std::launch::async, [&] {
            sequence_batch batch;
            while(!queryingDone.load() || workQueue.size_approx() > 0) {
                if(workQueue.try_dequeue_from_producer(ptok, batch)) {
                    process_sequence_batch(
                        batch, db, finalizeMtx,
                        bufsrc, bufupdate, bufsink);

                    // put batch storage back
                    storageQueue.enqueue(std::move(batch));
                }
                else {
                    std::this_thread::sleep_for(std::chrono::milliseconds{10});
                }
            }
        }));
    }

    query_id endId = startId;

    // read sequences from file
    try {
        sequence_pair_reader reader{filename1, filename2};
        reader.index_offset(startId);

        sequence_batch batch;

        while(reader.has_next()) {
            if(opt.numThreads > 1) {
                // get batch storage
                while(!storageQueue.try_dequeue(batch)) {
                    std::this_thread::sleep_for(std::chrono::milliseconds{10});
                }
            }

            if(batch.size() < opt.batchSize) { batch.resize(opt.batchSize); }

            // fill batch with sequences
            std::size_t i = 0;
            for(auto& query : batch) {
                if(i > queryLimit || !reader.has_next()) break;

                query.id = reader.next_header_and_data(query.header, query.seq1, query.seq2);
                query.groundTruth = nullptr;

                ++i;
            }

            if(i < batch.size()) { batch.resize(i); }

            if(!batch.empty()) {
                // either enqueue batch...
                if(opt.numThreads > 1) {
                    while(!workQueue.try_enqueue(ptok, std::move(batch))) {
                        std::this_thread::sleep_for(std::chrono::milliseconds{10});
                    };
                }
                // ... or process directly if single threaded
                else {
                    process_sequence_batch(
                        batch, db, finalizeMtx,
                        bufsrc, bufupdate, bufsink);
                }
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

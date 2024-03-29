/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2024 André Müller (muellan@uni-mainz.de)
 *                       & Robin Kobus  (kobus@uni-mainz.de)
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

#ifndef MC_DATABASE_QUERY_H_
#define MC_DATABASE_QUERY_H_


#include "batch_processing.h"
#include "candidate_generation.h"
#include "cmdline_utility.h"
#include "database.h"
#include "options.h"
#include "sequence_io.h"

#ifndef GPU_MODE
    #include "query_handler.h"
#else
    #include "query_batch.cuh"
#endif

#include <iostream>
#include <vector>


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
        seq1(std::move(s1)), seq2(std::move(s2))
    {}

    bool empty() const noexcept { return header.empty() || seq1.empty(); }

    query_id id = 0;
    std::string header;
    sequence seq1;
    sequence seq2;  // 2nd part of paired-end read
};



/*************************************************************************//**
 *
 * @brief process batch of results from database query
 *
 * @tparam Buffer        batch buffer object
 *
 * @tparam BufferUpdate  takes database matches of one query and a buffer;
 *                       must be thread-safe (only const operations on DB!)
 *
 *****************************************************************************/
#ifndef GPU_MODE
    template<class Buffer, class BufferUpdate>
    void query_host(
        const database& db,
        const query_options& opt,
        const std::vector<sequence_query>& batch,
        query_handler<location>& queryHandler,
        Buffer& resultsBuffer, BufferUpdate& update)
    {
        for (const auto& query : batch) {
            auto rules = make_candidate_generation_rules(
                query, opt.classify, db.target_sketching().winstride);

            db.query_host(query.seq1, query.seq2, queryHandler, opt.sketching, rules);

            update(resultsBuffer, query, queryHandler.allhits(), queryHandler.tophits());
        }
    }
#else
    template<class Buffer, class BufferUpdate>
    void query_gpu(
        const database& db,
        const query_options& opt,
        const std::vector<sequence_query>& sequenceBatch,
        query_batch<location>& queryBatch,
        part_id hostId,
        Buffer& batchBuffer, BufferUpdate& update,
        std::mutex& scheduleMtx)
    {
        for (const auto& query : sequenceBatch) {
            auto rules = make_candidate_generation_rules(
                query, opt.classify, db.target_sketching().winstride);

            if (!queryBatch.add_paired_read(hostId, query.seq1, query.seq2,
                                           opt.sketching, rules))
            {
                std::cerr << "query batch is too small for a single read!\n";
            }
        }

        if (queryBatch.host_data(hostId).num_windows() > 0)
        {
            {
                std::lock_guard<std::mutex> lock(scheduleMtx);
                db.query_gpu_async(queryBatch, hostId, opt.sketching, opt.classify.lowestRank);
            }
            queryBatch.host_data(hostId).wait_for_results();

            for (size_t s = 0; s < queryBatch.host_data(hostId).num_queries(); ++s) {
                update(batchBuffer, sequenceBatch[s],
                       queryBatch.host_data(hostId).allhits(s),
                       queryBatch.host_data(hostId).top_candidates(s));
            }

            queryBatch.host_data(hostId).clear();
        }
    }
#endif



 /*************************************************************************//**
 *
 * @brief queries database with batches of reads from ONE sequence source (pair)
 *        produces batch buffers with one match list per sequence
 *
 * @tparam BufferSource     returns a per-batch buffer object
 *
 * @tparam BufferUpdate     takes database matches of one query and a buffer;
 *                          must be thread-safe (only const operations on DB!)
 *
 * @tparam BufferSink       recieves buffer after batch is finished
 *
 * @tparam InfoCallback     prints messages
 *
 * @tparam ProgressHandler  prints progress messages
 *
 * @tparam ErrorHandler     handles exceptions
 *
 *****************************************************************************/
template<
    class BufferSource, class BufferUpdate, class BufferSink,
    class ErrorHandler
>
query_id query_batched(
    const std::string& filename1, const std::string& filename2,
    const database& db,
    const query_options& opt,
    query_id idOffset,
    BufferSource&& getBuffer, BufferUpdate&& update, BufferSink&& finalize,
    ErrorHandler&& handleErrors)
{
    if (opt.performance.queryLimit < 1) return idOffset;
    size_t queryLimit = opt.performance.queryLimit > 0 ?
                        size_t(opt.performance.queryLimit) :
                        std::numeric_limits<size_t>::max();

    std::mutex finalizeMtx;

    const unsigned numWorkers = opt.performance.numThreads - (opt.performance.numThreads > 1);

#ifndef GPU_MODE
    std::vector<query_handler<location>> queryHandlers;
    queryHandlers.resize(numWorkers);

#else
    std::vector<std::mutex> scheduleMtxs(opt.performance.replication);

    std::vector<query_batch<location>> queryBatches;
    queryBatches.reserve(opt.performance.replication);

    for (unsigned rep = 0; rep < opt.performance.replication; ++rep)
        queryBatches.emplace_back(
            opt.performance.batchSize,
            opt.performance.batchSize*opt.sketching.winlen,
            opt.sketching.sketchlen,
            opt.sketching.sketchlen*db.max_locations_per_feature(),
            opt.classify.maxNumCandidatesPerQuery,
            opt.output.analysis.showAllHits,
            numWorkers,
            db.num_parts(),
            rep);
#endif

    // get executor that runs classification in batches
    batch_processing_options<sequence_query> execOpt;
    execOpt.concurrency(1, numWorkers);
    execOpt.batch_size(opt.performance.batchSize);
    execOpt.queue_size(opt.performance.numThreads + 8);
    execOpt.on_error(handleErrors);
    execOpt.work_item_measure([&] (const auto& query) {
        using std::begin;
        using std::end;
        using std::distance;

        const numk_t kmerSize = opt.sketching.kmerlen;
        const size_t windowStride = opt.sketching.winstride;

        const size_t seqLength1 = distance(begin(query.seq1), end(query.seq1));
        const size_t seqLength2 = distance(begin(query.seq2), end(query.seq2));

        const window_id numWindows1 = (seqLength1 >= kmerSize) ?
                                      (seqLength1-kmerSize + windowStride) / windowStride : 0;
        const window_id numWindows2 = (seqLength2 >= kmerSize) ?
                                      (seqLength2-kmerSize + windowStride) / windowStride : 0;

        return std::max(window_id(1), numWindows1 + numWindows2);
    });

    batch_executor<sequence_query> executor {
        execOpt,
        // classifies a batch of input queries
        [&](int id, std::vector<sequence_query>& batch) {
            auto resultsBuffer = getBuffer();

#ifndef GPU_MODE
            query_host(db, opt, batch, queryHandlers[id], resultsBuffer, update);
#else
            // query batch to gpu and wait for results
            query_gpu(db, opt, batch,
                      queryBatches[id%opt.performance.replication], id/opt.performance.replication,
                      resultsBuffer, update, scheduleMtxs[id%opt.performance.replication]);
#endif

            std::lock_guard<std::mutex> lock(finalizeMtx);
            finalize(std::move(resultsBuffer));

            return true;
        }};

    std::size_t discardedCount = 0;
    std::size_t totalReadCount = 0;

    // read sequences from file
    try {
        sequence_pair_reader reader{filename1, filename2};
        reader.index_offset(idOffset);

        while (reader.has_next()) {
            if (queryLimit < 1) break;

            // get (ref to) next query sequence storage and fill it
            auto& query = executor.next_item();
            query.id = reader.next_header_and_data(query.header, query.seq1, query.seq2);
            ++totalReadCount;

            // read length filter
            while (query.seq1.size() < opt.minReadLength ||
                   query.seq1.size() > opt.maxReadLength)
            {
                ++discardedCount;
                if (!reader.has_next()) break;
                query.id = reader.next_header_and_data(query.header, query.seq1, query.seq2);
                ++totalReadCount;
            }

            --queryLimit;
        }

        idOffset = reader.index();
    }
    catch(std::exception& e) {
        handleErrors(e);
    }

    if (discardedCount > 0) {
        std::size_t readsUsedCount =  totalReadCount - discardedCount;
        std::cerr << "\nRead length filter ("
                << opt.minReadLength << " - " << opt.maxReadLength << ") active\n"
                << "total:     " << totalReadCount << "\n"
                << "discarded: " << discardedCount << " ("
                << (100.0 * double(discardedCount)/double(totalReadCount)) << "%)\n"
                << "used:      " << readsUsedCount << " ("
                << (100.0 * double(readsUsedCount)/double(totalReadCount)) << "%)\n";
    }

    return idOffset;
}




 /*************************************************************************//**
 *
 * @brief queries database with batches of reads from multiple sequence sources
 *
 * @tparam BufferSource     returns a per-batch buffer object
 *
 * @tparam BufferUpdate     takes database matches of one query and a buffer;
 *                          must be thread-safe (only const operations on DB!)
 *
 * @tparam BufferSink       recieves buffer after batch is finished
 *
 * @tparam InfoCallback     prints messages
 *
 * @tparam ProgressHandler  prints progress messages
 *
 * @tparam ErrorHandler     handles exceptions
 *
 *****************************************************************************/
template<
    class BufferSource, class BufferUpdate, class BufferSink,
    class InfoCallback, class ProgressHandler, class ErrorHandler
>
void query_database(
    const std::vector<std::string>& infilenames,
    const database& db,
    const query_options& opt,
    BufferSource&& bufsrc, BufferUpdate&& bufupdate, BufferSink&& bufsink,
    InfoCallback&& showInfo, ProgressHandler&& showProgress,
    ErrorHandler&& errorHandler)
{
    const size_t stride = opt.pairing == pairing_mode::files ? 1 : 0;
    const std::string nofile;
    query_id queryIdOffset = 0;

    // input filenames passed to sequence reader depend on pairing mode:
    // none     -> infiles[i], ""
    // sequence -> infiles[i], infiles[i]
    // files    -> infiles[i], infiles[i+1]

    for (size_t i = 0; i < infilenames.size(); i += stride+1) {
        // pair up reads from two consecutive files in the list
        const auto& fname1 = infilenames[i];

        const auto& fname2 = (opt.pairing == pairing_mode::none)
                             ? nofile : infilenames[i+stride];

        if (opt.pairing == pairing_mode::files) {
            showInfo(fname1 + " + " + fname2);
        } else {
            showInfo(fname1);
        }
        showProgress(infilenames.size() > 1 ? i/float(infilenames.size()) : -1);

        queryIdOffset = query_batched(fname1, fname2, db, opt, queryIdOffset,
                                     std::forward<BufferSource>(bufsrc),
                                     std::forward<BufferUpdate>(bufupdate),
                                     std::forward<BufferSink>(bufsink),
                                     errorHandler);
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
 * @tparam InfoCallback  prints status messages
 *
 *****************************************************************************/
template<
    class BufferSource, class BufferUpdate, class BufferSink, class InfoCallback
>
void query_database(
    const std::vector<std::string>& infilenames,
    const database& db,
    const query_options& opt,
    BufferSource&& bufsrc, BufferUpdate&& bufupdate, BufferSink&& bufsink,
    InfoCallback&& showInfo)
{
    query_database(infilenames, db, opt,
       std::forward<BufferSource>(bufsrc),
       std::forward<BufferUpdate>(bufupdate),
       std::forward<BufferSink>(bufsink),
       std::forward<InfoCallback>(showInfo),
       [] (float p) { show_progress_indicator(std::cerr, p); },
       [] (std::exception& e) { std::cerr << "FAIL: " << e.what() << '\n'; }
    );
}


} // namespace mc


#endif

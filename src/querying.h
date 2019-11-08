/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2019 André Müller (muellan@uni-mainz.de)
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

#ifndef MC_QUERYING_H_
#define MC_QUERYING_H_

#include <vector>
#include <iostream>

#include "sketch_database.h"
#include "sequence_io.h"
#include "cmdline_utility.h"
#include "query_options.h"
#include "cmdline_utility.h"
#include "batch_processing.h"
#include "query_batch.cuh"

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
    const database& db, const query_processing_options& opt,
    query_id idOffset,
    BufferSource&& getBuffer, BufferUpdate&& update, BufferSink&& finalize,
    ErrorHandler&& handleErrors)
{
    if(opt.queryLimit < 1) return idOffset;
    auto queryLimit = size_t(opt.queryLimit > 0 ? opt.queryLimit : std::numeric_limits<size_t>::max());

    std::mutex finalizeMtx;

    // each thread may have a gpu batch for querying
    std::vector<std::unique_ptr<query_batch<location>>> gpuBatches(opt.numThreads);

    // get executor that runs classification in batches
    batch_processing_options execOpt;
    execOpt.concurrency(opt.numThreads - 1);
    execOpt.batch_size(opt.batchSize);
    execOpt.queue_size(opt.numThreads > 1 ? opt.numThreads + 4 : 0);
    execOpt.on_error(handleErrors);

    batch_executor<sequence_query> executor {
        execOpt,
        // classifies a batch of input queries
        [&](int id, std::vector<sequence_query>& batch) {
            auto resultsBuffer = getBuffer();
            std::vector<database::matches_sorter> targetMatches(batch.size());
            std::vector<database::matches_sorter> targetMatchesTmp(batch.size());

            size_t outIndex = 0;

            if(!gpuBatches[id])
                gpuBatches[id] = std::make_unique<query_batch<location>>(
                    5, 100, db.query_sketcher().sketch_size()*db.max_locations_per_feature());

            for(size_t i = 0; i < batch.size(); ++i) {
                auto& seq = batch[i];

                auto numQueries = gpuBatches[id]->num_queries();

                if(!gpuBatches[id]->add_paired_read(seq.id, seq.seq1, seq.seq2, db.query_sketcher())) {
                    //could not add read, send full batch to gpu
                    if(gpuBatches[id]->num_queries() > 0) {
                        // std::cout << "send batch to gpu\n";
                        db.accumulate_matches_gpu(*(gpuBatches[id]), targetMatches, targetMatchesTmp, outIndex);
                        gpuBatches[id]->clear();

                        //try to add read again
                        if(!gpuBatches[id]->add_paired_read(seq.id, seq.seq1, seq.seq2, db.query_sketcher()))
                            std::cerr << "query batch is too small for a single read!" << std::endl;
                    }
                    else {
                        std::cerr << "query batch is too small for a single read!" << std::endl;
                    }
                }

                // std::cout << i << ". sequence: ";
                if(gpuBatches[id]->num_queries() != numQueries) {
                    // std::cout << gpuBatches[id]->num_queries() << " queries\n";
                }
                else {
                    // std::cout << "no queries\n";
                    database::match_locations locations;
                    update(resultsBuffer, seq, locations);
                    batch.erase(batch.begin() + i);
                }
            }

            if(gpuBatches[id]->num_queries() > 0) {
                // std::cout << "send final batch to gpu\n";
                db.accumulate_matches_gpu(*(gpuBatches[id]), targetMatches, targetMatchesTmp, outIndex);
                gpuBatches[id]->clear();
            }

            for(size_t i = 0; i < batch.size(); ++i) {
                auto& seq = batch[i];

                auto it = is_sorted_until(targetMatches[i].begin(), targetMatches[i].end());
                auto numSorted = std::distance(targetMatches[i].begin(), it);
                // std::cout << "sorted: " << numSorted << '\n';

                if(numSorted < 4064) {
                    std::cout << i << ". sorted: " << numSorted << '\n';

                    targetMatchesTmp[i].sort();

                    std::cout << i << ". targetMatchesTmp: ";
                    for(const auto& m : targetMatchesTmp[i]) {
                        if(m.tgt < std::numeric_limits<target_id>::max()) {
                            std::cout << m.tgt << ':' << m.win << ' ';
                        }
                        else {
                            break;
                        }
                    }
                    std::cout << '\n';


                    std::cout << i << ". targetMatches:    ";
                    for(const auto& m : targetMatches[i]) {
                        if(m.tgt < std::numeric_limits<target_id>::max()) {
                            std::cout << m.tgt << ':' << m.win << ' ';
                        }
                        else {
                            break;
                        }
                    }
                    std::cout << '\n';
                }

                auto it2 = std::find_if(targetMatches[i].begin(), targetMatches[i].end(), [](const auto& loc) {
                    return loc.tgt == std::numeric_limits<target_id>::max();
                });
                targetMatches[i].resize(std::distance(targetMatches[i].begin(), it2));

                update(resultsBuffer, seq, targetMatches[i].locations());
            }

            std::lock_guard<std::mutex> lock(finalizeMtx);
            finalize(std::move(resultsBuffer));
        }};

    // read sequences from file
    try {
        sequence_pair_reader reader{filename1, filename2};
        reader.index_offset(idOffset);

        while(reader.has_next()) {
            if(queryLimit < 1) break;

            // get (ref to) next query sequence storage and fill it
            auto& query = executor.next_item();
            query.id = reader.next_header_and_data(query.header, query.seq1, query.seq2);

            --queryLimit;
        }

        idOffset = reader.index();
    }
    catch(std::exception& e) {
        handleErrors(e);
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
    const query_processing_options& opt,
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
       [] (std::exception& e) { std::cerr << "FAIL: " << e.what() << '\n'; }
    );
}


} // namespace mc


#endif

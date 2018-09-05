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
    parallel_queue& queue,
    const std::string& filename1, const std::string& filename2,
    const database& db, const query_processing_options& opt,
    query_id& startId,
    BufferSource&& getBuffer, BufferUpdate&& update, BufferSink&& finalize)
{
    const auto load = 32 * queue.concurrency();

    std::mutex mtx;
    std::atomic<std::int_least64_t> queryLimit{opt.queryLimit};
    std::atomic<query_id> workId{startId};
    std::atomic<bool> done{false};

    std::mutex tipMtx;
    sequence_pair_reader::stream_positions pos{0,0};
    query_id qid{startId};

    while(!done && queryLimit > 0) {
        queue.enqueue_when([&]{return queue.unsafe_waiting() < load;}, [&] {
            sequence_pair_reader reader{filename1, filename2};

            auto buffer = getBuffer();
            match_locations matches;
            bool empty = true;

            for(std::size_t i = 0; i < opt.batchSize && queryLimit-- > 0; ++i) {
                sequence_pair_reader::stream_positions myPos;
                query_id myQid;

                // get most recent position and query id
                {
                    std::lock_guard<std::mutex> lock(tipMtx);
                    myPos = pos;
                    myQid = qid;
                }
                reader.seek(myPos);
                if(!reader.has_next()) break;
                reader.index_offset(myQid);

                // get work id and skip to this read
                auto wid = workId.fetch_add(1);
                if(myQid != wid) {
                    reader.skip(wid-myQid);
                }
                if(!reader.has_next()) break;
                auto seq = reader.next();

                // update most recent position and query id
                myPos = reader.tell();
                myQid = reader.index();
                if(myQid > qid) {// reading qid unsafe
                    std::lock_guard<std::mutex> lock(tipMtx);
                    if(myQid > qid) {// reading qid safe
                        pos = myPos;
                        qid = myQid;
                    }
                }

                if(!seq.first.header.empty()) {
                    empty = false;
                    matches.clear();

                    db.accumulate_matches(seq.first.data, matches);
                    db.accumulate_matches(seq.second.data, matches);

                    std::sort(matches.begin(), matches.end());

                    update(buffer,
                           sequence_query{seq.first.index,
                                          std::move(seq.first.header),
                                          std::move(seq.first.data),
                                          std::move(seq.second.data)},
                           matches );
                }
            }
            if(!empty) {
                std::lock_guard<std::mutex> lock(mtx);
                finalize(std::move(buffer));
            }

            if(!reader.has_next()) done = true;
        }); //enqueue
    }
    //wait for all enqueued tasks to finish
    queue.wait();

    startId = qid;
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

    const size_t stride = opt.pairing == pairing_mode::files ? 1 : 0;
    const std::string nofile;
    query_id readIdOffset = 0;

    for(size_t i = 0; i < infilenames.size(); i += stride+1) {
        //pair up reads from two consecutive files in the list
        const auto& fname1 = infilenames[i];

        const auto& fname2 = (opt.pairing == pairing_mode::none)
                             ? nofile : infilenames[i+stride];

        const auto progress = infilenames.size() > 1 ? i / float(infilenames.size()) : -1;
        if(opt.pairing == pairing_mode::files) {
            showStatus(fname1 + " + " + fname2, progress);
        } else {
            showStatus(fname1, progress);
        }

        try {
            query_batched(queue, fname1, fname2, db, opt, readIdOffset,
                          std::forward<BufferSource>(bufsrc),
                          std::forward<BufferUpdate>(bufupdate),
                          std::forward<BufferSink>(bufsink));
        }
        catch(std::exception& e) {
            showStatus(std::string("FAIL: ") + e.what(), progress);
        }
    }
}


} // namespace mc


#endif

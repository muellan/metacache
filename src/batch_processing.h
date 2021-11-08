#ifndef BATCH_PROCESSING_H_
#define BATCH_PROCESSING_H_


#include "../dep/queue/concurrentqueue.h"

#include <atomic>
#include <chrono>
#include <future>
#include <vector>


namespace mc {

// forward declarations
template<class WorkItem> class batch_executor;

/*************************************************************************//**
 *
 * @brief configuration for batch_executor
 *
 *****************************************************************************/
template<class WorkItem>
class batch_processing_options {
public:
    friend class batch_executor<WorkItem>;

    using error_handler   = std::function<void(std::exception&)>;
    using finalizer       = std::function<void(int)>;
    using item_measure    = std::function<size_t(const WorkItem&)>;

    batch_processing_options():
        numProducers_{1},
        numConsumers_{1},
        queueSize_{1},
        batchSize_{1},
        handleErrors_{[](std::exception&){}},
        finalize_{[](int){}},
        measureWorkItem_{[](const WorkItem&){ return 1; }}
    {}

    int num_producers()      const noexcept { return numProducers_; }
    int num_consumers()      const noexcept { return numConsumers_; }
    std::size_t batch_size() const noexcept { return batchSize_; }
    std::size_t queue_size() const noexcept { return queueSize_; }

    void concurrency(int p, int c) noexcept {
        numProducers_ = p >= 0 ? p : 1;
        numConsumers_ = c >= 0 ? c : 1;
    }
    void batch_size(std::size_t n) noexcept { batchSize_  = n > 0 ? n : 1; }
    void queue_size(std::size_t n) noexcept { queueSize_  = n > 0 ? n : 1; }

    void on_work_done(finalizer f)   { finalize_ = std::move(f); }
    void on_error(error_handler f)   { handleErrors_ = std::move(f); }
    void work_item_measure(item_measure f) { measureWorkItem_ = std::move(f); }

private:
    int numProducers_;
    int numConsumers_;
    std::size_t queueSize_;
    std::size_t batchSize_;
    error_handler handleErrors_;
    finalizer finalize_;
    item_measure measureWorkItem_;
};



/*************************************************************************//**
 *
 * @brief  single producer, multiple consumer parallel batch processing;
 *         uses batch storage recycling strategy to avoid frequent allocations;
 *         runs sequentially if concurrency is set to 0
 *
 * @tparam WorkItem
 *
 *****************************************************************************/
template<class WorkItem>
class batch_executor
{
public:

    using batch_type      = std::vector<WorkItem>;
    using batch_consumer  = std::function<bool(int,batch_type&)>;
    using error_handler   = typename batch_processing_options<WorkItem>::error_handler;
    using finalizer       = typename batch_processing_options<WorkItem>::finalizer;

private:
    struct batch_handler {
        std::size_t workCount_{0};
        std::size_t workSize_{0};
        batch_type batch_{};
        bool finalized_{false};
    };

public:
    // -----------------------------------------------------------------------
    /**
     * @param consume       processes a batch of work items
     * @param handleErrors  handles exeptions
     * @param finalize      runs after consumer is finished
     */
    batch_executor(batch_processing_options<WorkItem> opt,
                   batch_consumer consume)
    :
        param_{std::move(opt)},
        keepWorking_{param_.num_consumers()},
        storageQueue_{param_.queue_size()},
        workQueue_{param_.queue_size()},
        producerBatches_(param_.num_producers()),
        consumers_{},
        consume_{std::move(consume)}
    {
        for(auto& handler : producerBatches_)
            handler.batch_.resize(param_.batch_size());

        // fill batch storage queue with initial batches
        // we want to re-use the individual batches and all their members
        // => minimize/avoid batch-resizing during processing
        for(std::size_t i = 0; i < param_.queue_size(); ++i) {
            batch_type batch;
            batch.resize(param_.batch_size());
            storageQueue_.enqueue(std::move(batch));
        }

        // spawn consumer threads
        consumers_.reserve(param_.num_consumers());
        for(int consumerId = 0; consumerId < param_.num_consumers(); ++consumerId) {
            consumers_.emplace_back(std::async(std::launch::async, [&,consumerId] {
                batch_type batch;
                while(valid() || workQueue_.size_approx() > 0) {
                    if(workQueue_.try_dequeue(batch)) {
                        if(consume_(consumerId, batch)) {
                            // batch processed completely
                            // put batch storage back
                            storageQueue_.enqueue(std::move(batch));
                        }
                        else {
                            // work remaining in batch
                            // put batch back
                            workQueue_.enqueue(std::move(batch));
                            keepWorking_.fetch_sub(1);
                            break;
                        }
                    }
                    else {
                        std::this_thread::sleep_for(std::chrono::milliseconds{1});
                    }
                }
                param_.finalize_(consumerId);
            }));
        }
    }


    // -----------------------------------------------------------------------
    ~batch_executor() {
        try {
            for(int producerId = 0; producerId < param_.num_producers(); ++producerId) {
                finalize_producer(producerId);
            }

            // signal all consumers to finish as soon as no work is left
            keepWorking_.store(0);
            // wait until consumers are finished
            for(auto& consumer : consumers_) {
                if(consumer.valid()) consumer.get();
            }
        }
        catch(std::exception& e) {
            param_.handleErrors_(e);
        }
    }


    // -----------------------------------------------------------------------
    const batch_processing_options<WorkItem>&
    param() const noexcept {
        return param_;
    }


    // -----------------------------------------------------------------------
    /** @return false, if invalidated or destruction in progress */
    bool valid() const noexcept {
        return keepWorking_.load() > 0;
    }


    // -----------------------------------------------------------------------
    void invalidate() noexcept {
        return keepWorking_.store(0);
    }


    // -----------------------------------------------------------------------
    /** @brief  get reference to next work item */
    WorkItem& next_item(int producerId = 0) {
        auto& handler = producerBatches_[producerId];

        consume_batch_if_full(producerId);

        if(handler.workCount_ >= handler.batch_.size()) {
            handler.batch_.resize(handler.workCount_+1);
        }
        return handler.batch_[handler.workCount_++];
    }


    // -----------------------------------------------------------------------
    /** @brief  consume last batch if not empty */
    void finalize_producer(int producerId = 0) {
        auto& handler = producerBatches_[producerId];

        if(!handler.finalized_ && valid()) {
            consume_batch_if_full(producerId);

            // finalize last batch
            if(handler.workCount_ > 0) {
                handler.batch_.resize(handler.workCount_);
                consume_batch(producerId);
                handler.workCount_ = 0;
            }

            handler.finalized_ = true;
        }
    }


private:
    // -----------------------------------------------------------------------
    void consume_batch_if_full(int producerId) {
        auto& handler = producerBatches_[producerId];

        if(handler.workCount_ > 0) {
            const std::size_t lastSize = param_.measureWorkItem_(handler.batch_[handler.workCount_-1]);
            handler.workSize_ += lastSize;

            if(handler.workSize_ >= param_.batch_size()) {
                WorkItem lastItem;
                // if batch too big, remove last item
                if(handler.workSize_ > param_.batch_size())
                    lastItem = std::move(handler.batch_[--handler.workCount_]);

                handler.batch_.resize(handler.workCount_);

                consume_batch(producerId);

                // get new batch storage
                while(valid() && !storageQueue_.try_dequeue(handler.batch_)) {
                    std::this_thread::sleep_for(std::chrono::milliseconds{1});
                }

                // reinsert last item
                if(handler.workSize_ > param_.batch_size()) {
                    handler.batch_.front() = std::move(lastItem);
                    handler.workSize_ = lastSize;
                    handler.workCount_ = 1;
                }
                else {
                    handler.workSize_ = 0;
                    handler.workCount_ = 0;
                }
            }
        }
    }


    // -----------------------------------------------------------------------
    void consume_batch(int producerId) {
        auto& handler = producerBatches_[producerId];

        workQueue_.enqueue(std::move(handler.batch_));
    }


    // -----------------------------------------------------------------------
    using batch_queue = moodycamel::ConcurrentQueue<batch_type>;

    const batch_processing_options<WorkItem> param_;
    std::atomic_int keepWorking_;
    batch_queue storageQueue_;
    batch_queue workQueue_;
    std::vector<batch_handler> producerBatches_;
    std::vector<std::future<void>> consumers_;
    batch_consumer consume_;
};


} // namespace mc


#endif

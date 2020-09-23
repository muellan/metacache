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
    using abort_condition = std::function<bool(int)>;
    using finalizer       = std::function<void(int)>;
    using item_measure    = std::function<size_t(const WorkItem&)>;

    batch_processing_options():
        numProducers_{0},
        numConsumers_{0},
        queueSize_{1},
        batchSize_{1},
        handleErrors_{[](std::exception&){}},
        abortRequested_{[](int){ return false; }},
        finalize_{[](int){}},
        measureWorkItem_{[](const WorkItem&){ return 1; }}
    {}

    int num_producers()      const noexcept { return numProducers_; }
    int num_consumers()      const noexcept { return numConsumers_; }
    std::size_t batch_size() const noexcept { return batchSize_; }
    std::size_t queue_size() const noexcept { return queueSize_; }

    void concurrency(int p, int c) noexcept {
        numProducers_ = p >= 0 ? p : 1;
        numConsumers_ = c >= 0 ? c : 0;
    }
    void batch_size(std::size_t n) noexcept { batchSize_  = n > 0 ? n : 1; }
    void queue_size(std::size_t n) noexcept { queueSize_  = n > 0 ? n : 1; }

    void on_work_done(finalizer f)   { finalize_ = std::move(f); }
    void on_error(error_handler f)   { handleErrors_ = std::move(f); }
    void abort_if(abort_condition f) { abortRequested_ = std::move(f); }
    void work_item_measure(item_measure f) { measureWorkItem_ = std::move(f); }

private:
    int numProducers_;
    int numConsumers_;
    std::size_t queueSize_;
    std::size_t batchSize_;
    error_handler handleErrors_;
    abort_condition abortRequested_;
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
    using batch_consumer  = std::function<void(int,batch_type&)>;
    using error_handler   = typename batch_processing_options<WorkItem>::error_handler;
    using abort_condition = typename batch_processing_options<WorkItem>::abort_condition;
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
        keepWorking_{true},
        storageQueue_{param_.queue_size()},
        workQueue_{param_.queue_size()},
        producerBatches_(param_.num_producers()),
        consumers_{},
        consume_{std::move(consume)}
    {
        for(auto& handler : producerBatches_)
            handler.batch_.resize(32);

        if(param_.num_consumers() > 0) {
            // fill batch storage queue with initial batches
            // we want to re-use the individual batches and all their members
            // => minimize/avoid batch-resizing during processing
            for(std::size_t i = 0; i < param_.queue_size(); ++i) {
                batch_type batch;
                batch.resize(32);
                storageQueue_.enqueue(std::move(batch));
            }

            // spawn consumer threads
            consumers_.reserve(param_.num_consumers());
            for(int consumerId = 0; consumerId < param_.num_consumers(); ++consumerId) {
                consumers_.emplace_back(std::async(std::launch::async, [&,consumerId] {
                    batch_type batch;
                    validate(consumerId);
                    while(valid() || workQueue_.size_approx() > 0) {
                        if(workQueue_.try_dequeue(batch)) {
                            consume_(consumerId, batch);
                            // put batch storage back
                            storageQueue_.enqueue(std::move(batch));
                            validate(consumerId);
                        }
                        else {
                            std::this_thread::sleep_for(std::chrono::milliseconds{1});
                        }
                    }
                    param_.finalize_(consumerId);
                }));
            }
        }
    }


    // -----------------------------------------------------------------------
    ~batch_executor() {
        try {
            for(int producerId = 0; producerId < param_.num_producers(); ++producerId) {
                finalize_producer(producerId);
            }

            if(!consumers_.empty()) {
                // signal all consumers to finish as soon as no work is left
                keepWorking_.store(false);
                // wait until consumers are finished
                for(auto& consumer : consumers_) {
                    if(consumer.valid()) consumer.get();
                }
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
    /** @return false, if abort condition was met or destruction in progress */
    bool valid() const noexcept {
        return keepWorking_.load();
    }


    // -----------------------------------------------------------------------
    /** @brief  get reference to next work item */
    WorkItem& next_item(int producerId = 0) {
        auto& handler = producerBatches_[producerId];

        handler.finalized_ = false;

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

        if(!handler.finalized_) {
            consume_batch_if_full(producerId);

            // finalize last batch
            if(handler.workCount_ > 0) {
                handler.batch_.resize(handler.workCount_);
                consume_batch(producerId);
                handler.workCount_ = 0;
            }

            if(consumers_.empty()) {
                param_.finalize_(producerId);
            }

            handler.finalized_ = true;
        }
    }


private:
    // -----------------------------------------------------------------------
    void consume_batch_if_full(int producerId) {
        auto& handler = producerBatches_[producerId];

        if(handler.workCount_ > 0) {
            std::size_t lastSize = param_.measureWorkItem_(handler.batch_[handler.workCount_-1]);
            handler.workSize_ += lastSize;

            if(handler.workSize_ >= param_.batch_size()) {
                WorkItem lastItem;
                // if batch too big, remove last item
                if(handler.workSize_ > param_.batch_size())
                    std::swap(lastItem, handler.batch_[--handler.workCount_]);

                handler.batch_.resize(handler.workCount_);

                consume_batch(producerId);

                // get new batch storage
                if(!consumers_.empty()) {
                    while(!storageQueue_.try_dequeue(handler.batch_)) {
                        std::this_thread::sleep_for(std::chrono::milliseconds{1});
                    }
                }

                // reinsert last item
                if(handler.workSize_ > param_.batch_size()) {
                    std::swap(lastItem, handler.batch_.front());
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

        // either enqueue if using consumer threads...
        if(!consumers_.empty()) {
            workQueue_.enqueue(std::move(handler.batch_));
        }
        // ... or consume directly if only producer threads
        else {
            validate(producerId);
            if(valid()) consume_(producerId, handler.batch_);
        }
    }


    // -----------------------------------------------------------------------
    void validate(int id) {
        if(param_.abortRequested_(id)) {
            keepWorking_.store(false);
        }
    }


    // -----------------------------------------------------------------------
    using batch_queue = moodycamel::ConcurrentQueue<batch_type>;

    const batch_processing_options<WorkItem> param_;
    std::atomic_bool keepWorking_;
    batch_queue storageQueue_;
    batch_queue workQueue_;
    std::vector<batch_handler> producerBatches_;
    std::vector<std::future<void>> consumers_;
    batch_consumer consume_;
};


} // namespace mc


#endif

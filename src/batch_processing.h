#ifndef BATCH_PROCESSING_H_
#define BATCH_PROCESSING_H_

#include <vector>
#include <atomic>
#include <future>
#include <chrono>

#include "../dep/queue/concurrentqueue.h"


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
    using abort_condition = std::function<bool()>;
    using finalizer       = std::function<void(int)>;
    using item_measure    = std::function<size_t(const WorkItem&)>;

    batch_processing_options():
        numWorkers_{0},
        queueSize_{1},
        batchSize_{1},
        handleErrors_{[](std::exception&){}},
        abortRequested_{[]{ return false; }},
        finalize_{[](int){}},
        measureWorkItem_{[](const WorkItem&){ return 1; }}
    {}

    int concurrency()        const noexcept { return numWorkers_; }
    std::size_t batch_size() const noexcept { return batchSize_; }
    std::size_t queue_size() const noexcept { return queueSize_; }

    void concurrency(int n)        noexcept { numWorkers_ = n >= 0 ? n : 0; }
    void batch_size(std::size_t n) noexcept { batchSize_  = n > 0 ? n : 1; }
    void queue_size(std::size_t n) noexcept { queueSize_  = n > 0 ? n : 1; }

    void on_work_done(finalizer f)   { finalize_ = std::move(f); }
    void on_error(error_handler f)   { handleErrors_ = std::move(f); }
    void abort_if(abort_condition f) { abortRequested_ = std::move(f); }
    void work_item_measure(item_measure f) { measureWorkItem_ = std::move(f); }

private:
    int numWorkers_;
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


    // -----------------------------------------------------------------------
    /**
     * @param consume       processes a batch of work items
     * @param handleErrors  handles exeptions
     * @param finalize      runs after worker is finished
     */
    batch_executor(batch_processing_options<WorkItem> opt,
                   batch_consumer consume)
    :
        param_{std::move(opt)},
        keepWorking_{true},
        currentWorkCount_{0},
        currentWorkSize_{0},
        currentBatch_{},
        storageQueue_{param_.queue_size()},
        workQueue_{param_.queue_size()},
        prodToken_{workQueue_},
        consume_{std::move(consume)},
        workers_{}
    {
        currentBatch_.resize(32);

        if(param_.concurrency() > 0) {
            // fill batch storage queue with initial batches
            // we want to re-use the individual batches and all their members
            // => minimize/avoid batch-resizing during processing
            for(std::size_t i = 0; i < param_.queue_size(); ++i) {
                batch_type batch;
                batch.resize(32);
                storageQueue_.enqueue(std::move(batch));
            }

            // spawn consumer threads
            workers_.reserve(param_.concurrency());
            for(int i = 0; i < param_.concurrency(); ++i) {
                workers_.emplace_back(std::async(std::launch::async, [&,i] {
                    batch_type batch;
                    validate();
                    while(valid() || workQueue_.size_approx() > 0) {
                        if(workQueue_.try_dequeue_from_producer(prodToken_, batch)) {
                            consume_(i, batch);
                            // put batch storage back
                            storageQueue_.enqueue(std::move(batch));
                            validate();
                        }
                        else {
                            std::this_thread::sleep_for(std::chrono::milliseconds{1});
                        }
                    }
                    param_.finalize_(i);
                }));
            }
        }
    }


    // -----------------------------------------------------------------------
    ~batch_executor() {
        try {
            consume_batch_if_full();
            // finalize last batch
            if(currentWorkCount_ > 0) {
                currentBatch_.resize(currentWorkCount_);
                consume_current_batch();
            }

            if(workers_.empty()) {
                param_.finalize_(0);
            }
            else {
                // signal all workers to finish as soon as no work is left
                keepWorking_.store(false);
                // wait until workers are finished
                for(auto& worker : workers_) {
                    if(worker.valid()) worker.get();
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
    WorkItem& next_item() {
        consume_batch_if_full();

        if(currentWorkCount_ >= currentBatch_.size()) {
            currentBatch_.resize(currentWorkCount_+1);
        }
        return currentBatch_[currentWorkCount_++];
    }


private:
    // -----------------------------------------------------------------------
    void consume_batch_if_full() {
        if(currentWorkCount_ > 0) {
            std::size_t lastSize = param_.measureWorkItem_(currentBatch_[currentWorkCount_-1]);
            currentWorkSize_ += lastSize;

            if(currentWorkSize_ >= param_.batch_size()) {
                WorkItem lastItem;
                // if batch too big, remove last item
                if(currentWorkSize_ > param_.batch_size())
                    std::swap(lastItem, currentBatch_[--currentWorkCount_]);

                currentBatch_.resize(currentWorkCount_);

                consume_current_batch();

                // get new batch storage
                if(!workers_.empty()) {
                    while(!storageQueue_.try_dequeue(currentBatch_)) {
                        std::this_thread::sleep_for(std::chrono::milliseconds{1});
                    }
                }

                // reinsert last item
                if(currentWorkSize_ > param_.batch_size()) {
                    std::swap(lastItem, currentBatch_.front());
                    currentWorkSize_ = lastSize;
                    currentWorkCount_ = 1;
                }
                else {
                    currentWorkSize_ = 0;
                    currentWorkCount_ = 0;
                }
            }
        }
    }


    // -----------------------------------------------------------------------
    void consume_current_batch() {
        // either enqueue if multi-threaded...
        if(!workers_.empty()) {
            workQueue_.enqueue(prodToken_, std::move(currentBatch_));
        }
        // ... or consume directly if single-threaded
        else {
            validate();
            if(valid()) consume_(0, currentBatch_);
        }
    }


    // -----------------------------------------------------------------------
    void validate() {
        if(param_.abortRequested_()) {
            keepWorking_.store(false);
        }
    }


    // -----------------------------------------------------------------------
    using batch_queue = moodycamel::ConcurrentQueue<batch_type>;

    const batch_processing_options<WorkItem> param_;
    std::atomic_bool keepWorking_;
    std::size_t currentWorkCount_;
    std::size_t currentWorkSize_;
    batch_type currentBatch_;
    batch_queue storageQueue_;
    batch_queue workQueue_;
    moodycamel::ProducerToken prodToken_;
    batch_consumer consume_;
    std::vector<std::future<void>> workers_;
};


} // namespace mc


#endif

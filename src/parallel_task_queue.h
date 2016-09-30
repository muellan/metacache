#ifndef MC_PARALLEL_TASK_QUEUE_H_
#define MC_PARALLEL_TASK_QUEUE_H_

#include <deque>
#include <vector>
#include <algorithm>
#include <iostream>
#include <thread>
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <functional>


namespace mc {


/*****************************************************************************
 *
 * @brief  thread that executes a single task and waits after completion
 *         until it is needed again
 *
 *****************************************************************************/
template<class Task>
class task_thread
{
    enum class status {busy = 0, idle = 1, terminate = 2 };

public:
    //---------------------------------------------------------------
    using task_type = Task;


    //---------------------------------------------------------------
    explicit
    task_thread():
        mutables_{}, wakeup_{}, isIdle_{},
        status_{status::idle},
        task_{},
        thread_{ [&] { work(); }}
    {}

    //---------------------------------------------------------------
    template<
        class Arg, class... Args, class =
        typename std::enable_if<!std::is_same<typename std::decay<Arg>::type,task_thread>::value,Arg>::type
    >
    explicit
    task_thread(Arg&& arg, Args&&... args):
        mutables_{}, wakeup_{}, isIdle_{},
        status_{status::idle},
        task_{std::forward<Arg>(arg), std::forward<Args>(args)...},
        thread_{ [&] { work(); }}
    {}


    //---------------------------------------------------------------
    ~task_thread() {
        //signal thread to terminate
        status_.store(status::terminate);
        wakeup_.notify_all();
        isIdle_.notify_all();

        wait_until_available();

        thread_.join();
    }


    //---------------------------------------------------------------
    task_thread(const task_thread& source):
        mutables_{}, wakeup_{}, isIdle_{},
        status_{status::idle},
        task_{source.task_},
        thread_{ [&] { work(); }}
    {}

    //---------------------------------------------------------------
    task_thread(task_thread&& source)
        noexcept(noexcept(task_type(std::declval<task_type>())))
    :
        mutables_{}, wakeup_{}, isIdle_{},
        status_{status::idle},
        task_{std::move(source.task_)},
        thread_{ [&] { work(); }}
    {}


    //---------------------------------------------------------------
    task_thread& operator = (const task_thread& source)
    {
        wait_until_available();

        std::lock_guard<std::mutex> lock(mutables_);
        task_ = source.task_;

        return *this;
    }

    //---------------------------------------------------------------
    task_thread& operator = (task_thread&& source)
        noexcept(noexcept(std::declval<task_type>() = source.task_))
    {
        wait_until_available();

        std::lock_guard<std::mutex> lock(mutables_);
        task_ = std::move(source.task_);

        return *this;
    }


    //---------------------------------------------------------------
    bool available() const noexcept {
        return bool(status_.load());
    }


    //---------------------------------------------------------------
    /**
     * @brief block execution of calling thread until task is finished
     */
    void wait_until_available()
    {
        std::unique_lock<std::mutex> lock(mutables_);

        while(status_.load() == status::busy) {
            isIdle_.wait(lock);
        }
    }


    //---------------------------------------------------------------
    const task_type&
    unsafe_task() const & noexcept {
        return task_;
    }
    task_type&
    unsafe_task() & noexcept {
        return task_;
    }
    task_type&&
    unsafe_task() && noexcept {
        return std::move(task_);
    }


    //---------------------------------------------------------------
    bool operator () ()
    {
        if(status_.load() != status::idle) return false;
        std::lock_guard<std::mutex> lock{mutables_};
        status_.store(status::busy);

        wakeup_.notify_all();
        return true;
    }

    //---------------------------------------------------------------
    bool operator () (const task_type& task)
    {
        if(status_.load() != status::idle) return false;
        std::lock_guard<std::mutex> lock{mutables_};
        status_.store(status::busy);

        task_ = task;

        wakeup_.notify_all();
        return true;
    }

    //---------------------------------------------------------------
    bool operator () (task_type&& task)
    {
        if(status_.load() != status::idle) return false;
        std::lock_guard<std::mutex> lock{mutables_};
        status_.store(status::busy);

        task_ = std::move(task);

        wakeup_.notify_all();
        return true;
    }

    //---------------------------------------------------------------
    template<class Arg, class... Args>
    bool operator () (Arg&& arg, Args&&... args)
    {
        if(status_.load() != status::idle) return false;
        std::lock_guard<std::mutex> lock{mutables_};
        status_.store(status::busy);

        task_ = task_type{std::forward<Arg>(arg),
                          std::forward<Args>(args)...};

        wakeup_.notify_all();
        return true;
    }


private:

    //---------------------------------------------------------------
    void work()
    {
         while(!abort_requested()) {
            wait_until_needed();

            if(!abort_requested()) {
                task_();
                status_.store(status::idle);
                isIdle_.notify_all();
            }
         }
         isIdle_.notify_all();
    }


    //---------------------------------------------------------------
    bool abort_requested() const noexcept {
        return (status_.load() == status::terminate);
    }


    //---------------------------------------------------------------
    void wait_until_needed() {
        std::unique_lock<std::mutex> lock{mutables_};

        //block execution of calling thread until woken up
        while(status_.load() == status::idle) {
            wakeup_.wait(lock);
        }
    }


    //---------------------------------------------------------------
    mutable std::mutex mutables_;
    std::condition_variable wakeup_;
    std::condition_variable isIdle_;
    std::atomic<status> status_;
    task_type task_;
    std::thread thread_;
};




/*****************************************************************************
 *
 * @brief    simple queue that holds tasks and executes them in parallel
 *
 * @tparam   Task : assumed to be cheap to copy or move
 *                  (e.g. std::function, std::reference_wrapper, ...)
 *                  Task::operator()() must be defined
 *
 *****************************************************************************/
template<class Task>
class parallel_task_queue
{
public:
    //---------------------------------------------------------------
    using task_type = Task;


private:
    //---------------------------------------------------------------
    class task_executor {
    public:
        explicit
        task_executor(parallel_task_queue* callback = nullptr,
                      task_type&& task = task_type{}) noexcept
        :
            queue_{callback},
            task_{std::move(task)}
        {}

        void operator () () {
            task_();
            --queue_->running_;
        }

    private:
        parallel_task_queue* queue_;
        task_type task_;
    };

    friend class task_executor;


public:
    //---------------------------------------------------------------
    explicit
    parallel_task_queue(unsigned int concurrency =
        std::thread::hardware_concurrency())
    :
        mutables_{},
        active_{true}, hasWaiting_{false},
        running_{0},
        waiting_{},
        workers_(concurrency),
        isDone_{},
        scheduler_{ [&] { schedule(); }}
    {}

    //-----------------------------------------------------
    parallel_task_queue(const parallel_task_queue&) = delete;
    parallel_task_queue(parallel_task_queue&&) = delete;


    //-----------------------------------------------------
    ~parallel_task_queue()
    {
        clear();
        //make sure that scheduler wakes up and terminates
        active_.store(false);
        //wait for scheduler to terminate
        scheduler_.join();
    }


    //---------------------------------------------------------------
    parallel_task_queue& operator = (const parallel_task_queue&) = delete;
    parallel_task_queue& operator = (parallel_task_queue&&) = delete;


    //---------------------------------------------------------------
    void
    enqueue(const task_type& t)
    {
        std::lock_guard<std::mutex> lock(mutables_);
        waiting_.push_back(t);
        hasWaiting_.store(true);
    }
    //-----------------------------------------------------
    void
    enqueue(task_type&& t)
    {
        std::lock_guard<std::mutex> lock(mutables_);
        waiting_.push_back(std::move(t));
        hasWaiting_.store(true);
    }
    //-----------------------------------------------------
    template<class InputIterator>
    void
    enqueue(InputIterator first, InputIterator last)
    {
        std::lock_guard<std::mutex> lock(mutables_);
        waiting_.insert(begin(waiting_), first, last);
        hasWaiting_.store(true);
    }
    //-----------------------------------------------------
    void
    enqueue(std::initializer_list<task_type> il)
    {
        std::lock_guard<std::mutex> lock(mutables_);
        waiting_.insert(waiting_.end(), il);
        hasWaiting_.store(true);
    }
    //-----------------------------------------------------
    template<class... Args>
    void
    enqueue_emplace(Args&&... args)
    {
        std::lock_guard<std::mutex> lock(mutables_);
        waiting_.emplace_back(std::forward<Args>(args)...);
        hasWaiting_.store(true);
    }


    /****************************************************************
     * @brief removes first waiting task that compares equal to 't'
     */
    bool
    try_remove(const task_type& t)
    {
        std::lock_guard<std::mutex> lock(mutables_);
        auto it = std::find(waiting_.begin(), waiting_.end(), t);
        if(it != waiting_.end()) {
            waiting_.erase(it);
            if(waiting_.empty()) hasWaiting_.store(false);
            return true;
        }
        return false;
    }


    //---------------------------------------------------------------
    void
    clear() {
        std::lock_guard<std::mutex> lock(mutables_);
        waiting_.clear();
        hasWaiting_.store(false);
    }


    //---------------------------------------------------------------
    unsigned int
    concurrency() const noexcept {
        return workers_.size();
    }

    //---------------------------------------------------------------
    bool
    empty() const noexcept {
        return !hasWaiting_.load();
    }


    //-----------------------------------------------------
    /// @return number of waiting tasks
    std::size_t
    waiting() const noexcept {
        std::lock_guard<std::mutex> lock(mutables_);
        return waiting_.size();
    }
    std::size_t
    unsafe_waiting() const noexcept {
        return waiting_.size();
    }

    //-----------------------------------------------------
    /// @return number of running tasks
    std::size_t
    running() const noexcept {
        return running_.load();
    }

    //-----------------------------------------------------
    /// @return true, if all threads are working on a task
    bool
    busy() const noexcept {
        return running_.load() >= int(workers_.size());
    }

    //-----------------------------------------------------
    bool
    complete() const noexcept {
        std::lock_guard<std::mutex> lock(mutables_);
        return empty() && (running() < 1);
    }


    //---------------------------------------------------------------
    /**
     * @brief block execution of calling thread until all tasks are complete
     */
    void wait()
    {
        std::unique_lock<std::mutex> lock{mutables_};
        while(!empty() || running()) isDone_.wait(lock);
    }


private:

    //---------------------------------------------------------------
    void try_assign_tasks()
    {
        for(auto& worker : workers_) {
            if(worker.available()) {
                std::lock_guard<std::mutex> lock{mutables_};
                if(waiting_.empty()) {
                    hasWaiting_.store(false);
                    return;
                }
                if(worker(this, std::move(waiting_.front()))) {
                    ++running_;
                    waiting_.pop_front();
                    if(waiting_.empty()) {
                        hasWaiting_.store(false);
                        return;
                    }
                    else if(busy()) {
                        return;
                    }
                }
            }
        }
    }

    //---------------------------------------------------------------
    /// @brief this will run in a separate, dedicated thread
    void schedule()
    {
         while(active_.load()) {
             if(!empty() && !busy()) {
                 try_assign_tasks();
             }
             else if(running() < 1) {
                 std::lock_guard<std::mutex> lock{mutables_};
                 if(empty() && (running() < 1)) {
                     isDone_.notify_all();
                 }
             }
         }
    }


    //---------------------------------------------------------------
    mutable std::mutex mutables_;
    std::atomic_bool active_;
    std::atomic_bool hasWaiting_;
    std::atomic_int running_;
    std::deque<task_type> waiting_;
    std::vector<task_thread<task_executor>> workers_;
    std::condition_variable isDone_;
    std::thread scheduler_;
};



//-------------------------------------------------------------------
/// @brief convenience alias
using parallel_function_queue = parallel_task_queue<std::function<void()>>;



} // namespace mc


#endif

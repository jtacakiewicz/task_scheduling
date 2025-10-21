#ifndef EMP_THREAD_POOL_HPP
#define EMP_THREAD_POOL_HPP
#include <condition_variable>
#include <functional>
#include <queue>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>

struct Worker;
class TaskQueue {
    std::queue< std::function<void()> > m_tasks;
    std::mutex m_mutex;
    std::atomic<size_t> m_remaining_tasks = 0;
    std::condition_variable m_cond;

    void getTask(std::function<void()>& target_callback) {
        std::lock_guard<std::mutex> lock_guard{m_mutex};
        if (m_tasks.empty()) {
            return;
        }
        target_callback = std::move(m_tasks.front());
        m_tasks.pop();
    }
    void wait() {
        std::unique_lock<std::mutex> lk(m_mutex);
        m_cond.wait(lk);
    }
    void wakeUpAll() {
        m_mutex.unlock();
        m_cond.notify_all();
    }
    void workDone() {
        m_remaining_tasks--;
    }
public:
    template<typename TCallback>
    void addTask(TCallback&& callback) {
        std::lock_guard<std::mutex> lock_guard{m_mutex};
        m_tasks.push(std::forward<TCallback>(callback));
        m_cond.notify_all();
        m_remaining_tasks++;
    }

    void waitForCompletion() const {
        while (m_remaining_tasks > 0) {
            std::this_thread::yield();
        }
    }
    friend Worker;
    ~TaskQueue() {
    }
};
class Worker
{
    uint32_t                m_id      = 0;
    std::thread             m_thread;
    std::function<void()>   m_task    = nullptr;
    std::atomic<bool>       m_running = true;
    TaskQueue*              m_queue   = nullptr;
public:

    Worker() = default;

    Worker(TaskQueue& queue, uint32_t id)
        : m_id{id}
        , m_queue{&queue}
    {
        m_thread = std::thread([this](){
            run();
        });
    }

    void run() {
        while (m_running) {
            m_queue->getTask(m_task);
            if (m_task == nullptr) {
                m_queue->wait();
            } else {
                m_task();
                m_queue->workDone();
                m_task = nullptr;
            }
        }
    }
    void stop() {
        m_running = false;
        m_queue->wakeUpAll();
        m_thread.join();
    }
};
class ThreadPool {
    uint32_t            m_thread_count = 0;
    TaskQueue           m_queue;
    std::vector<std::unique_ptr<Worker>> m_workers;
public:

    explicit ThreadPool(uint32_t thread_count = std::thread::hardware_concurrency())
        : m_thread_count{thread_count}
    {
        m_workers.reserve(thread_count);
        for (uint32_t i{thread_count}; i--;) {
            m_workers.emplace_back(std::make_unique<Worker>(m_queue, static_cast<uint32_t>(m_workers.size())));
        }
    }

    ~ThreadPool()
    {
        for (auto& worker : m_workers) {
            worker.get()->stop();
        }
    }

    template<typename TCallback>
    void addTask(TCallback&& callback)
    {
        m_queue.addTask(std::forward<TCallback>(callback));
    }

    void waitForCompletion() const
    {
        m_queue.waitForCompletion();
    }

    template<typename TCallback>
    void dispatch(uint32_t element_count, TCallback&& callback)
    {
        const uint32_t batch_size = element_count / m_thread_count;
        for (uint32_t i{0}; i < m_thread_count; ++i) {
            const uint32_t start = batch_size * i;
            const uint32_t end   = start + batch_size;
            addTask([start, end, &callback](){ callback(start, end); });
        }

        if (batch_size * m_thread_count < element_count) {
            const uint32_t start = batch_size * m_thread_count;
            callback(start, element_count);
        }
    }
};
#endif // EMP_THREAD_POOL_HPP

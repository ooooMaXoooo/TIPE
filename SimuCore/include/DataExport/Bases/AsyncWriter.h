#pragma once

#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <memory>
#include <filesystem>
// #include <vector>
// #include <string>

#include <DataExport/Bases/Writable.h>


class AsyncWriter {
public:
    explicit AsyncWriter(const std::string_view& directory);

    ~AsyncWriter();

    virtual void enqueue(std::shared_ptr<Writable> data, const std::string filename);

protected:
    std::filesystem::path m_directory;
    std::queue<std::pair<std::shared_ptr<Writable>, std::string>> m_queue;
    std::mutex m_mutex;
    std::condition_variable m_cv;
    std::atomic<bool> m_stop;
    std::thread m_worker;

    virtual void process();
};

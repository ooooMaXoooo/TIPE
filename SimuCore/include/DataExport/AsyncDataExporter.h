#pragma once

#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <memory>
#include <filesystem>
#include <string>
#include <stdexcept>
#include <fstream>

#include "Writable.h"


class AsyncDataExporter {
public :
    explicit AsyncDataExporter();

    ~AsyncDataExporter();

    void enqueue(std::shared_ptr<Writable> data, std::filesystem::path path);

protected:
    using data_type = std::pair<std::shared_ptr<Writable>, std::string>;
    std::queue<data_type> m_queue;

    std::mutex m_mutex;
    std::condition_variable m_cv;
    std::atomic<bool> m_stop;
    std::thread m_worker;

    void process();
    
    void writeData(std::ofstream& file, const Writable& data);
};
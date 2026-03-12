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

class Writable;
class AsyncDataExporter;



class Writable {
public:
    virtual std::string string() const {
        throw std::runtime_error("The function Writable::GetData needs to be override !");
    }
};

std::ofstream& operator<<(std::ofstream& fs, const Writable& writable);


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
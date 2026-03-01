#include <pch.h>
#include <DataExport/Bases/AsyncWriter.h>

AsyncWriter::AsyncWriter(const std::string_view &directory)
    : m_stop(false)
{
    m_directory = directory;
    if (!std::filesystem::is_directory(m_directory)) {
        throw std::invalid_argument("directory is not a directory");
    }
    std::filesystem::create_directories(m_directory);

    m_worker = std::thread([this] { this->process(); });
}

AsyncWriter::~AsyncWriter() {
    {
        std::unique_lock lock(m_mutex);
        m_stop = true;
    }
    m_cv.notify_one();
    if (m_worker.joinable())
        m_worker.join();
}

void AsyncWriter::enqueue(std::shared_ptr<Writable> data,const std::string filename){
    {
        std::unique_lock lock(m_mutex);
        m_queue.push({data, filename});
    }
    m_cv.notify_one();
}

void AsyncWriter::process(){
    throw std::runtime_error("AsyncWriter shouldn't be used, it is a base class");
}
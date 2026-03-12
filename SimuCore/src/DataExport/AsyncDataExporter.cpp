#include <pch.h>
#include <DataExport/AsyncDataExporter.h>

AsyncDataExporter::AsyncDataExporter()
    : m_stop(false)
{
    m_worker = std::thread([this] { this->process(); });
}

AsyncDataExporter::~AsyncDataExporter() {
    {
        std::unique_lock lock(m_mutex);
        m_stop = true;
    }
    m_cv.notify_one();
    if (m_worker.joinable())
        m_worker.join();
}

void AsyncDataExporter::enqueue(std::shared_ptr<Writable> data, std::filesystem::path path) {
    std::filesystem::create_directories(path.parent_path());

    {
        std::unique_lock lock(m_mutex);
        m_queue.push({data, path.string()});
    }
    m_cv.notify_one();
}

void AsyncDataExporter::process() {
    std::shared_ptr<std::ofstream> file = nullptr;

    while (true) {
        std::shared_ptr<Writable> data;
        {
            std::unique_lock lock(m_mutex);
            m_cv.wait(lock, [this] { return !m_queue.empty() || m_stop; });

            if (m_stop && m_queue.empty())
                break;

            data = m_queue.front().first;
            file = std::make_shared<std::ofstream>(m_queue.front().second);

            if (!file->is_open()) {
                throw std::runtime_error("Impossible d'ouvrir le fichier demande !");
            }

            m_queue.pop();
        }

        if (data) {
            writeData(*file, *data);

            file->close();
            file.reset();
        }
    }
}

void AsyncDataExporter::writeData(std::ofstream& file, const Writable& data) {
    file << data;
}

std::ofstream& operator<<(std::ofstream& fs, const Writable& writable)
{
    fs << writable.string();
    return fs;
}

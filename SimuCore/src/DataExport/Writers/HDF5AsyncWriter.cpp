#include <pch.h>
#include <DataExport/Writers/HDF5AsyncWriter.h>

HDF5AsyncWriter::HDF5AsyncWriter(const std::string_view &directory)
    : AsyncWriter(directory)
{}

void HDF5AsyncWriter::process() {
    std::unique_ptr<H5::H5File> file_ptr;

    while (true) {
        std::shared_ptr<Writable> data;
        {
            std::unique_lock lock(m_mutex);
            m_cv.wait(lock, [this] { return !m_queue.empty() || m_stop; });

            if (m_stop && m_queue.empty())
                break;

            data = m_queue.front().first;

            file_ptr = std::make_unique<H5::H5File>(m_queue.front().second, H5F_ACC_TRUNC);
            m_queue.pop();
        }

        if (data) {
            write_data(*file_ptr, *data);
            file_ptr->close();
            file_ptr.reset();
        }
    }
}

void HDF5AsyncWriter::write_data(H5::H5File &file, const Writable &data) {
    throw std::runtime_error("Implement a specific method in a child class");
}
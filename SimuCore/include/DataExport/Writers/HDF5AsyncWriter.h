#pragma once

#include <H5Cpp.h>

#include <DataExport/Bases/AsyncWriter.h>
#include <DataExport/Bases/Writable.h>

#include <string_view>

class HDF5AsyncWriter : public AsyncWriter {
    public:
    HDF5AsyncWriter(const std::string_view& directory);

    protected :

    void process() override;

    virtual void write_data(H5::H5File& file, const Writable& data);
};
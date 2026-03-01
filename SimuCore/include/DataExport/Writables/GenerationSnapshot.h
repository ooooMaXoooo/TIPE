#pragma once

#include <DataExport/Bases/Writable.h>

class GenerationSnapshot : public Writable {

    virtual std::string GetData() const override {
        throw std::runtime_error("The function GenerationSnapshot::GetData needs to be created !");
    }


};
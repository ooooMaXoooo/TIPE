#pragma once

#include <string>
#include <stdexcept>

class Writable {
    public :
    Writable() {}

    virtual std::string GetData() const {
        throw std::runtime_error("The function Writable::GetData needs to be override !");
    }
};
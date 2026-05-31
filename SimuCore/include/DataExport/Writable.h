#pragma once

#include <string>
#include <fstream>

class Writable {
public:
    virtual std::string string() const = 0;
};

std::ofstream& operator<<(std::ofstream& fs, const Writable& writable);
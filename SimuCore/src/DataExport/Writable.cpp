#include <pch.h>
#include <DataExport/Writable.h>

std::ofstream& operator<<(std::ofstream& fs, const Writable& writable)
{
    fs << writable.string();
    return fs;
}
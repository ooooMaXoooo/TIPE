#pragma once

#include <DataExport/Writable.h>

#include <vector>

class ClustersData : public Writable {
public:
	ClustersData(const std::vector<std::vector<int>>& clusters);

	std::string string() const override;

private:
	std::vector<std::vector<int>> m_clusters;
};
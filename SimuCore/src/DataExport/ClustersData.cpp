#include <pch.h>
#include <DataExport/ClustersData.h>

ClustersData::ClustersData(const std::vector<std::vector<int>>& clusters) :
	m_clusters(clusters)
{}

std::string ClustersData::string() const {
	std::ostringstream oss;

	oss << m_clusters.size() << '\n';

	for (const std::vector<int>& cluster : m_clusters) {
		int i = 0;
		int n = cluster.size() - 1;
		for (; i < n; i++) {
			oss << cluster[i] << ';';
		}
		oss << cluster[i] << '\n';
	}

	return oss.str();
}
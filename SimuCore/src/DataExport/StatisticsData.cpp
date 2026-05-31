#include <pch.h>
#include <DataExport/StatisticsData.h>

StatisticsData::StatisticsData(size_t NClusters, double best_score, double worst_score, double mean_score, int* kindsCount) :
	m_NClusters(NClusters),
	m_best_score(best_score),
	m_worst_score(worst_score),
	m_mean_score(mean_score),
	m_kindsCount(kindsCount)
{}

std::string StatisticsData::string() const {
	std::ostringstream oss;
	oss << std::setprecision(std::numeric_limits<double>::max_digits10);

	oss << m_NClusters << '\n'
		<< m_best_score << '\n'
		<< m_worst_score << '\n'
		<< m_mean_score << '\n';

	for (int i = 0; i < s_kindsCountSize; ++i) {
		oss << m_kindsCount[i] << '\n';
	}

	return oss.str();
}

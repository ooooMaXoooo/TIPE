#include <pch.h>
#include <DataExport/StatisticsData.h>

StatisticsData::StatisticsData(size_t NClusters, double best_score, double worst_score, double mean_score, std::array<int, 12> kindsCount) :
	m_NClusters(NClusters),
	m_best_score(best_score),
	m_worst_score(worst_score),
	m_mean_score(mean_score) {

	#pragma omp critical
	{
		for (int i = 0; i < 12; i++) {
			m_kindsCount[i] = kindsCount[i];
		}
	}
}


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

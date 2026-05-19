#pragma once

#include <DataExport/AsyncDataExporter.h>

class StatisticsData : public Writable {
public :
	StatisticsData(
		size_t NClusters,
		double best_score,
		double worst_score,
		double mean_score,
		int* kindsCount
	);

	std::string string() const override;

private:
	/** On veut :
	* - le nombre de patates (clusters) de la population
	* - meilleur score de la population
	* - pire score de la population
	* - score moyen de la population
	* - le nombre de fusée ayant un même état, pour chaque état, i.e la taille de chaque état
	*/

	size_t m_NClusters;

	double m_best_score;
	double m_worst_score;
	double m_mean_score;

	int* m_kindsCount;
	static constexpr int s_kindsCountSize = 12;

};
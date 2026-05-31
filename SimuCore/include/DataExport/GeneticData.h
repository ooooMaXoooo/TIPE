#pragma once

#include <DataExport/Writable.h>

class GeneticData : public Writable {
	public:
	GeneticData() = default;

	std::string string() const override;

private :
	size_t m_generationId;

	double m_bestFitness;
	// double m_averageFitness; // pas encore implémenté

	size_t m_populationSize;

};
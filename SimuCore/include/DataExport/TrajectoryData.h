#pragma once

#include <DataExport/Writable.h>
#include <string>
#include <glm/glm.hpp>

class TrajectoryData : public Writable {
public:
	TrajectoryData(
		const std::vector<glm::dvec3>& trajectory,
		int dimension = 2);

	std::string string() const override;

private:
	std::vector<glm::dvec3> m_trajectory;
	int m_dimension;

private:
	std::string formatVector(const glm::dvec3& vec) const;
};
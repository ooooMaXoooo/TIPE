#pragma once

#include <DataExport/AsyncDataExporter.h>
#include <pch.h>

class TrajectoryData : public Writable {
public:
	TrajectoryData(
		const std::vector<glm::dvec3>& trajectory,
		int dimension = 2)
	: m_trajectory(trajectory),
		m_dimension(dimension)
	{
		if (dimension != 2 && dimension != 3) {
			throw std::invalid_argument("Dimension must be 2 or 3");
		}
	}

	std::string string() const override {
		std::ostringstream oss;
		for (const glm::dvec3& vec : m_trajectory) {
			oss << formatVector(vec) << '\n';
		}
		return oss.str();
	}

private:
	std::vector<glm::dvec3> m_trajectory;
	int m_dimension;

private:
	std::string formatVector(const glm::dvec3& vec) const {
		std::ostringstream oss;
		if (m_dimension == 2) {
			oss << vec.x << ';' << vec.y;
			return oss.str();
		}
		else
		{
			oss << vec.x << ';' << vec.y << ';' << vec.z;
		}
		return oss.str();
	}
};
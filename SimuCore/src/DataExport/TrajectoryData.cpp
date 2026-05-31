#include <pch.h>
#include <DataExport/TrajectoryData.h>

TrajectoryData::TrajectoryData(const std::vector<glm::dvec3>& trajectory, int dimension) :
	m_trajectory(trajectory),
	m_dimension(dimension)
{
	if (dimension != 2 && dimension != 3) {
		throw std::invalid_argument("Dimension must be 2 or 3");
	}
}

std::string TrajectoryData::string() const {
	std::ostringstream oss;
	oss << std::setprecision(std::numeric_limits<double>::max_digits10);

	for (const glm::dvec3& vec : m_trajectory) {
		oss << formatVector(vec) << '\n';
	}
	return oss.str();
}

std::string TrajectoryData::formatVector(const glm::dvec3& vec) const {
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

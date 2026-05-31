#pragma once

#include <DataExport/Writable.h>
#include <string>
#include <cstdint>
#include <limits>
#include <utility>
#include <glm/fwd.hpp>


class RocketData : public Writable {
public:
	RocketData(
		size_t id,
		double tof,
		double maxTime,
		glm::dvec3 initialImpulsion,
		glm::dvec3 initialPosition,
		glm::dvec3 initialVelocity,
		size_t finalPlanet_startIndice,
		uint8_t startPlanet,
		uint8_t finalPlanet,
		uint8_t nbImpulsions,
		std::vector<std::vector<uint32_t>> genome,
		uint8_t dimension = 2
	);

	~RocketData();

	std::string string() const override;

	void RegisterNewCartesianPosition(double x, double y);
	void RegisterNewCartesianPosition(glm::dvec3 pos);
	void RegisterNewPolarPosition(double R, double theta);
	void RegisterNewImpulsion(const glm::dvec3& impulsion, double time);

private:
	uint8_t m_dimension;

	size_t m_id;

	// en rad
	double m_thetaMin = std::numeric_limits<double>::max();
	// en rad
	double m_thetaMax = std::numeric_limits<double>::lowest();
	// en rad
	double m_thetaMean = 0.0;

	// en rad
	double m_TotalTheta = 0.0;
	// en rad
	double m_LastTheta = 0.0;

	// en UA
	double m_Rmin = std::numeric_limits<double>::max();
	// en UA
	double m_Rmax = std::numeric_limits<double>::lowest();
	// en UA
	double m_Rmean = 0.0;

	// en jours
	double m_tof;
	// en jours
	double m_maxTime = 0.0;

	// en km/s
	glm::dvec3 m_initialImpulsion;
	// en (km/s, jours)
	std::vector<std::pair<glm::dvec3, double>> m_effectivesImpulsions;

	// en UA
	glm::dvec3 m_initialPosition;
	// en km/s
	glm::dvec3 m_initialVelocity;

	size_t m_finalPlanet_startIndice;

	uint8_t m_startPlanet;
	uint8_t m_finalPlanet;

	uint8_t m_nbImpulsions;
	uint8_t m_currentImpulsionIndex = 0;

	size_t m_NbRegisteredPositions = 0;


	std::vector<std::vector<uint32_t>> m_genome; // pour stocker le génome de l'individu génétique associé à cette fusée, pour pouvoir le simuler et avoir sa trajectoire en python

private:
	std::string formatVector(const glm::dvec3& vec) const;
	std::string formatImpulsions() const;
	std::string formatGenome() const;

	void RegisterNewTheta(double theta);
	void RegisterNewR(double R);

public:
	inline size_t GetId() const noexcept { return m_id; }
	inline void SetTof(double tof) noexcept { m_tof = tof; }

	friend double distance(const RocketData& r1, const RocketData& r2);
};
#include <pch.h>
#include <DataExport/RocketData.h>

#include <SimuCore/constants.h>
#include <cmath>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <glm/fwd.hpp>

RocketData::RocketData(size_t id, double tof, double maxTime, glm::dvec3 initialImpulsion, glm::dvec3 initialPosition, glm::dvec3 initialVelocity, size_t finalPlanet_startIndice, uint8_t startPlanet, uint8_t finalPlanet, uint8_t nbImpulsions, std::vector<std::vector<uint32_t>> genome, uint8_t dimension)
	:
	m_dimension(dimension),
	m_id(id),
	m_tof(tof),
	m_maxTime(maxTime),
	m_initialImpulsion(initialImpulsion),
	m_initialPosition(initialPosition),
	m_initialVelocity(initialVelocity),
	m_finalPlanet_startIndice(finalPlanet_startIndice),
	m_startPlanet(startPlanet),
	m_finalPlanet(finalPlanet),
	m_nbImpulsions(nbImpulsions),
	m_genome(genome)
{
	assert(dimension == 2 && "We need polar coordinates to use theta and R, so dimension must be 2");
	assert(m_nbImpulsions > 1 && "There must be at least 2 impulsions to store the initial impulsion and at least one effective impulsion");

	m_effectivesImpulsions.resize(nbImpulsions - 1); // la première impulsion étant stockée en dehors du tableau, on n'a besoin que de nbImpulsions-1 cases pour stocker les impulsions suivantes.
	for (size_t i = 0; i < nbImpulsions-1; i++) {
		m_effectivesImpulsions[i] = std::make_pair(glm::dvec3(0.0), 0.0);
	}
}

RocketData::~RocketData() {
}

std::string RocketData::string() const {
	std::ostringstream oss;

	//oss << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
	oss << std::setprecision(std::numeric_limits<double>::max_digits10);

	oss << m_id << '\n'
		<< m_tof << '\n'
		<< m_maxTime << '\n'
		<< m_thetaMin << '\n' << m_thetaMax << '\n' << m_thetaMean / m_NbRegisteredPositions << '\n'
		<< m_TotalTheta / (2*SimuCore::constants::PI) << '\n'
		<< m_Rmin << '\n' << m_Rmax << '\n' << m_Rmean / m_NbRegisteredPositions << '\n'
		<< static_cast<int>(m_startPlanet) << '\n' << static_cast<int>(m_finalPlanet) << '\n'
		<< formatVector(m_initialPosition) << '\n'
		<< formatVector(m_initialVelocity) << '\n'
		<< m_finalPlanet_startIndice << '\n'
		<< formatVector(m_initialImpulsion) << '\n'
		<< static_cast<int>(m_nbImpulsions) << '\n'
		<< formatImpulsions() << '\n'
		<< formatGenome();

	return oss.str();
}

void RocketData::RegisterNewCartesianPosition(double x, double y) {
	double norm = std::sqrt(x * x + y * y);
	double theta = std::atan2(y, x);
	RegisterNewPolarPosition(norm, theta);
}

void RocketData::RegisterNewCartesianPosition(glm::dvec3 pos) {
	RegisterNewCartesianPosition(pos.x, pos.y);
}

void RocketData::RegisterNewPolarPosition(double R, double theta) {
	RegisterNewTheta(theta);
	RegisterNewR(R);

	m_NbRegisteredPositions++;
}

void RocketData::RegisterNewImpulsion(const glm::dvec3& impulsion, double time) {
	if (m_currentImpulsionIndex < m_nbImpulsions) {
		m_effectivesImpulsions[m_currentImpulsionIndex].first = impulsion;
		m_effectivesImpulsions[m_currentImpulsionIndex].second = time;
		m_currentImpulsionIndex++;
	}
	else {
		std::cerr << "Error: trying to register more impulsions than the specified number of impulsions" << std::endl;
		throw std::runtime_error("Too many impulsions registered");
	}
}

std::string RocketData::formatVector(const glm::dvec3& vec) const {
	std::ostringstream oss;
	oss << vec.x << ';' << vec.y;
	return oss.str();
}

std::string RocketData::formatImpulsions() const {
	std::ostringstream oss;
	for (size_t i = 0; i < m_nbImpulsions-1; i++) {
		oss << formatVector(m_effectivesImpulsions[i].first) << '\n' << m_effectivesImpulsions[i].second;
		if (i < m_nbImpulsions - 1) {
			oss << '\n';
		}
	}
	return oss.str();
}

std::string RocketData::formatGenome() const {
	std::ostringstream oss;
	
	// on affiche le nombre de chromosomes et le nombre de gènes par chromosome pour simplifier la lecture des données
	oss << m_genome.size() << '\n'
		<< (m_genome.empty() ? 0 : m_genome[0].size()) << '\n';

	for (size_t i = 0; i < m_genome.size(); i++) {
		for (size_t j = 0; j < m_genome[i].size(); j++) {
			oss << m_genome[i][j];
			if (j < m_genome[i].size() - 1) {
				oss << ';';
			}
		}
		if (i < m_genome.size() - 1) {
			oss << '\n';
		}
	}

	return oss.str();
}

void RocketData::RegisterNewTheta(double theta) {
	if (theta < m_thetaMin) {
		m_thetaMin = theta;
	}
	if (theta > m_thetaMax) {
		m_thetaMax = theta;
	}
	m_thetaMean += theta;

	m_TotalTheta += theta - m_LastTheta;
	m_LastTheta = theta;
}

void RocketData::RegisterNewR(double R) {
	if (R < m_Rmin) {
		m_Rmin = R;
	}
	if (R > m_Rmax) {
		m_Rmax = R;
	}
	m_Rmean += R;
}

double distance(const RocketData& r1, const RocketData& r2) {
	DistInfoRocketData info1, info2;

	constexpr double PI2 = 2 * SimuCore::constants::PI;

	info1.rMin		= r1.m_Rmin;
	info1.rMax		= r1.m_Rmax;
	info1.rMean		= r1.m_Rmean / r1.m_NbRegisteredPositions;
	info1.thetaMin	= r1.m_thetaMin;
	info1.thetaMax	= r1.m_thetaMax;
	info1.thetaMean	= r1.m_thetaMean / r1.m_NbRegisteredPositions;
	info1.nbImpuls	= r1.m_currentImpulsionIndex;
	info1.nbTurns	= r1.m_TotalTheta / PI2;

	info2.rMin		= r2.m_Rmin;
	info2.rMax		= r2.m_Rmax;
	info2.rMean		= r2.m_Rmean / r2.m_NbRegisteredPositions;
	info2.thetaMin	= r2.m_thetaMin;
	info2.thetaMax	= r2.m_thetaMax;
	info2.thetaMean = r2.m_thetaMean / r2.m_NbRegisteredPositions;
	info2.nbImpuls	= r2.m_currentImpulsionIndex;
	info2.nbTurns	= r2.m_TotalTheta / PI2;

	return distance(info1, info2);
}

double distance(const DistInfoRocketData& r1, const DistInfoRocketData& r2) {
	constexpr double coefs[7] = { 100, 250, 25, 100, 25, 500, 10 };

	double dist_rMin = 0, dist_rMax = 0, dist_rMean = 0, dist_theta = 0, dist_thetaMean = 0, dist_turns = 0, dist_nbImpulsions = 0;

	dist_nbImpulsions = std::abs(r1.nbImpuls - r2.nbImpuls);

	// les rayons devraient être en UA

	dist_rMin = std::abs(r1.rMin - r2.rMin);
	dist_rMax = std::abs(r1.rMax - r2.rMax);
	dist_rMean = std::abs(r1.rMean - r2.rMean);

	// les angles en radians
	dist_thetaMean = std::abs(r1.thetaMean - r2.thetaMean);
	dist_theta = std::abs(std::abs(r1.thetaMax - r1.thetaMin) - std::abs(r2.thetaMax - r2.thetaMin));

	dist_turns = std::abs(r1.nbTurns - r2.nbTurns);

	double distance = coefs[0] * dist_rMin + coefs[1] * dist_rMax + coefs[2] * dist_rMean + coefs[3] * dist_theta + coefs[4] * dist_thetaMean + coefs[5] * dist_turns + coefs[6] * dist_nbImpulsions;

	return distance;
}
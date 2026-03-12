#pragma once

#include <DataExport/AsyncDataExporter.h>
#include <pch.h>

class PhysicsData : public Writable {
public:
	PhysicsData(
		double simulation_time,
		double time_step,
		uint8_t start_planet_index,
		uint8_t final_planet_index,
		const glm::dvec3& start_planet_start_position,
		const glm::dvec3& final_planet_start_position,
		const glm::dvec3& start_planet_final_position,
		const glm::dvec3& final_planet_final_position,
		const glm::dvec3& final_rocket_velocity,
		const std::vector<double>& impulse_times,
		const std::vector<glm::dvec3>& impulse_vectors,
		int dimension = 2) :
		m_dimension(dimension),
		m_simulation_time(simulation_time),
		m_time_step(time_step),
		m_start_planet_index(start_planet_index),
		m_final_planet_index(final_planet_index),
		m_start_planet_start_position(start_planet_start_position),
		m_final_planet_start_position(final_planet_start_position),
		m_start_planet_final_position(start_planet_final_position),
		m_final_planet_final_position(final_planet_final_position),
		m_final_rocket_velocity(final_rocket_velocity),
		m_impulse_times(impulse_times),
		m_impulse_vectors(impulse_vectors)
	{
		if (dimension != 2 && dimension != 3) {
			throw std::invalid_argument("Dimension must be 2 or 3");
		}

		m_number_of_impulsions = impulse_times.size();
	}

	std::string string() const override {
		std::ostringstream oss;

		/* On veut afficher :
		*	- Position de départ de la planète de départ
		*	- Position de départ de la planète d'arrivée
		*	
		*	- Position d'arrivée de la planète de départ
		*	- Position d'arrivée de la planète d'arrivée
		* 
		*	- le dernier vecteur vitesse de la fusée
		* 
		*	- le nombre d'impulsions
		*	- les instants des impulsions
		*	- les vecteurs d'impulsion
		* 
		*	- temps de simulation
		*	- pas de temps
		* 
		*	- numéro de la planète de départ
		*	- numéro de la planète d'arrivée 
		*/

		oss << m_simulation_time << '\n' 
			<< m_time_step << '\n' 
			<< static_cast<int>(m_start_planet_index) << '\n' 
			<< static_cast<int>(m_final_planet_index) << '\n' 
			<< m_number_of_impulsions << '\n' 
			<< formatVector(m_start_planet_start_position) << '\n'
			<< formatVector(m_final_planet_start_position) << '\n'
			<< formatVector(m_start_planet_final_position) << '\n'
			<< formatVector(m_final_planet_final_position) << '\n'
			<< formatVector(m_final_rocket_velocity) << '\n';

		for (size_t i = 0; i < m_number_of_impulsions; ++i) {
				oss << m_impulse_times[i] << '\n';
				oss << formatVector(m_impulse_vectors[i]) << '\n';
		}

		return oss.str();
	}

private:
	double m_dimension;

	double m_simulation_time;
	double m_time_step;

	uint8_t m_start_planet_index;
	uint8_t m_final_planet_index;

	size_t m_number_of_impulsions;

	glm::dvec3 m_start_planet_start_position;
	glm::dvec3 m_final_planet_start_position;

	glm::dvec3 m_start_planet_final_position;
	glm::dvec3 m_final_planet_final_position;

	glm::dvec3 m_final_rocket_velocity;

	std::vector<double> m_impulse_times;
	std::vector<glm::dvec3> m_impulse_vectors;


private :
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
#pragma once

#include <DataExport/Writable.h>
#include <string>
#include <glm/glm.hpp>

class PhysicsData : public Writable {
public:

	/// <summary>
	/// Classe pour empacter les données physiques utiles à une simulation
	/// </summary>
	/// <param name="simulation_time"> en jours </param>
	/// <param name="time_step"> en secondes </param>
	/// <param name="start_planet_index"> uint8_t </param>
	/// <param name="final_planet_index"> uint8_t </param>
	/// <param name="start_planet_start_position"> en UA </param>
	/// <param name="final_planet_start_position"> en UA </param>
	/// <param name="start_planet_final_position"> en UA </param>
	/// <param name="final_planet_final_position"> en UA </param>
	/// <param name="final_rocket_velocity"> en km/s </param>
	/// <param name="impulse_times"> jours </param>
	/// <param name="impulse_vectors"> km/s </param>
	/// <param name="tof"> en jours </param>
	/// <param name="delta_v"> en km/s </param>
	/// <param name="etat_lie"> bool </param>
	/// <param name="r_min"> en km </param>
	/// <param name="r_max"> en km </param>
	/// <param name="dimension"> uint16_t </param>
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
		double tof,
		double delta_v,
		bool etat_lie = false,
		double r_min = -1,
		double r_max = -1,
		uint16_t dimension = 2);

	std::string string() const override;

private:
	double m_dimension;

	// en jours
	double m_simulation_time;
	// en secondes
	double m_time_step;

	
	uint8_t m_start_planet_index;
	uint8_t m_final_planet_index;

	size_t m_number_of_impulsions;

	// en UA
	glm::dvec3 m_start_planet_start_position;
	//en UA
	glm::dvec3 m_final_planet_start_position;

	// en UA
	glm::dvec3 m_start_planet_final_position;
	// en UA
	glm::dvec3 m_final_planet_final_position;

	// en km/s
	glm::dvec3 m_final_rocket_velocity;

	// en jours
	std::vector<double> m_impulse_times;

	// en km/s
	std::vector<glm::dvec3> m_impulse_vectors;

	bool m_etat_lie;
	double m_r_min;
	double m_r_max;

	double m_tof;
	double m_delta_v;

private:

	std::string formatVector(const glm::dvec3& vec) const;
};
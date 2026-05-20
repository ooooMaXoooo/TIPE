#pragma once

#include <pch.h>
#include <SimuCore/structures/structures.h>
#include <GaLib/genetic.hpp>

#include "Patate.h"

namespace SimuCore {
	namespace Statistics {

		template <size_t NbImpulsions>
		using ConfigType = genetic::Config<double, uint32_t, 2 * NbImpulsions + 1, 2>;

		template <size_t NbImpulsions>
		using Ind = typename genetic::Individu<ConfigType<NbImpulsions>>;

		template <size_t NbImpulsions>
		using Pop = typename std::vector<Ind<NbImpulsions>>;

		void rotate_rocket(SimuCore::Structures::Rocket* rocket, SimuCore::Systems::AdaptedSystem& system) {

			glm::dvec3 pos_terre = system.GetStartPlanet_StartPosition();
			double angle_rad = - std::arg(std::complex<double>(pos_terre.x, pos_terre.y)); // Angle en radians
			
			// Rotation autour de l'axe z
			double cos_angle = std::cos(angle_rad);
			double sin_angle = std::sin(angle_rad);

			glm::dvec3 new_position;
			{
				new_position.x = rocket->position.x * cos_angle - rocket->position.y * sin_angle;
				new_position.y = rocket->position.x * sin_angle + rocket->position.y * cos_angle;
				new_position.z = rocket->position.z; // Pas de changement pour la composante z
			}
			rocket->position = new_position;

			glm::dvec3 new_velocity;
			{
				new_velocity.x = rocket->velocity.x * cos_angle - rocket->velocity.y * sin_angle;
				new_velocity.y = rocket->velocity.x * sin_angle + rocket->velocity.y * cos_angle;
				new_velocity.z = rocket->velocity.z; // Pas de changement pour la composante z
			}
			rocket->velocity = new_velocity;

			std::vector<std::pair<SimuCore::Structures::Impulsion, double>> new_impulsions;
			std::vector<std::pair<SimuCore::Structures::Impulsion, double>> impulsions = rocket->getImpulsions();
			new_impulsions.reserve(impulsions.size());

			for (auto& [impulsion, time] : impulsions) {
				glm::dvec3 new_deltaV;
				{
					new_deltaV.x = impulsion.GetDeltaV_vec().x * cos_angle - impulsion.GetDeltaV_vec().y * sin_angle;
					new_deltaV.y = impulsion.GetDeltaV_vec().x * sin_angle + impulsion.GetDeltaV_vec().y * cos_angle;
					new_deltaV.z = impulsion.GetDeltaV_vec().z; // Pas de changement pour la composante z
				}

				new_impulsions.emplace_back(SimuCore::Structures::Impulsion(new_deltaV), time);
			}

			rocket->setImpulsions(new_impulsions);
		}

		template <size_t NbImpulsions>
		std::vector<std::pair<SimuCore::Structures::Rocket, size_t>> population_to_rockets(const Pop<NbImpulsions>& population, const SimuCore::Systems::AdaptedSystem& system) {
			const size_t pop_size = population.size();
			std::vector<std::pair<SimuCore::Structures::Rocket, size_t>> rockets_and_indices;
			rockets_and_indices.reserve(pop_size);

			for (size_t i = 0; i < pop_size; i++) {
				auto& ind = population[i];
				
				auto copy_system = system;
				auto [rocket, state] = IndividualToRocket(ind.to_real_vectors(), copy_system);

				
				// on remet les fusées dans un état où la Terre est sur l'axe x pour que la distance à la planète finale soit plus facilement interprétable et pour que les clusters soient plus facilement interprétables
				rotate_rocket(&rocket, copy_system);
				rockets_and_indices.emplace_back(rocket, i);
			}

			return rockets_and_indices;
		}

		using Cluster = std::vector<int>;

		template <size_t NbImpulsions>
		std::vector<Cluster> HAC(const Pop<NbImpulsions>& population, const SimuCore::Systems::AdaptedSystem& system) {
			std::vector<std::pair<SimuCore::Structures::Rocket, size_t>> rockets_and_indices = population_to_rockets<NbImpulsions>(population, system);

			constexpr double dissimilarity_threshold = 12; // à ajuster en fonction des résultats

			// on initialise les patates avec une fusée chacune
			std::vector<Patate> patates;
			patates.reserve(rockets_and_indices.size());

			for (const auto& [rocket, indice] : rockets_and_indices) {
				Patate patate(5); // TODO : à modifier
				patate.AddRocket(&rocket, indice);

				patates.emplace_back(patate);
			}

			// tant que la dissimilitude entre les patates est inférieure au seuil, on fusionne les deux patates les plus proches
			double dissim = std::numeric_limits<double>::min();

			while (dissim < dissimilarity_threshold && patates.size() > 1) {

				size_t idx1 = 0, idx2 = 0;
				dissim = std::numeric_limits<double>::max();

				// on parcourt toutes les paires de patates pour trouver les deux patates les plus proches
				for (size_t i = 0; i < patates.size(); i++) {
					for (size_t j = i + 1; j < patates.size(); j++) {
						double d = dissimilarity(patates[i], patates[j]);
						if (d < dissim) {
							dissim = d;
							idx1 = i;
							idx2 = j;
						}
					}
				}

				if (dissim < dissimilarity_threshold) {
					patates[idx1].Merge(patates[idx2]);
					patates.erase(patates.begin() + idx2);
				}
			}

			// on retourne les clusters composés des pointeurs vers les individus de la population
			std::vector<Cluster> clusters;

			for (const auto& patate : patates) {
				Cluster cluster;
				cluster.reserve(patate.size());

				for (const auto& [rocket_ptr, indice] : patate.getRockets()) {
					cluster.push_back(indice);
				}
				clusters.emplace_back(cluster);
			}
			return clusters;
		}

	}; // namespace Statistics
}; // namespace SimuCore
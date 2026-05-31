#pragma once

#include "pch.h"
#include "SimuCore/utility.h"
#include "SimuCore/constants.h"
#include <SimuCore\structures\structures.h>
#include <SimuCore\theory\formula.h>

namespace SimuCore {
	struct Info {
		double potential;
		double kinetic;
		glm::dvec3 position;

		Info(glm::dvec3 pos, double energy_kin, double energy_pot) {
			position = pos;
			potential = energy_pot;
			kinetic = energy_kin;
		}
	};

	using MovementInfo = std::vector<std::pair < std::string, std::vector<Info>>>;

	namespace Integrator {

		using Entity = SimuCore::Structures::Entity;
		using Planet = SimuCore::Structures::Planet;
		using Rocket = SimuCore::Structures::Rocket;

		using Integrator = std::function<void(std::vector<Entity>&, double)>;


		enum class IntegratorMethod: uint8_t {
			EULER,
			VERLET,
			VELOCITY_VERLET,
			RUNGE_KUTTA_4
		};

		MovementInfo simulation(double duration, double dt,
			const std::vector<Planet>& planets,
			const Rocket& rocket);

		void IntegrateStep(std::vector<Entity*>& entities, double dt, double t_simu, IntegratorMethod integratorMethod = IntegratorMethod::VELOCITY_VERLET);

		void VelocityVerletIntegrator(std::vector<Entity*>& entities, double dt, double t_simu);

		void CalculateForces(std::vector<Entity*>& entities, double t_simu);


		using RocketState = SimuCore::Systems::AdaptedSystem::RocketState;

		/// <summary>
		/// Simule l'évolution d'une fusée sur une durée donnée et renvoie son état final.
		/// </summary>
		/// <param name="dt">Pas de temps de la simulation (secondes).</param>
		/// <param name="max_time">Durée maximale de la simulation (jours).</param>
		/// <param name="mass">Masse de la fusée (kg).</param>
		/// <param name="initial_pos">Position initiale de la fusée (UA).</param>
		/// <param name="initial_velocity">Vitesse initiale de la fusée (km/s).</param>
		/// <param name="final_planet_pos_indice">Indice de la position finale relative à la planète cible.</param>
		/// <param name="nb_impulsions">Nombre d'impulsions/manœuvres prévues.</param>
		/// <param name="impulsions">Liste des vecteurs d'impulsion appliqués (km/s), un par manœuvre.</param>
		/// <param name="dates">Liste des instants (jours) d'application des impulsions, correspondant à 'impulsions'.</param>
		/// <param name="startPlanet">Planète de départ (SimuCore::Systems::PlanetsName).</param>
		/// <param name="finalPlanet">Planète cible/finale (SimuCore::Systems::PlanetsName).</param>
		/// <param name="arretPasNeutre">Si vrai, arrête la simulation lorsqu'un pas neutre est atteint (comportement dépendant de l'implémentation).</param>
		/// <param name="positions">Pointeur optionnel vers un vecteur pour enregistrer les positions en coordonnées (UA) au cours de la simulation.</param>
		/// <returns>État final de la fusée (SimuCore::Systems::AdaptedSystem::RocketState) après la simulation.</returns>
		SimuCore::Systems::AdaptedSystem::RocketState Simulate(double dt, double max_time, double mass, glm::dvec3 initial_pos, glm::dvec3 initial_velocity, size_t final_planet_pos_indice, uint8_t nb_impulsions, std::vector<glm::dvec3> impulsions, std::vector<double> dates, SimuCore::Systems::PlanetsName startPlanet, SimuCore::Systems::PlanetsName finalPlanet, bool arretPasNeutre=false, std::vector<glm::dvec2>* positions=nullptr);
	}
};

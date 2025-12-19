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
	}
};

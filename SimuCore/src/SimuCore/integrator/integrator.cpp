#include <pch.h>
#include <SimuCore/integrator/integrator.h>

namespace SimuCore::Integrator {

	MovementInfo simulation(double duration, double dt,
		const std::vector<Planet>& planets,
		const Rocket& rocket) {

		const size_t NB_ITERATION = static_cast<const size_t>(duration / dt);
		const size_t NB_PLANETS = planets.size();

		double time = 0.0;

		std::vector<std::unique_ptr<Entity>> entities;
		entities.reserve(1 + NB_PLANETS);

		entities.push_back(std::make_unique<Rocket>(rocket));
		for (int i = 0; i < NB_PLANETS; i++) {
			entities.push_back(std::make_unique<Planet>(planets[i]));
		}

		////////////////////////////////////////
		////////////////  INIT  ////////////////
		////////////////////////////////////////

		MovementInfo res;
		res.reserve(1 + NB_PLANETS);

		std::vector<Info> rocket_info;
		rocket_info.reserve(NB_ITERATION);

		res.push_back({ "rocket" ,  rocket_info });

		for (size_t i = 0; i < NB_PLANETS; i++) {
			std::vector<Info> planet_info;
			planet_info.reserve(NB_ITERATION);

			res.push_back({ planets[i].Name() ,  planet_info });
		}


		////////////////////////////////////////
		////////////////  SIMU  ////////////////
		////////////////////////////////////////

		for (size_t k = 1; k < NB_ITERATION; k++) {
			// on simule l'étape actuelle pour toute les entités

			// 1) on calcul les forces
			for (size_t i = 0; i < NB_PLANETS + 1; i++) {

				// on vérifie que l'entité i est encore "vivante"
				if (!entities[i]->IsAlive(time)) continue;

				for (size_t j = i + 1; j < NB_PLANETS + 1; j++) {
					// on vérifie que l'entité j est encore "vivante"
					if (!entities[j]->IsAlive(time)) continue;

					// on calcul la force d'attraction de j sur i
					glm::dvec3 attractionForce = forceAttractionGrav(*entities[j], *entities[i]);

					// on applique les 2 forces
					entities[i]->forces += attractionForce;
					entities[j]->forces -= attractionForce;
				}
				// 2) on peut bouger la planète i
				entities[i]->UpdateFirstPart(dt);
			}

			// 3) on refait la même chose pour la seconde partie de l'algo LeapFrog
			for (size_t i = 0; i < NB_PLANETS + 1; i++) {
				// on vérifie que l'entité i est encore "vivante"
				if (!entities[i]->IsAlive(time)) continue;

				for (size_t j = i + 1; j < NB_PLANETS + 1; j++) {
					// on vérifie que l'entité j est encore "vivante"
					if (!entities[j]->IsAlive(time)) continue;

					// on calcul la force d'attraction de j sur i
					glm::dvec3 attractionForce = forceAttractionGrav(*entities[j], *entities[i]);

					// on applique les 2 forces
					entities[i]->forces += attractionForce;
					entities[j]->forces -= attractionForce;
				}
				// 2) on peut bouger la planète i
				entities[i]->UpdateSecondPart(dt);
			}
			
			time += dt;

			// 4) on sauvegarde les données
			for (size_t i = 0; i < NB_PLANETS + 1; i++) {
				// on calcul l'énergie cinétique et potentielle
				double potential = 0.0;
				for (size_t j = 0; j < entities.size(); ++j) {
					if (j == i) continue;
					const glm::dvec3 r = entities[j]->position - entities[i]->position;
					double dist = glm::length(r) + constants::softening;
					potential -= constants::G * entities[i]->mass * entities[j]->mass / dist;
				}

				const double speed = glm::length(entities[i]->velocity);
				res[i].second.push_back(Info(entities[i]->position, 0.5 * entities[i]->mass * speed * speed, potential));
			}
		}

		return res;
	}

	void IntegrateStep(std::vector<Entity*>& entities, double dt, double t_simu, IntegratorMethod integratorMethod) {
		switch (integratorMethod)
		{
		case SimuCore::Integrator::IntegratorMethod::EULER:
			throw std::invalid_argument("Euler integrator not implemented yet.");
			break;
		case SimuCore::Integrator::IntegratorMethod::VERLET:
			throw std::invalid_argument("Verlet integrator not implemented yet.");
			break;
		case SimuCore::Integrator::IntegratorMethod::VELOCITY_VERLET:
			VelocityVerletIntegrator(entities, dt, t_simu);
			break;
		case SimuCore::Integrator::IntegratorMethod::RUNGE_KUTTA_4:
			throw std::invalid_argument("Runge-Kutta 4 integrator not implemented yet.");
			break;
		default:
			break;
		}
	}

	void VelocityVerletIntegrator(std::vector<Entity*>& entities, double dt, double t_simu) {
		// Étape 0 : Calculer les forces F_n (à l'instant t_n)
		CalculateForces(entities, t_simu);

		// Étape 1 : Mise à jour de la demi-vitesse (v_{n+1/2}) et de la position (x_{n+1})
		for (size_t i = 0; i < entities.size(); i++) {
			if (!entities[i]->IsAlive(t_simu)) continue;
			entities[i]->UpdateFirstPart(dt);
		}

		// Application des Impulsions pour la Fusée
		// t_simu est le temps à la fin de ce pas de temps (t_{n+1})
		Rocket* rocket = dynamic_cast<Rocket*>(entities[0]);
		if (rocket != nullptr) {
			rocket->ApplyImpulsions(t_simu, dt);
		}


		// Étape 2 : Calculer les nouvelles forces F_{n+1} aux nouvelles positions x_{n+1}
		CalculateForces(entities, t_simu);

		// Étape 3 : Mise à jour de la pleine vitesse v_{n+1}
		for (size_t i = 0; i < entities.size(); i++) {
			if (!entities[i]->IsAlive(t_simu)) continue;
			entities[i]->UpdateSecondPart(dt);
		}
	}

	void CalculateForces(std::vector<Entity*>& entities, double t_simu) {
		// 1. Réinitialiser toutes les forces avant de calculer les nouvelles.
		for (Entity* entity : entities) {
			entity->forces = glm::dvec3(0);
		}

		// 2. Calculer les interactions gravitationnelles
		for (size_t i = 0; i < entities.size(); i++) {
			if (!entities[i]->IsAlive(t_simu)) continue;

			for (size_t j = i + 1; j < entities.size(); j++) {
				if (!entities[j]->IsAlive(t_simu)) continue;

				glm::dvec3 attractionForce = forceAttractionGrav(*entities[j], *entities[i]);

				entities[i]->forces += attractionForce;
				entities[j]->forces -= attractionForce;
			}
		}
	}
};
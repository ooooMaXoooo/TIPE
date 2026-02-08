#pragma once
#include <pch.h>
#include <SimuCore/Units/all.h>


namespace SimuCore {
	namespace Structures {

		struct Entity {
			glm::dvec3 position;	// position de l'entité (AU)
			glm::dvec3 velocity;	// vitesse de l'entité (km/s)
			double mass;			// masse de l'entité (kg)
			glm::dvec3 forces;		// forces appliquées à l'entité (kN) (kg*km/s²)


			/**
			 * @brief Constructeur de la classe Entity.
			 * @param m Masse de la planète. (kg)
			 * @param p0 Position initiale (glm::dvec3). (AU)
			 * @param v0 Vitesse initiale (glm::dvec3). (km/s)
			 */
			Entity(double m, glm::dvec3 p0 = glm::dvec3(0, 0, 0), glm::dvec3 v0 = glm::dvec3(0, 0, 0)) {
				mass = m;
				position = p0;
				velocity = v0;
				forces = glm::dvec3(0);
			}

			Entity& operator=(const Entity& other) {
				if (this != &other) {
					this->position = other.position;
					this->velocity = other.velocity;
					this->mass = other.mass;
					this->forces = other.forces;
				}
				return *this;
			}

			/// <summary>
			/// Première partie de l'intégration de la position et de la vitesse de l'entité, selon le schéma de Verlet à 2 étapes.
			/// </summary>
			/// <param name="dt">pas de temps en s</param>
			virtual void UpdateFirstPart(double dt) {
				glm::dvec3 acceleration = forces / mass;								// a = F/m (en km/s²)
				velocity += 0.5 * dt * acceleration;									// a * dt/2 (en km/s)
				position += static_cast<double>((1._km_to_AU) * dt) * velocity;			// dt * v en km

				forces = glm::dvec3(0);
			}

			/// <summary>
			/// Seconde partie de l'intégration de la position et de la vitesse de l'entité, selon le schéma de Verlet à 2 étapes.
			/// </summary>
			/// <param name="dt">pas de temps en s</param>
			virtual void UpdateSecondPart(double dt) {
				glm::dvec3 acceleration = forces / mass;	// F en kN (kg*km/s²)   |  m en kg       = a en km/s²
				velocity += 0.5 * dt * acceleration;		// a en km/s²		    | dt en s        = a*dt en km/s

				forces = glm::dvec3(0);
			}

			virtual bool IsAlive(double current_time) const { return true; };
		};
	}
}
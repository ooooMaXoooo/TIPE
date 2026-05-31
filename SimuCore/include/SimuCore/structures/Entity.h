#pragma once
#include <pch.h>
#include <SimuCore/Units/all.h>


namespace SimuCore {
	namespace Structures {

		struct Entity {
			// position de l'entité (AU)
			glm::dvec3 position;

			// vitesse de l'entité (km/s)
			glm::dvec3 velocity;

			// masse de l'entité (kg)
			double mass;

			// forces appliquées à l'entité (kN) (kg*km/s²)
			glm::dvec3 forces;


			/**
			 * @brief Constructeur de la classe Entity.
			 * @param m Masse de la planète. (kg)
			 * @param p0 Position initiale (glm::dvec3). (AU)
			 * @param v0 Vitesse initiale (glm::dvec3). (km/s)
			 */
			Entity(double m, glm::dvec3 p0 = glm::dvec3(0, 0, 0), glm::dvec3 v0 = glm::dvec3(0, 0, 0));

			Entity& operator=(const Entity& other);

			/// <summary>
			/// Première partie de l'intégration de la position et de la vitesse de l'entité, selon le schéma de Verlet à 2 étapes.
			/// </summary>
			/// <param name="dt">pas de temps en s</param>
			virtual void UpdateFirstPart(double dt);

			/// <summary>
			/// Seconde partie de l'intégration de la position et de la vitesse de l'entité, selon le schéma de Verlet à 2 étapes.
			/// </summary>
			/// <param name="dt">pas de temps en s</param>
			virtual void UpdateSecondPart(double dt);

			virtual bool IsAlive(double current_time) const { return true; };

			virtual void Rotate(double theta);
		};
	}
}
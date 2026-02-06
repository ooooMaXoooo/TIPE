#pragma once
#include <pch.h>


namespace SimuCore {
	namespace Structures {

		struct Entity {
			glm::dvec3 position;
			glm::dvec3 velocity;
			double mass;
			glm::dvec3 forces;


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

			virtual void UpdateFirstPart(double dt) {
				glm::dvec3 acceleration = forces / mass;
				velocity += 0.5 * dt * acceleration;
				position += dt * velocity;

				forces = glm::dvec3(0);
			}

			virtual void UpdateSecondPart(double dt) {
				glm::dvec3 acceleration = forces / mass;
				velocity += 0.5 * dt * acceleration;

				forces = glm::dvec3(0);
			}

			virtual bool IsAlive(double current_time) const { return true; };
		};
	}
}
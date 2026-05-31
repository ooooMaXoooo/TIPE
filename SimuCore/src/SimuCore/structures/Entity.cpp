#include "pch.h"
#include "structures/Entity.h"

namespace SimuCore {
	namespace Structures {

		Entity::Entity(double m, glm::dvec3 p0, glm::dvec3 v0) {
			mass = m;
			position = p0;
			velocity = v0;
			forces = glm::dvec3(0);
		}

		Entity& Entity::operator=(const Entity& other) {
			if (this != &other) {
				this->position = other.position;
				this->velocity = other.velocity;
				this->mass = other.mass;
				this->forces = other.forces;
			}
			return *this;
		}

		void Entity::UpdateFirstPart(double dt) {
			glm::dvec3 acceleration = forces / mass;								// a = F/m (en km/s²)
			velocity += 0.5 * dt * acceleration;									// a * dt/2 (en km/s)
			position += static_cast<double>((1._km_to_AU) * dt) * velocity;			// dt * v en km

			forces = glm::dvec3(0);
		}

		void Entity::UpdateSecondPart(double dt) {
			glm::dvec3 acceleration = forces / mass;	// F en kN (kg*km/s²)   |  m en kg       = a en km/s²
			velocity += 0.5 * dt * acceleration;		// a en km/s²		    | dt en s        = a*dt en km/s

			forces = glm::dvec3(0);
		}

		void Entity::Rotate(double theta) {
			// Rotation autour de l'axe z
			double cos_angle = std::cos(theta);
			double sin_angle = std::sin(theta);

			glm::dvec3 new_position;
			{
				new_position.x = position.x * cos_angle - position.y * sin_angle;
				new_position.y = position.x * sin_angle + position.y * cos_angle;
			}
			position = new_position;

			glm::dvec3 new_velocity;
			{
				new_velocity.x = velocity.x * cos_angle - velocity.y * sin_angle;
				new_velocity.y = velocity.x * sin_angle + velocity.y * cos_angle;
			}
			velocity = new_velocity;
		}

	}; // namespace Structures
} // namespace SimuCore
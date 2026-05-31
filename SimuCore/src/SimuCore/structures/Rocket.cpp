#include <pch.h>
#include <SimuCore\structures\Rocket.h>


namespace SimuCore {
	namespace Structures {

		Impulsion operator-(const Impulsion& I1, const Impulsion& I2) {
			return Impulsion(I1.velocity - I2.velocity);
		}

		Rocket operator-(const Rocket& r1, const Rocket& r2) {

			std::vector<std::pair<Impulsion, double>> impulsions;
			impulsions.reserve(r1.m_Impulsions.size());

			for (int i = 0; i < r1.m_Impulsions.size(); i++) {
				auto& [imp1, date1] = r1.m_Impulsions[i];
				auto& [imp2, date2] = r2.m_Impulsions[i];

				impulsions.emplace_back(imp1 - imp2, date1 - date2);
			}


			return Rocket(r1.lifetime, impulsions, r1.mass, r1.m_Vitesse_ejection_gaz, r1.position - r2.position, r1.velocity - r2.velocity);
		}

		Rocket::Rocket(double _lifetime,
			std::vector<std::pair<Impulsion, double>> impulsions,
			double m,
			double vitesse_ejection_gaz,
			glm::dvec3 p0,
			glm::dvec3 v0)
			:
			Entity(m, p0, v0),
			lifetime(_lifetime),
			m_Impulsions(impulsions),
			acceleration(0),
			m_Vitesse_ejection_gaz(vitesse_ejection_gaz)
		{}

		Rocket& Rocket::operator=(const Rocket& other) {
			if (this != &other) {
				// Appeler l'opérateur d'affectation de la classe de base
				Entity::operator=(other);
				// Copier les membres spécifiques de Rocket
				this->lifetime = other.lifetime;
				this->acceleration = other.acceleration;
				this->m_Impulsions = other.m_Impulsions;
				this->m_Vitesse_ejection_gaz = other.m_Vitesse_ejection_gaz;
			}
			return *this;
		}

		void Rocket::UpdateFirstPart(double dt) {
			glm::dvec3 acceleration = forces / mass;	// F en kN (kg*km/s²)   |  m en kg       = a en km/s²
			velocity += 0.5f * dt * acceleration;		// a en km/s²		    | dt en s        = a*dt en km/s
			position += static_cast<double>(1.0_km_to_AU) * dt * velocity; // dt en s, v en km/s, position en AU

			//std::cout << __FUNCSIG__ << '\n';

			forces = glm::dvec3(0);
		}

		void Rocket::ApplyImpulsions(double t_current, double dt) {
			t_current = days_to_seconds(t_current); // en secondes

			for (const auto& pair : m_Impulsions) {
				const auto& t_impulsion = days_to_seconds(pair.second);	// en secondes
				const auto& impulsion = pair.first;		// en km/s

				// Si l'impulsion t_impulsion se trouve dans l'intervalle [t_current, t_current+dt[
				if (t_current <= t_impulsion && t_impulsion < t_current + dt) {
					// Appliquer l'impulsion au v_{n+1/2} (la vitesse qui vient d'être calculée)
					impulsion.Apply(this);
				}
			}
		}

		void Rocket::UpdateSecondPart(double dt) {
			glm::dvec3 _acceleration = forces / mass;	// F en kN (kg*km/s²)   |  m en kg       = a en km/s²
			velocity += 0.5f * dt * _acceleration;		// a en km/s²		    | dt en s        = a*dt en km/s

			acceleration = glm::length(_acceleration);	// en km/s²
			//std::cout << __FUNCSIG__ << '\n';

			forces = glm::dvec3(0);
		}

		double Rocket::getDeltaV(double time) const {
			double deltaV = 0;
			for (auto& [impuls, t] : m_Impulsions) {
				if (t <= time)
					deltaV += impuls.Length(); // en km/s
			}

			return deltaV;
		}

		void Rocket::Rotate(double cos_theta, double sin_theta) {
			glm::dvec3 new_position;
			{
				new_position.x = position.x * cos_theta - position.y * sin_theta;
				new_position.y = position.x * sin_theta + position.y * cos_theta;
			}
			position = new_position;

			glm::dvec3 new_velocity;
			{
				new_velocity.x = velocity.x * cos_theta - velocity.y * sin_theta;
				new_velocity.y = velocity.x * sin_theta + velocity.y * cos_theta;
			}
			velocity = new_velocity;

			for (auto& [impulsion, _] : m_Impulsions) {
				const auto& dv = impulsion.GetDeltaV_vec();
				impulsion.setImpulsion({
					dv.x * cos_theta - dv.y * sin_theta,
					dv.x * sin_theta + dv.y * cos_theta,
					dv.z
				});
			}
		}

		void Rocket::Rotate(double theta) {
			// Rotation autour de l'axe z
			double cos_angle = std::cos(theta);
			double sin_angle = std::sin(theta);
			Rotate(cos_angle, sin_angle);
		}

		void Rocket::Rotate(const glm::dvec3& pos) {
			double norm = std::sqrt(pos.x * pos.x + pos.y * pos.y);
			// cos(-θ) = cos(θ) = x/r,  sin(-θ) = -y/r
			Rotate(pos.x / norm, -pos.y / norm);
		}

		glm::dvec3 Rocket::GetInitialImpulsion() const {
			return m_Impulsions.empty() ? glm::dvec3(0) : m_Impulsions[0].first.GetDeltaV_vec();
		}

	}; // namespace Structures
}; // namespace SimuCore
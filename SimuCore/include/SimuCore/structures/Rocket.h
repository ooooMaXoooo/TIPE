#pragma once
#include <pch.h>
#include <SimuCore\structures\Entity.h>
#include <SimuCore\utility.h>

namespace SimuCore {
	namespace Structures {

		class Impulsion {
			glm::dvec3 velocity;

		public:
			Impulsion(glm::dvec3 v) : velocity(v) {}

			void Apply(Entity* entity) const {
				entity->velocity += velocity;
			}

			double Length() const { return glm::length(velocity); }
		};

		struct Rocket : public Entity {
		private:
			std::vector<std::pair<Impulsion, double>> m_Impulsions; // en m.s-1
			double m_Vitesse_ejection_gaz; // en m.s-1

		public:
			double lifetime; // en s
			double acceleration; // en m.s-2

			/// <summary>
			/// 
			/// </summary>
			/// <param name="_lifetime"> temps de simulation de la fusée en jour </param>
			/// <param name="impulsions"> listes des impulsions (modelisé par ajout de vitesse) en km.s-2 </param>
			/// <param name="m"> masse de la fusée en kg </param>
			/// <param name="vitesse_ejection_gaz"> en km.s-1 </param>
			/// <param name="p0"> position initiale en km </param>
			/// <param name="v0"> vitesse initiale en km.s-1 </param>
			Rocket(double _lifetime, std::vector<std::pair<Impulsion, double>> impulsions,
				double m, double vitesse_ejection_gaz, glm::dvec3 p0 = glm::dvec3(0, 0, 0), glm::dvec3 v0 = glm::dvec3(0, 0, 0))
				: Entity(m, p0, v0), lifetime(daysInSeconds(_lifetime)), m_Impulsions(impulsions), acceleration(0), m_Vitesse_ejection_gaz (vitesse_ejection_gaz) {
			}

			Rocket& operator=(const Rocket& other) {
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

			virtual void UpdateFirstPart(double dt) override {
				glm::dvec3 acceleration = forces / mass;
				velocity += 0.5f * dt * acceleration;
				position += dt * velocity;

				//std::cout << __FUNCSIG__ << '\n';

				forces = glm::dvec3(0);
			}

			void ApplyImpulsions(double t_current, double dt) {
				for (const auto& pair : m_Impulsions) {
					const auto& t_impulsion = pair.second;
					const auto& impulsion = pair.first;

					// Si l'impulsion t_impulsion se trouve dans l'intervalle [t_current, t_current+dt[
					if (t_current <= t_impulsion && t_impulsion < t_current+dt) {
						// Appliquer l'impulsion au v_{n+1/2} (la vitesse qui vient d'être calculée)
						impulsion.Apply(this);
					}
				}
			}

			virtual void UpdateSecondPart(double dt) override {
				glm::dvec3 _acceleration = forces / mass;
				velocity += 0.5f * dt * _acceleration;

				acceleration = glm::length(_acceleration);
				//std::cout << __FUNCSIG__ << '\n';

				forces = glm::dvec3(0);
			}

			virtual bool IsAlive(double current_time) const override { return current_time < lifetime; }


			double getDeltaV() const {
				double deltaV = 0;
				for (auto& [impuls, t] : m_Impulsions) {
					deltaV += impuls.Length();
				}

				return deltaV;
			}

			double getDeltaM() const {
				double coef_alpha = std::exp(getDeltaV() / m_Vitesse_ejection_gaz);

				return mass * ((coef_alpha - 1) / coef_alpha); // on ne change jamais la masse car il n'y a aucun effet sur le PFD
			}

			void setImpulsions(const std::vector<std::pair<Impulsion, double>>& impulsions) {
				m_Impulsions = impulsions;
			}
		};
	};
};
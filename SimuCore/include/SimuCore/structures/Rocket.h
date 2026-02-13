#pragma once
#include <pch.h>
#include <SimuCore\structures\Entity.h>
#include <SimuCore\utility.h>

namespace SimuCore {
	namespace Structures {

		class Impulsion {
			glm::dvec3 velocity;	// en km/s

		public:

			/// <summary>
			/// Constructeur qui initialise le membre velocity avec le vecteur 3D fourni.
			/// </summary>
			/// <param name="v"> velocity (km/s)</param>
			Impulsion(glm::dvec3 v) : velocity(v) {}

			void Apply(Entity* entity) const {
				entity->velocity += velocity;
			}

			/// <summary>
			/// Renvoie la norme du vecteur membre velocity.
			/// </summary>
			/// <returns>La norme du vecteur velocity, retournée comme double. (en km/s) </returns>
			double Length() const { return glm::length(velocity); }

			/// <summary>
			/// 
			/// </summary>
			/// <returns> km/s </returns>
			const glm::dvec3& GetDeltaV_vec() const noexcept { return velocity; }
		};

		struct Rocket : public Entity {
		private:
			std::vector<std::pair<Impulsion, double>> m_Impulsions; // une impulsion et un instant en km/s et en jours
			double m_Vitesse_ejection_gaz; // en km/s

		public:
			double lifetime; // en jours
			double acceleration; // en km/s²

			/// <summary>
			/// Constructeur de la classe Rocket, qui initialise les membres de la classe Entity (position, vitesse, masse) ainsi que les membres spécifiques à la classe Rocket (lifetime, impulsions, vitesse d'éjection des gaz).
			/// </summary>
			/// <param name="_lifetime"> temps de simulation de la fusée. (en jours)</param>
			/// <param name="impulsions"> listes des impulsions (modelisé par ajout de vitesse) en km/s et en jours </param>
			/// <param name="m"> masse de la fusée en kg </param>
			/// <param name="vitesse_ejection_gaz"> en km/s </param>
			/// <param name="p0"> position initiale en UA </param>
			/// <param name="v0"> vitesse initiale en km/s </param>
			Rocket(
				double _lifetime,
				std::vector<std::pair<Impulsion, double>> impulsions,
				double m, 
				double vitesse_ejection_gaz,
				glm::dvec3 p0 = glm::dvec3(0, 0, 0),
				glm::dvec3 v0 = glm::dvec3(0, 0, 0)
			)
				:
				Entity(m, p0, v0),
				lifetime(_lifetime),
				m_Impulsions(impulsions),
				acceleration(0),
				m_Vitesse_ejection_gaz (vitesse_ejection_gaz)
			{
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

			/// <summary>
			/// Première partie de l'intégration de la position et de la vitesse de l'entité, selon le schéma de Verlet à 2 étapes.
			/// </summary>
			/// <param name="dt">pas de temps en s</param>
			virtual void UpdateFirstPart(double dt) override {
				glm::dvec3 acceleration = forces / mass;	// F en kN (kg*km/s²)   |  m en kg       = a en km/s²
				velocity += 0.5f * dt * acceleration;		// a en km/s²		    | dt en s        = a*dt en km/s
				position += static_cast<double>(1.0_km_to_AU) * dt * velocity; // dt en s, v en km/s, position en AU

				//std::cout << __FUNCSIG__ << '\n';

				forces = glm::dvec3(0);
			}

			/// <summary>
			/// Parcourt la liste d'impulsions et applique à l'objet courant celles dont le temps d'application t_impulsion est dans l'intervalle [t_current, t_current + dt[.
			/// </summary>
			/// <param name="t_current">Temps courant (en jours) </param>
			/// <param name="dt">Durée du pas de temps ; l'intervalle considéré est [t_current, t_current + dt[ (dt en secondes) .</param>
			void ApplyImpulsions(double t_current, double dt) {
				t_current = days_to_seconds(t_current); // en secondes

				for (const auto& pair : m_Impulsions) {
					const auto& t_impulsion = days_to_seconds(pair.second);	// en secondes
					const auto& impulsion = pair.first;		// en km/s

					// Si l'impulsion t_impulsion se trouve dans l'intervalle [t_current, t_current+dt[
					if (t_current <= t_impulsion && t_impulsion < t_current+dt) {
						// Appliquer l'impulsion au v_{n+1/2} (la vitesse qui vient d'être calculée)
						impulsion.Apply(this);
					}
				}
			}

			/// <summary>
			/// Seconde partie de l'intégration de la position et de la vitesse de l'entité, selon le schéma de Verlet à 2 étapes.
			/// 
			/// Attention les forces doivent être RECALCULEE après UpdateFirstPart
			/// </summary>
			/// <param name="dt">pas de temps en s</param>
			virtual void UpdateSecondPart(double dt) override {
				glm::dvec3 _acceleration = forces / mass;	// F en kN (kg*km/s²)   |  m en kg       = a en km/s²
				velocity += 0.5f * dt * _acceleration;		// a en km/s²		    | dt en s        = a*dt en km/s

				acceleration = glm::length(_acceleration);	// en km/s²
				//std::cout << __FUNCSIG__ << '\n';

				forces = glm::dvec3(0);
			}

			/// <summary>
			/// Indique si l'objet est encore actif à l'instant fourni en comparant current_time avec sa durée de vie (lifetime).
			/// </summary>
			/// <param name="current_time">en jours</param>
			/// <returns>true si current_time est inférieur à lifetime (l'objet est encore vivant), false sinon.</returns>
			virtual bool IsAlive(double current_time) const override { return current_time < lifetime; }


			/// <summary>
			/// Calcule la variation de vitesse totale (ΔV) en sommant les longueurs des vecteurs d'impulsion stockés dans m_Impulsions.
			/// </summary>
			/// <returns>La variation de vitesse totale (somme des longueurs des impulsions). (en km/s)</returns>
			double getDeltaV() const {
				double deltaV = 0;
				for (auto& [impuls, t] : m_Impulsions) {
					deltaV += impuls.Length(); // en km/s
				}

				return deltaV;
			}


			/// <summary>
			/// Calcule la variation de masse (Δm) nécessaire pour produire le Δv obtenu par getDeltaV(), en utilisant coef_alpha = exp(getDeltaV() / m_Vitesse_ejection_gaz). La méthode est const et n'altère pas l'état de l'objet.
			/// </summary>
			/// <returns>La variation de masse (Δm) calculée comme mass * ((coef_alpha - 1) / coef_alpha), où coef_alpha = exp(getDeltaV() / m_Vitesse_ejection_gaz). (en kg)</returns>
			double getDeltaM() const {
				double coef_alpha = std::exp(getDeltaV() / m_Vitesse_ejection_gaz);

				return mass * ((coef_alpha - 1) / coef_alpha); // on ne change jamais la masse car il n'y a aucun effet sur le PFD
			}

			/// <summary>
			/// Remplace la collection interne d'impulsions par le vecteur fourni.
			/// </summary>
			/// <param name="impulsions">Vecteur de paires (Impulsion, double) représentant les impulsions et la valeur numérique associée ; le contenu remplace l'état interne m_Impulsions.(en km/s et en jours)</param>
			void setImpulsions(const std::vector<std::pair<Impulsion, double>>& impulsions) {
				m_Impulsions = impulsions;
			}

			/// <summary>
			/// Renvoi une référence constante au tableau des impulsions
			/// </summary>
			/// <returns> km/s et jours </returns>
			const std::vector<std::pair<Impulsion, double>>& getImpulsions() const noexcept { return m_Impulsions; }
		};
	};
};
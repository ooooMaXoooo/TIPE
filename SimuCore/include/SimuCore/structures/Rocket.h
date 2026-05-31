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

			inline void Apply(Entity* entity) const { entity->velocity += velocity;}

			/// <summary>
			/// Renvoie la norme du vecteur membre velocity.
			/// </summary>
			/// <returns>La norme du vecteur velocity, retournée comme double. (en km/s) </returns>
			inline double Length() const { return glm::length(velocity); }

			/// <summary>
			/// 
			/// </summary>
			/// <returns> km/s </returns>
			inline const glm::dvec3& GetDeltaV_vec() const noexcept { return velocity; }

			inline void setImpulsion(const glm::dvec3& v) { velocity = v; }

			friend Impulsion operator-(const Impulsion& I1, const Impulsion& I2);
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
			);

			Rocket& operator=(const Rocket& other);

			friend Rocket operator-(const Rocket& r1, const Rocket& r2);


			/// <summary>
			/// Première partie de l'intégration de la position et de la vitesse de l'entité, selon le schéma de Verlet à 2 étapes.
			/// </summary>
			/// <param name="dt">pas de temps en s</param>
			virtual void UpdateFirstPart(double dt) override;

			/// <summary>
			/// Parcourt la liste d'impulsions et applique à l'objet courant celles dont le temps d'application t_impulsion est dans l'intervalle [t_current, t_current + dt[.
			/// </summary>
			/// <param name="t_current">Temps courant (en jours) </param>
			/// <param name="dt">Durée du pas de temps ; l'intervalle considéré est [t_current, t_current + dt[ (dt en secondes) .</param>
			void ApplyImpulsions(double t_current, double dt);

			/// <summary>
			/// Seconde partie de l'intégration de la position et de la vitesse de l'entité, selon le schéma de Verlet à 2 étapes.
			/// 
			/// Attention les forces doivent être RECALCULEE après UpdateFirstPart
			/// </summary>
			/// <param name="dt">pas de temps en s</param>
			virtual void UpdateSecondPart(double dt) override;

			/// <summary>
			/// Indique si l'objet est encore actif à l'instant fourni en comparant current_time avec sa durée de vie (lifetime).
			/// </summary>
			/// <param name="current_time">en jours</param>
			/// <returns>true si current_time est inférieur à lifetime (l'objet est encore vivant), false sinon.</returns>
			inline virtual bool IsAlive(double current_time) const override { return current_time < lifetime; }


			/// <summary>
			/// Calcule la variation de vitesse totale (ΔV) en sommant les longueurs des vecteurs d'impulsion stockés dans m_Impulsions.
			/// </summary>
			/// <param name="time"> en jours</param>
			/// <returns>La variation de vitesse totale (somme des longueurs des impulsions). (en km/s)</returns>
			double getDeltaV(double time) const;


			/// <summary>
			/// Remplace la collection interne d'impulsions par le vecteur fourni.
			/// </summary>
			/// <param name="impulsions">Vecteur de paires (Impulsion, double) représentant les impulsions et la valeur numérique associée ; le contenu remplace l'état interne m_Impulsions.(en km/s et en jours)</param>
			inline void setImpulsions(const std::vector<std::pair<Impulsion, double>>& impulsions) { m_Impulsions = impulsions;}

			/// <summary>
			/// Renvoi une référence constante au tableau des impulsions
			/// </summary>
			/// <returns> km/s et jours </returns>
			inline const std::vector<std::pair<Impulsion, double>>& getImpulsions() const noexcept { return m_Impulsions; }

			void Rotate(double cos_theta, double sin_theta);

			/// <summary>
			/// Effectue une rotation selon l'angle fourni.
			/// </summary>
			/// <param name="theta">Angle de rotation en radians.</param>
			void Rotate(double theta) override;

			/// <summary>
			/// Effectue une rotation sur la fusée d'un angle permettant d'aligner le vecteur pos sur l'axe x.
			/// </summary>
			/// <param name="pos"></param>
			void Rotate(const glm::dvec3& pos);

			glm::dvec3 GetInitialImpulsion() const;
		};		
	}; // namespace Structures
}; // namespace SimuCore
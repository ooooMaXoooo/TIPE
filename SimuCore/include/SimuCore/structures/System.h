#pragma once

#include "Planet.h"
#include "Entity.h"
#include "Rocket.h"

#include <array>
#include <string>




namespace SimuCore {
	namespace Systems {
		enum class PlanetsName : uint8_t {
			Mercure,
			Venus,
			Terre,
			Mars,
			Jupiter,
			Saturne,
			Uranus,
			Neptune
		};

		class AdaptedSystem {
			using Entity = SimuCore::Structures::Entity;
			using Planet = SimuCore::Structures::Planet;
			using Rocket = SimuCore::Structures::Rocket;

			inline static constexpr uint8_t m_NbPlanets = 9;
			static const std::array<Planet, m_NbPlanets> m_planets; // fusée + soleil + planète depart + planète arrivée

			Planet m_startPlanet;
			Planet m_finalPlanet;

			Rocket m_rocket;
			Planet m_sun;

			double m_time = 0;

			double m_startAngle;
			double m_finalAngle;

			double m_MaxTime;
			double m_deltaTime;

			using vecteur = glm::dvec3;

			const enum class RocketState : uint8_t {
				DEAD,
				NEUTRAL,
				VALID
			};

			const enum class ObjectName : uint8_t {
				SUN,
				START,
				FINAL
			};


		public:
			using Real = double;
			inline static constexpr double m_LowScore = -1e6;

			AdaptedSystem();

			AdaptedSystem(PlanetsName start_planet, PlanetsName final_planet, double start_angle, double final_angle, Rocket rocket, double max_duration, double dt_seconds);

			/**
			* @brief Simule le système à partir de son état courant
			*
			* @param state Une fonction booléenne, prenant en paramètre une fusée, et indiquant si son état est valide
			* @return l'état de la fusée à la fin de la simulation
			*/
			RocketState Run(std::function<RocketState(const Rocket&) > state); // TODO : Revoir la fonction de simulation pour modification
			void Reset();

			Real Score(const std::vector<std::vector<Real>>& individu);

			Real RingSize_meter() const;

		private:
			void InitPlanet(bool is_start_planet, PlanetsName name); // TODO AMELIORER L'IMPLEM
			void SetAnglePlanet(Planet& planet, double theta);
			void SetPlanetSpeed(Planet& planet, double theta);

			Real HandleScoreValidState() const;
			Real HandleScoreNeutralState() const;

			bool rocket_collide_with(ObjectName name) const;

		};
	};
};




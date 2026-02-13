#pragma once


#include <pch.h>

#include "Planet.h"
#include "Entity.h"
#include "Rocket.h"

#include <SimuCore/utility.h>




namespace SimuCore {
	const enum class GenerationState : uint8_t {
		VALID,
		INVALID_GENES_OUT_OF_BOUNDS_POSITION,
		INVALID_GENES_OUT_OF_BOUNDS_VELOCITY,
		INVALID_GENES_OUT_OF_BOUNDS_IMPULSION
	};

	namespace Systems {
		enum PlanetsName : uint8_t { /// l'indexation suppose que m_planets est initialisé dans le bon ordre
			Mercure = 1,
			Venus,
			Terre,
			Mars,
			Jupiter,
			Saturne,
			Uranus,
			Neptune
		};

		struct PlanetInfo {
			double distance_to_sun = 0; // en UA
			double angular_velocity = 0; // en rad/s
			double muPlanet = 0; // en m^3/s²
			size_t nb_iterations_orbit = 0; // nombre d'itérations pour faire une orbite complète (en fonction du pas de temps)
		};

		SimuCore::Structures::Planet getPlanetFromName(PlanetsName name);

		class AdaptedSystem {
			using Entity = SimuCore::Structures::Entity;
			using Planet = SimuCore::Structures::Planet;
			using Rocket = SimuCore::Structures::Rocket;

			Rocket m_rocket;

			inline static PlanetsName m_start_planet = PlanetsName::Terre;
			inline static PlanetsName m_final_planet = PlanetsName::Mars;

			double m_time = 0;							// en jours
			inline static double s_MaxTime = 50;		// en jours
			inline static double s_deltaTime = 1000;	// en secondes

			using vecteur = glm::dvec3;

			const enum class ObjectName : uint8_t {
				SUN,
				START,
				FINAL
			};

			static std::vector<vecteur> s_startPlanet_positions; // en UA (AU en anglais)
			static std::vector<vecteur> s_finalPlanet_positions; // en UA (AU en anglais)

			size_t m_start_planet_start_indice = 0;
			size_t m_final_planet_start_indice = 0;

			static inline bool s_initialized = false; // pour éviter de réinitialiser les positions des planètes à chaque fois que le constructeur est appelé, ce qui est coûteux en temps de calcul

			/* Pas implémenté pour le moment, mais pourrait être utile pour faire du multi-objective optimization ou pour faire du curriculum learning --> pas mon commentaire, se renseigner
			* Intéressant à faire pour que la gestion du score soit plus claire et plus facile à faire évoluer.
			struct ScoreComponents {
				double distance_to_final_planet = 0;
				double speed_at_final_planet = 0;
				double fuel_consumed = 0;
				double time_taken = 0;
			};
			*/
			

			static PlanetInfo s_start_planet_info;
			static PlanetInfo s_final_planet_info;

			constexpr static double m_max_acceleration = 5 * constants::g;

			double m_SolarSystemBound = 25; // en UA


		public:
			using Real = double;
			inline static constexpr double m_LowestScore = -1e9;
			inline static constexpr double m_CstScore = 1e2;

			inline static constexpr uint8_t m_NbPlanets = 9; // soleil + 8 planètes
			static const std::array<Planet, m_NbPlanets> m_planets;

			const enum class RocketState : uint8_t {
				DEAD_TOUCH_START_PLANET_HIGH_SPEED,
				DEAD_TOUCH_START_PLANET_LOW_SPEED,
				DEAD_TOUCH_SUN_HIGH_SPEED,
				DEAD_TOUCH_SUN_LOW_SPEED,
				DEAD_ACCELERATION_TOO_HIGH,
				DEAD_GET_TOO_FAR,
				DEAD_TOUCH_FINAL_PLANET_HIGH_SPEED,
				DEAD_TOUCH_FINAL_PLANET_LOW_SPEED,
				NEUTRAL,
				VALID
			};

			AdaptedSystem();

			/// <summary>
			/// constructeur de la classe AdaptedSystem
			/// </summary>
			/// <param name="start_planet"></param>
			/// <param name="final_planet"></param>
			/// <param name="rocket"></param>
			/// <param name="max_duration"> en jours </param>
			/// <param name="dt_seconds"> en secondes </param>
			AdaptedSystem(PlanetsName start_planet, PlanetsName final_planet, Rocket rocket, double max_duration, double dt_seconds);

			AdaptedSystem(const AdaptedSystem& sys);

			AdaptedSystem& operator=(const AdaptedSystem& sys) = delete;

			/// <summary>
			/// Effectue l'initialisation nécessaire pour préparer le composant ou le système à l'utilisation.
			/// </summary>
			void Initialize();

			/**
			* @brief Simule le système à partir de son état courant
			*
			* @param state Une fonction booléenne, prenant en paramètre une fusée, et indiquant si son état est valide
			* @return l'état de la fusée à la fin de la simulation
			*/
			RocketState Run(std::function<RocketState()> state); // TODO : Revoir la fonction de simulation pour modification

			/// <summary>
			/// Réinitialise le système à l'état à un instant donné, qui deviendra l'instant initial de la simulation.
			/// Ainsi, à la première mise à jour, le système sera à l'instant t, et m_time vaudra 0.
			/// </summary>
			/// <param name="t">l'instant auquel le système sera réinitialisé (en jours) </param>
			void Reset(double t);

			Real Score(const std::vector<std::vector<Real>>& genome);

			/// <summary>
			/// La taille de l'anneau dans lequel les fusées sont générées
			/// </summary>
			/// <returns>en km</returns>
			Real RingSize_meter() const noexcept;


			///////////////////////////////////////////////////////
			// Getters et Setters
			///////////////////////////////////////////////////////


			const Planet& getStartPlanet() const noexcept { return m_planets[static_cast<size_t>(m_start_planet)]; }
			const Planet& getFinalPlanet() const noexcept { return m_planets[static_cast<size_t>(m_final_planet)]; }
			const Planet& getSun() const noexcept  { return m_planets[0]; }

			/// <summary>
			/// Renvoie une référence constante au delta_time
			/// </summary>
			/// <returns> en secondes </returns>
			const double& getDeltaTime() const noexcept { return s_deltaTime; }

			/// <summary>
			/// Renvoie une référence constante à la durée maximale de la simulation.
			/// </summary>
			/// <returns> en jours </returns>
			const double& getMaxTime() const noexcept { return s_MaxTime; }

			/// <summary>
			/// Renvoie une référence constante à la date courante de la simulation.
			/// </summary>
			/// <returns> en jours </returns>
			const double& getCurrentTime() const noexcept { return m_time; }

			
			/// <summary>
			/// Renvoie la vitesse angulaire de la planète de départ.
			/// </summary>
			/// <returns>rad/s</returns>
			double getStartPlanetAngularVelocity() const noexcept { return s_start_planet_info.angular_velocity; }

			/// <summary>
			/// Renvoie la vitesse angulaire de la planète d'arrivée.
			/// </summary>
			/// <returns>rad/s</returns>
			double getFinalPlanetAngularVelocity() const noexcept { return s_final_planet_info.angular_velocity; }

			/// <summary>
			/// Retourne une référence constante vers le vecteur contenant les positions de la planète de départ.
			/// </summary>
			/// <returns>UA</returns>
			static const std::vector<vecteur>& getStartPlanetPositions() { return s_startPlanet_positions; }

			/// <summary>
			/// Retourne une référence constante vers le vecteur contenant les positions de la planète d'arrivée.
			/// </summary>
			/// <returns>UA</returns>
			static const std::vector<vecteur>& getFinalPlanetPositions() { return s_finalPlanet_positions; }

			static const PlanetInfo& getStartPlanetInfo() noexcept { return s_start_planet_info; }
			static const PlanetInfo& getFinalPlanetInfo() noexcept { return s_final_planet_info; }

			size_t getStartPlanetStartIndice() const noexcept { return m_start_planet_start_indice; }
			size_t getFinalPlanetStartIndice() const noexcept { return m_final_planet_start_indice; }

			size_t getStartPlanetPositionIndice() const noexcept { return static_cast<size_t>(daysInSeconds(m_time) / s_deltaTime) % s_start_planet_info.nb_iterations_orbit; }
			size_t getFinalPlanetPositionIndice() const noexcept { return static_cast<size_t>(daysInSeconds(m_time) / s_deltaTime) % s_final_planet_info.nb_iterations_orbit; }

			static bool IsInitialized() noexcept { return s_initialized; }

			void SetRocket(Rocket rocket) { m_rocket = rocket; }
			
			///////// other

			size_t getMaxIterations() const { return static_cast<size_t>(s_MaxTime / s_deltaTime); }

			RocketState Rocket_state() const;

			const char* TypeOfTrajectory(double score) const;
			const char* TypeOfTrajectory(RocketState state) const;

			void GetRocketTrajectory(std::vector<glm::dvec3>&);

			friend RocketState GetRocketState(const Rocket& rocket, AdaptedSystem& system);

		private:
			void InitPlanet(bool is_start_planet, PlanetsName name); // TODO AMELIORER L'IMPLEM
					

			bool rocket_collide_with(ObjectName name) const;


			///////////////////////////////////////////////////////
			// Physics
			///////////////////////////////////////////////////////

			/// <summary>
			/// Calcul la trajectoire d'un astre à une distance d du soleil, à une vitesse angulaire w
			/// </summary>
			/// <param name="w">rad/s</param>
			/// <param name="d">UA</param>
			/// <param name="roll">the angle around the x-axis (rad) </param>
			/// <param name="picth">the angle around the y-axis (rad)</param>
			/// <param name="positions">le tableau dans lequel écrire les positions</param>
			void Calculate_planet_trajectory(
				double w, double d,
				double roll, double pitch,
				std::vector<vecteur>* positions,
				const PlanetInfo& planet_info) const;


			///////////////////////////////////////////////////////
			// Fonctions de score
			///////////////////////////////////////////////////////

			Real HandleScoreValidState() const;
			Real HandleScoreNeutralState() const;
			Real HandleScoreInvalidGenerationState(GenerationState gen_state) const;

			/*
			* attractor_pos : position de l'attracteur (en UA)
			* attractor_mu : paramètre gravitationnel de l'attracteur (en m^3/s²)
			* object_pos : position d'un objet attiré gravitationnellement par l'attracteur (en UA)
			* object_mass : en kg
			*
			* Renvoi des kN (kg*km/s²)
			*/
			glm::dvec3 ComputeAttractionForce(glm::dvec3 attractor_pos, double attractor_mu,
				glm::dvec3 object_pos, double object_mass) const;

		}; // class AdaptedSystem


		AdaptedSystem::RocketState GetRocketState(const Structures::Rocket& rocket, AdaptedSystem& system);

	}; // namespace Systems


















	template <typename Real>
	std::pair<SimuCore::Structures::Rocket, GenerationState> IndividualToRocket(
		const std::vector<std::vector<Real>>& genome,
		SimuCore::Systems::AdaptedSystem& system
	) {
		// vérifier que Real est bien un type flottant
		static_assert(std::is_floating_point_v<Real>,
			"Real type must be a floating-point type (float, double, long double)");

		if (!Systems::AdaptedSystem::IsInitialized()) {
			std::runtime_error("Le système doit être initialisé avant de pouvoir convertir un genome en fusée.");
			std::abort();
		}

		SimuCore::Structures::Rocket rocket(0, std::vector<std::pair<SimuCore::Structures::Impulsion, double>>(), 700000, 4.5);

		// Comme les individus sont codés dans des vecteurs de coordonnées comprise entre 
		// [-system.RingSize_meter(), system.RingSize_meter()], on créer une fonction pour convertir
		//  cet intervalle en un intervalle [min, max]

		const double max_reels_genes = system.RingSize_meter(); // en km

		auto convert = [max_reels_genes](double min, double max) {
			return [max_reels_genes, min, max](double x) {
				return convertIntervals(-max_reels_genes, max_reels_genes, min, max, x);
				};
			};
		// la fonction convert est une fonction qui :
		// - prend en paramètre min et max
		// - retourne une fonction qui prend en paramètre x
		//		et qui convertit x de [-max_reels_genes, max_reels_genes] en [min, max]
		//
		// Autrement dit, elle sert à convertir les gènes de l'individu en des valeurs réelles comprises entre min et max.

		// ramène un gène dans [0, 1] |   norm01 = convert(0, 1) ?
		/*auto norm01 = [max_reels_genes](double x) {
			return (x + max_reels_genes) / (2.0 * max_reels_genes);
			};
			*/
		auto norm01 = convert(0, 1);

		///  --> à vérifier (Ok ?)
		auto convert_to_sphere_radius_1 = [norm01](
			const glm::dvec3& gene
			) {
				
				double u = norm01(gene.x); // ~~> correspond à r
				double v = norm01(gene.y); // ~~> correspond à theta
				double w = norm01(gene.z); // ~~> correspond à phi

				double r = std::cbrt(u);
				double costh = 2.0 * v - 1.0;
				double sinth = std::sqrt(1.0 - costh * costh);
				double phi = 2.0 * constants::PI * w;

				return glm::dvec3(
					r * sinth * std::cos(phi),
					r * sinth * std::sin(phi),
					r * costh
				);
			};
		

		auto convert_to_circle_radius_1 = [norm01](
			const glm::dvec2& gene
			) -> glm::dvec2 {
				// gene.x correspond à r et est généré uniformément
				// gene.y correspond à theta est est généré uniformément

				const double test = norm01(gene.x);
				const double r = std::sqrt(test);
				const double theta = 2 * constants::PI * norm01(gene.y);

				return glm::dvec2(
					r * std::cos(theta),
					r * std::sin(theta)
				);
			};


		/**

		Comme on force la première impulsion à t=0s, mais que l'on souhaite quand même pouvoir explorer différentes configurations
		initiales, on code les temps des impulsions suivantes comme des écarts de temps entre chaque impulsion.
		De plus, cela implique qu'il faut faire en sorte que le système soit réinitialisé à la bonne configuration initiale.
		**/

		// On réinitialise le système à l'intant t1 ("t1 devient 0").

		// pour cela, on met le temps t1 entre 0 et <le temps qu'il faut aux deux planètes pour revenir à la même configuration (on reste en 2D pour le moment)>
		double temps_1ere_impulsion = -1; // en jours
		{
			const SimuCore::Structures::Planet& sun = system.getSun();
			double omega_planete_initiale = system.getStartPlanetAngularVelocity();		// rad/s
			double omega_planete_finale = system.getFinalPlanetAngularVelocity();		// rad/s
			double temps_1ere_impulsion_max = std::abs(2 * constants::PI / (omega_planete_initiale - omega_planete_finale)); // en secondes

			// on veut remplacer l'intervalle [min_real, max_real] en [0, temps_1ere_impulsion_max]
			temps_1ere_impulsion = convert(0, seconds_to_days(temps_1ere_impulsion_max))(genome[2][0]);
		}

		system.Reset(temps_1ere_impulsion);











		/**		gestion de la position initiale qui va conditionner la transformation des vecteurs en impulsions
			En effet, on souhaite que la vitesse initiale (i.e la première impulsion) soit
			de norme supérieure à la vitesse de libération de la planète de départ.
			Et cette dernière dépend de la position de la fusée, qui est la position initiale,
			à l'instant initial.
		*/

		// On récupère la position initiale encore codée dans l'individu, qui est donc différente de
		// la position initiale réelle.

		
		glm::dvec2 position_initiale = convert_to_circle_radius_1( // en km
			glm::dvec2( 
				genome[0][0],
				genome[0][1]
			)
		);

		// Pour des raisons de non-priorisation de certaines directions, on impose que la position initiale soit dans un anneau (3D)
		//  (dont le rayon minimal est l'altitude de l'exosphère de la planète de départ
		//	 et dont le rayon maximal est l'altitude maximale autorisée par la planète de départ).

		{ // check pour enlever les cas dégénérés
			double norme_position_initiale = glm::length(position_initiale);
			if (norme_position_initiale < SimuCore::constants::epsilon
				|| norme_position_initiale > 1 // n'arrive normalement jamais
				) {
				// On retourne une fusée "invalide" et un indicateur pour annoncer que l'individu est invalide car la position initiale
				// n'est pas dans la sphère.
				if (norme_position_initiale > 1) { std::abort(); }

				return { rocket, GenerationState::INVALID_GENES_OUT_OF_BOUNDS_POSITION };
			}
		}



		const auto& start_planet = system.getStartPlanet();

		// On souhaite que les futurs positions soient dans un anneau de largeur maxOrbitRadius - minOrbitRadius.
		// Pour mettre la fusée sur l'anneau, la position initiale donnée doit être ajoutée avec 
		//	un vecteur colinéaire à lui même, de norme égale au rayon minimal de l'anneau.
		// On fera attention à mettre la bonne largeur dès le début

		// transformation 1er vecteur en position initiale
		{
			glm::dvec2 direction = glm::normalize(position_initiale);
			position_initiale *= (start_planet.maxOrbitRadius() - start_planet.minOrbitRadius());
			position_initiale += (start_planet.minOrbitRadius() * direction);
		}

		// A cet instant, la position initiale est dans l'anneau, mais elle n'est pas encore centrée sur la planète de départ.
		// On profite de ce moment pour faire la transformation sur la vitesse initiale, qui a besoin de la distance entre la
		// planète de départ et la position initiale pour calculer la vitesse de libération, et ainsi s'assurer que la vitesse
		// initiale soit dans le bon anneau.
		// Si on faisait la transformation de la vitesse initiale après avoir centré la position initiale sur la planète de
		// départ, on aurait un problème : le vecteur position_initiale serait la position de la fusée par rapport au soleil
		// et non pas par rapport à la planète de départ.


		/**		gestion de la vitesse initiale (1ère impulsion)
			La vitesse initiale doit être de norme supérieure à la vitesse de libération
			de la planète de départ, pour que la fusée puisse s'échapper de l'attraction de la planète.
			Et ainsi, on peut au moins espérer atteindre la planète d'arrivée.

			De plus, pour les mêmes raisons que pour la position initiale, on souhaite que la vitesse initiale soit
			distribuée uniformément dans l'espace (3D ou 2D).
		**/

		glm::dvec2 vitesse_initiale = convert_to_circle_radius_1( // km/s
			glm::dvec2(
				genome[1][0],
				genome[1][1]
			)
		);

		{ // traitement des cas dégénérés
			double norme_vitesse_initiale = glm::length(vitesse_initiale);
			if (norme_vitesse_initiale < SimuCore::constants::epsilon
				|| norme_vitesse_initiale > 1) {

				// On retourne une fusée "invalide" et un indicateur pour annoncer que l'individu est invalide car la vitesse initiale
				// n'est pas dans la sphère.
				if (norme_vitesse_initiale > 1) { std::abort(); }

				return { rocket, GenerationState::INVALID_GENES_OUT_OF_BOUNDS_VELOCITY };
			}
		}
		

		// On s'assure que la vitesse initiale soit dans un anneau (3D ou 2D) dont le rayon minimal
		// est la vitesse de libération de la planète de départ.

		// vesc​(r) ≤∥v0​∥≤ 2vesc​(r)
		{
			glm::dvec2 direction = glm::normalize(vitesse_initiale);

			const double v_liberation = start_planet.extractionVelocity(glm::length(position_initiale)); // km/s
			vitesse_initiale *= v_liberation;
			vitesse_initiale += (direction * v_liberation);
		}

		// maintenant que la vitesse initiale est dans le bon anneau, on peut la mettre dans la fusée.
		// Cepedant, la vitesse initiale va être ajouté grâce à la première impulsion de la fusée.
		// Nous allons donc simplement préparer cette impulsion, puis elle sera appliquée lors du lancement de la fusée.

		// Toutefois, il faut penser à ajouter la vitesse de la planète de départ pour avoir la vitesse dans le référentiel héliocentrique.
		// De plus, il faut aussi prendre en compte le fait que la fusée est supposée en orbite circulaire à l'instant 0- TODO

		const glm::dvec3& start_planet_initial_position = // en UA
			system.getStartPlanetPositions()[system.getStartPlanetStartIndice()];

		const glm::dvec3 start_planet_initial_velocity = glm::length(start_planet.velocity) * glm::normalize( // km/s
			glm::dvec3(
				- start_planet_initial_position[1],
				start_planet_initial_position [0],
				0 // 2D
			)
		);

		glm::dvec3 orbit_velocity = glm::dvec3(0); // km/s

		{
			// on récupère la position de la fusée autour de la planète de départ.
			// on en déduit le vecteur u_theta
			const glm::dvec3 u_theta = glm::normalize(glm::dvec3(
					- position_initiale.y,
					position_initiale.x,
					0
				)
			);

			// d'après la cahier de TIPE, on connait la norme de la vitesse d'un objet en orbite circulaire autour d'un astre.
			const double norme = meters_per_seconds_to_kilometers_per_seconds(1) * std::sqrt( // en km/s
				start_planet.getMu() / kilometers_to_meters(
					glm::length(position_initiale)
				)
			);

			orbit_velocity = norme * u_theta;
		}

		rocket.velocity = start_planet_initial_velocity + orbit_velocity; // on ajoute la vitesse de la planète de départ pour avoir la vitesse absolue. TODO


		// maintenant que la vitesse initiale est dans son anneau, on peut finaliser la position initiale de la fusée.
		// il suffit de translater l'anneau autour de la planète de départ
		position_initiale += (static_cast<double>(1.0_AU_to_km) * glm::dvec2(
			start_planet_initial_position.x,
			start_planet_initial_position.y
		)); // anneau centré sur la planète initiale

		// le calcul de la position initiale est terminé. On peut la mettre dans la fusée.
		rocket.position = (static_cast<double>(1.0_km_to_AU) * glm::dvec3(
			position_initiale.x,
			position_initiale.y,
			0
		));





		size_t nombre_impulsions = (genome.size() - 1) / 2; // c.f. cahier TIPE (2 vecteurs par impulsion + 1 vecteur pour p0)
		std::vector<std::pair<Structures::Impulsion, double>> toutes_les_impulsions; // km/s et en jours
		toutes_les_impulsions.reserve(nombre_impulsions);
		toutes_les_impulsions.emplace_back(Structures::Impulsion(glm::dvec3(
			vitesse_initiale.x,
			vitesse_initiale.y,
			0
		)), constants::epsilon); // la première impulsion est toujours appliquée à t=0s. On met un temps très proche de 0 pour éviter les problèmes numériques.

		


		// On peut maintenant traiter et ajouter les autres impulsions
		// Les individus sont composés de vecteurs générés aléatoirement dans un pavé droit, où chaque composante
		// est comprise entre [min_real, max_real].
		// Or on souhaite que les impulsions de la fusée soient de "direction aléatoire uniforme" (sphère != pavé droit)
		// Il va donc falloir utiliser la même méthode que précédemment pour la vitesse et la position.
		//
		// De plus, d'une part, une impulsion est la donnée d'un vecteur vitesse à ajouter au vecteur vitesse de notre fusée, 
		// et de l'instant auquel on souhaite que cette impulsion soit effectuée.
		// et d'autre part, les impulsions sont ordonnées.
		// Cela nous a poussé à ce que l'individu ne stocke pas directement les impulsions, mais des écarts de temps entre chaque impulsion.
		//
		// Cependant, à cause de cette contrainte d'ordre, il est possible que la somme des écarts de temps soit plus grande que le temps 
		// de simulation maximal autorisé.
		// Pour pallier ce problème, on va changer la signification des écarts de temps :
		// au lieu d'être des écarts de temps absolus, ils seront des écarts de temps "relatifs", c'est-à-dire
		// qu'on convertira l'écart de temps i, noté dt_i, en une nouvelle valeur, notée dt'_i, telle que :
		//		dt'_i = convert(0, m_MaxTime - somme_{j=1}^{i-1} dt'_j)(dt_i)
		// Ainsi, on s'assure que la somme des écarts de temps est toujours inférieure à m_MaxTime, que 
		// l'on garde le bon nombre d'impulsions, et que l'ordre des impulsions est respecté.
		// 
		// Cela permet donc de ne pas à renvoyer une fusée invalide si la somme des écarts de temps dépasse m_MaxTime.
		// Et cela a pour conséquence directe que la fonction de score sera beacoup moins plate, et donc plus facile à optimiser.
		// En effet, avant cette modification, si la somme des écarts de temps dépassait m_MaxTime,
		// la fusée était invalide, et donc le score était toujours le même (score très faible).


		{
			double somme_des_temps = 0; // ~~> somme_{j=1}^{i-1} dt'_j   |   (en jours)

			for (size_t i = 1; i < nombre_impulsions; i++) {
				// genome[2 * i + 2][0] correspond à l'écart de temps entre l'impulsion i et i-1
				// on remet instant entre 0 et (m_MaxTime - somme_des_temps)
				double instant = convert(0, system.getMaxTime() - somme_des_temps)(genome[2 * i + 2][0]) + constants::epsilon; // ~~> dt'_i et on ajoute une petite valeur epsilon pour éviter les problèmes numériques (en jours)
				somme_des_temps += instant; // mise à jour de la somme des temps

				glm::dvec2 impulsion = convert_to_circle_radius_1( // en km/s
					glm::dvec2(
						genome[2 * i + 1][0],
						genome[2 * i + 1][1]
					)
				); // en km/s

				// on s'assure que l'impulsion impulsion est généré dans la sphère :
				if (glm::length(impulsion) > 1) {
					// On retourne une fusée "invalide" et un indicateur pour annoncer que l'individu est invalide car une impulsion
					// n'est pas dans la sphère.

					return { rocket, GenerationState::INVALID_GENES_OUT_OF_BOUNDS_IMPULSION };
				}

				{
					double r0 = AU_to_meters(glm::length(start_planet_initial_position)); // m
					double v_orb = std::sqrt(constants::G * system.getSun().mass / r0); // m/s

					double alpha = std::clamp(0.3 / nombre_impulsions, 0.02, 0.15);

					/*
					Plus le coef alpha est grand, plus les impulsions pourront être forte. Ainsi la formulate reflète :
						des manoeuvres grosses mais rares, ou petites mais nombreuses.
					*/


					double dv_max = alpha * meters_per_seconds_to_kilometers_per_seconds(v_orb); // km/s

					impulsion *= dv_max; /*TODO : vérifier qu'Homman est faisable*/
				}

				const auto& [impul_precedente, temps_precedent] = toutes_les_impulsions[i - 1]; // impul_prec en km/s et temps_prec en jours
				instant += temps_precedent; // on convertit la variable instant en temps absolu

				toutes_les_impulsions.emplace_back(Structures::Impulsion(glm::dvec3(
					impulsion.x,
					impulsion.y,
					0
				)), instant);
			}
		}

		rocket.setImpulsions(toutes_les_impulsions);
		rocket.acceleration = 0; // une acceleration valide pour que le RocketState initial soit sûr d'être valide. (Pas utiliser dans les calculs)
		
		return { rocket, GenerationState::VALID };
	}
}; // namespace SimuCore
#pragma once


#include <pch.h>

#include "Planet.h"
#include "Entity.h"
#include "Rocket.h"





namespace SimuCore {
	const enum class GenerationState : uint8_t {
		VALID,
		INVALID_GENES_OUT_OF_BOUNDS_POSITION,
		INVALID_GENES_OUT_OF_BOUNDS_VELOCITY,
		INVALID_GENES_OUT_OF_BOUNDS_IMPULSION
	};

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

		SimuCore::Structures::Planet getPlanetFromName(PlanetsName name);

		class AdaptedSystem {
			using Entity = SimuCore::Structures::Entity;
			using Planet = SimuCore::Structures::Planet;
			using Rocket = SimuCore::Structures::Rocket;


			Planet m_startPlanet;
			Planet m_finalPlanet;

			Rocket m_rocket;
			Planet m_sun;

			double m_time = 0;

			double m_startAngle;
			double m_finalAngle;

			const double m_MaxTime;
			const double m_deltaTime;

			using vecteur = glm::dvec3;

			const enum class ObjectName : uint8_t {
				SUN,
				START,
				FINAL
			};


		public:
			using Real = double;
			inline static constexpr double m_LowestScore = -1e9;
			inline static constexpr double m_CstScore = 1e8;

			inline static constexpr uint8_t m_NbPlanets = 9;
			static const std::array<Planet, m_NbPlanets> m_planets; // fusée + soleil + planète depart + planète arrivée

			const enum class RocketState : uint8_t {
				VALID,
				NEUTRAL,
				DEAD_TOUCH_PLANET_HIGH_SPEED,
				DEAD_TOUCH_PLANET_LOW_SPEED,
				DEAD_ACCELERATION_TOO_HIGH
			};

			AdaptedSystem();

			/// <summary>
			/// constructeur de la classe AdaptedSystem
			/// </summary>
			/// <param name="start_planet"></param>
			/// <param name="final_planet"></param>
			/// <param name="start_angle"></param>
			/// <param name="final_angle"></param>
			/// <param name="rocket"></param>
			/// <param name="max_duration"> en jours </param>
			/// <param name="dt_seconds"> en secondes </param>
			AdaptedSystem(PlanetsName start_planet, PlanetsName final_planet, double start_angle, double final_angle, Rocket rocket, double max_duration, double dt_seconds);

			AdaptedSystem(const AdaptedSystem& sys);

			AdaptedSystem& operator=(const AdaptedSystem& sys) = delete;

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
			/// <param name="t">l'instant auquel le système sera réinitialisé</param>
			void Reset(double t);

			Real Score(const std::vector<std::vector<Real>>& individu);

			Real RingSize_meter() const;


			///////////////////////////////////////////////////////
			// Getters et Setters
			///////////////////////////////////////////////////////

			const Planet& getStartPlanet() const { return m_startPlanet; }
			const Planet& getFinalPlanet() const { return m_finalPlanet; }
			const Planet& getSun() const { return m_sun; }
			const double& getDeltaTime() const { return m_deltaTime; }
			const double& getMaxTime() const { return m_MaxTime; }
			const double& getCurrentTime() const { return m_time; }
			
			///////// other

			size_t getMaxIterations() const { return static_cast<size_t>(m_MaxTime / m_deltaTime); }

			RocketState Rocket_state() const;

			friend RocketState GetRocketState(const Rocket& rocket, AdaptedSystem& system);


		private:
			void InitPlanet(bool is_start_planet, PlanetsName name); // TODO AMELIORER L'IMPLEM
			void SetAnglePlanet(Planet& planet, double theta);
			void SetPlanetSpeed(Planet& planet, double theta);
					

			bool rocket_collide_with(ObjectName name) const;


			///////////////////////////////////////////////////////
			// Fonctions de score
			///////////////////////////////////////////////////////

			Real HandleScoreValidState() const;
			Real HandleScoreNeutralState() const;
			Real HandleScoreInvalidGenerationState(GenerationState gen_state) const;

		}; // class AdaptedSystem


		AdaptedSystem::RocketState GetRocketState(const Structures::Rocket& rocket, AdaptedSystem& system);

	}; // namespace Systems


















	template <typename Real>
	std::pair<SimuCore::Structures::Rocket, GenerationState> IndividualToRocket(
		const std::vector<std::vector<Real>>& individu,
		SimuCore::Systems::AdaptedSystem& system
	) {
		// vérifier que Real est bien un type flottant
		static_assert(std::is_floating_point_v<Real>,
			"Real type must be a floating-point type (float, double, long double)");


		// faire la conversion

		SimuCore::Structures::Rocket rocket(0, std::vector<std::pair<SimuCore::Structures::Impulsion, double>>(), 700000, 4.5);

		// Comme les individus sont codés dans des vecteurs de coordonnées comprise entre 
		// [-system.RingSize_meter(), system.RingSize_meter()], on créer une fonction pour convertir
		//  cet intervalle en un intervalle [min, max]

		const double max_reels_genes = system.RingSize_meter();

		auto convert = [max_reels_genes](double min, double max) {
			return [max_reels_genes, min, max](double x) {
				return convertIntervals(-max_reels_genes, max_reels_genes, -max, max, x);
			};
		};

		// la fonction convert est une fonction qui :
		// - prend en paramètre min et max
		// - retourne une fonction qui prend en paramètre x
		//		et qui convertit x de [-max_reels_genes, max_reels_genes] en [min, max]
		//
		// Autrement dit, elle sert à convertir les gènes de l'individu en des valeurs réelles comprises entre min et max.


		/**		gestion de la position initiale qui va conditionner la transformation des vecteurs en impulsions
			En effet, on souhaite que la vitesse initiale (i.e la première impulsion) soit
			de norme supérieure à la vitesse de libération de la planète de départ.
			Et cette dernière dépend de la position de la fusée, qui est la position initiale,
			à l'instant initial.
		*/

		// On récupère la position initiale encore codée dans l'individu, qui est donc différente de
		// la position initiale réelle.
		glm::dvec3 position_initiale = glm::dvec3(
			individu[0][0],
			individu[0][1],
			//individu[0][2] // 3D
			0 // 2D
		);

		// Pour des raisons de non-priorisation de certaines directions, on impose que la position initiale soit dans un anneau (3D)
		//  (dont le rayon minimal est l'altitude de l'exosphère de la planète de départ
		//	 et dont le rayon maximal est déduit de l'intervalle dans lequel les coordonnées des gènes sont).

		// Comme les gènes sont "uniformément distribués" dans un cube centré en 0 de coté 2*max_reels_genes, afin de s'assurer
		// d'avoir une distribution uniforme dans l'anneau, on impose que la position initiale codée dans l'individu
		// soit dans une sphère. Comme le cube a un coté de 2*max_reels_genes, on va choisir un rayon de sphère
		// égal à max_reels_genes.
		if (glm::length(position_initiale) > max_reels_genes) {
			// On retourne une fusée "invalide" et un indicateur pour annoncer que l'individu est invalide car la position initiale
			// n'est pas dans la sphère.

			// TODO : gérer ça proprement
			return { rocket, GenerationState::INVALID_GENES_OUT_OF_BOUNDS_POSITION};
		}

		// Pour mettre la fusée sur l'anneau, la position initiale donnée doit être ajoutée avec 
		//	un vecteur colinéaire à lui même, de norme égale au rayon minimal de l'anneau.

		// transformation 1er vecteur en position initiale
		{
			// on récupère un vecteur colinéaire au vecteur position_initiale
			glm::dvec3 direction = glm::normalize(position_initiale);
			position_initiale += system.getStartPlanet().minOrbitRadius() * direction;
		}

		// maintenant que la position initiale est dans le bon anneau, il suffit de translater l'anneau autour de la planète de départ
		position_initiale += system.getStartPlanet().position; // anneau centré sur la planète initiale

		// le calcul de la position initiale est terminé. On peut la mettre dans la fusée.
		rocket.position = position_initiale;


		/**		gestion de la vitesse initiale (1ère impulsion)
			La vitesse initiale doit être de norme supérieure à la vitesse de libération
			de la planète de départ, pour que la fusée puisse s'échapper de l'attraction de la planète.
			Et ainsi, on peut au moins espérer atteindre la planète d'arrivée.

			De plus, pour les mêmes raisons que pour la position initiale, on souhaite que la vitesse initiale soit
			distribuée uniformément dans l'espace (3D ou 2D).
		**/

		glm::dvec3 vitesse_initiale = glm::dvec3(
			individu[1][0],
			individu[1][1],
			//individu[1][2] // 3D
			0 // 2D
		);

		if (glm::length(vitesse_initiale) > max_reels_genes) {
			// On retourne une fusée "invalide" et un indicateur pour annoncer que l'individu est invalide car la vitesse initiale
			// n'est pas dans la sphère.

			// TODO : gérer ça proprement
			return { rocket, GenerationState::INVALID_GENES_OUT_OF_BOUNDS_VELOCITY };
		}

		// On utilise le même algorithme que pour la position initiale, pour s'assurer que la vitesse initiale soit dans
		// un anneau (3D ou 2D) dont le rayon minimal est la vitesse de libération de la planète de départ.
		{
			glm::dvec3 direction = glm::normalize(vitesse_initiale);
			vitesse_initiale += system.getStartPlanet().extractionVelocity(glm::length(position_initiale)) * direction;
		}

		// maintenant que la vitesse initiale est dans le bon anneau, on peut la mettre dans la fusée.
		// Cepedant, la vitesse initiale va être ajouté grâce à la première impulsion de la fusée.
		// Nous allons donc simplement préparer cette impulsion, puis elle sera appliquée lors du lancement de la fusée.

		// Toutefois, il faut penser à ajouter la vitesse de la planète de départ pour avoir la vitesse absolue dans le système.
		rocket.velocity = system.getStartPlanet().velocity; // on ajoute la vitesse de la planète de départ pour avoir la vitesse absolue.


		size_t nombre_impulsions = (individu.size() - 1) / 2; // c.f. cahier TIPE (2 vecteurs par impulsion + 1 vecteur pour p0)
		std::vector<std::pair<Structures::Impulsion, double>> toutes_les_impulsions;
		toutes_les_impulsions.reserve(nombre_impulsions);
		toutes_les_impulsions.emplace_back(Structures::Impulsion(vitesse_initiale), constants::epsilon); // la première impulsion est toujours appliquée à t=0s. On met un temps très proche de 0 pour éviter les problèmes numériques.

		/**

		Comme on force la première impulsion à t=0s, mais que l'on souhaite quand même pouvoir explorer différentes configurations
		initiales, on code les temps des impulsions suivantes comme des écarts de temps entre chaque impulsion.
		De plus, cela implique qu'il faut faire en sorte que le système soit réinitialisé à la bonne configuration initiale.

		**/

		// On réinitialise le système à l'intant t1 ("t1 devient 0").

		// pour cela, on met le temps t1 entre 0 et <le temps qu'il faut aux deux planètes pour revenir à la même configuration (on reste en 2D pour le moment)>
		double temps_1ere_impulsion = individu[2][0];
		{
			const SimuCore::Structures::Planet& sun = system.getSun();
			double omega_planete_initiale = system.getStartPlanet().getAngularVelocity(sun);
			double omega_planete_finale = system.getFinalPlanet().getAngularVelocity(sun);;
			double temps_1ere_impulsion_max = std::abs(2 * constants::PI / (omega_planete_initiale - omega_planete_finale));

			// on veut remplacer l'intervalle [min_real, max_real] en [0, temps_1ere_impulsion_max]
			temps_1ere_impulsion = convert(0, temps_1ere_impulsion_max)(temps_1ere_impulsion);
		}

		system.Reset(temps_1ere_impulsion);

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
			double somme_des_temps = 0; // ~~> somme_{j=1}^{i-1} dt'_j

			for (size_t i = 1; i < nombre_impulsions; i++) {
				// individu[2 * i + 2][0] correspond à l'écart de temps entre l'impulsion i et i-1
				// on remet instant entre 0 et (m_MaxTime - somme_des_temps)
				double instant = convert(0, system.getMaxTime() - somme_des_temps)(individu[2 * i + 2][0]) + constants::epsilon; // ~~> dt'_i et on ajoute une petite valeur epsilon pour éviter les problèmes numériques
				somme_des_temps += instant; // mise à jour de la somme des temps

				glm::dvec3 impulsion = glm::dvec3{
					individu[2 * i + 1][0],
					individu[2 * i + 1][1],
					//individu[2 * i + 1][2] // 3D
					0 // 2D
				};

				const auto& [impul_precedente, temps_precedent] = toutes_les_impulsions[i - 1];
				instant += temps_precedent; // on convertit la variable instant en temps absolu

				// on s'assure que l'impulsion impulsion est généré dans la sphère :
				if (glm::length(impulsion) > max_reels_genes) {
					// On retourne une fusée "invalide" et un indicateur pour annoncer que l'individu est invalide car une impulsion
					// n'est pas dans la sphère.

					// TODO : gérer ça proprement
					return { rocket, GenerationState::INVALID_GENES_OUT_OF_BOUNDS_IMPULSION };
				}

				toutes_les_impulsions.emplace_back(Structures::Impulsion(impulsion), instant);
			}
		}

		rocket.setImpulsions(toutes_les_impulsions);
		rocket.acceleration = 0; // une acceleration valide pour que le RocketState initial soit sûr d'être valide. (Pas utiliser dans les calculs)
		
		return { rocket, GenerationState::VALID };
	}
}; // namespace SimuCore




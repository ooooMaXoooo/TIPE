#include <pch.h>
#include <SimuCore/structures/System.h>

#include <SimuCore/constants.h>
#include <SimuCore/Units/all.h>

namespace SimuCore::Systems {
	static constexpr inline double convert_speed = static_cast<double>(1.0_m_per_s__to__km_per_s);

	std::array<SimuCore::Structures::Planet, AdaptedSystem::m_NbPlanets> const AdaptedSystem::m_planets = {
		// (Nom, Masse (kg), Rayon physique (km), Exobase (km), Alt. max (anneau)(km), Position initiale (UA), Vitesse initiale (km/s))

		// Le Soleil n'est pas une cible d'orbite finale. Rayon/Altitude = 0.
		SimuCore::Structures::Planet("Soleil",  1.989e30, 696340, 0, -696340),

		// Planètes avec Rayon Physique et Altitude Cible de 200 km (200e3 m)
		SimuCore::Structures::Planet("Mercure", 3.3011e23, 2439.7,    0,     -75, glm::dvec3(  57.91e9_m_to_AU, 0, 0), convert_speed * glm::dvec3(0, 47.87e3, 0)),
		SimuCore::Structures::Planet("Venus",   4.8675e24, 6051.8,  175,   10921, glm::dvec3( 108.21e9_m_to_AU, 0, 0), convert_speed * glm::dvec3(0, 35.02e3, 0)),
		SimuCore::Structures::Planet("Terre",   5.9722e24, 6371.0,  450,   19611, glm::dvec3(  149.6e9_m_to_AU, 0, 0), convert_speed * glm::dvec3(0, 29.78e3, 0)),
		SimuCore::Structures::Planet("Mars",    6.4171e23, 3389.5,  200,    9556, glm::dvec3( 227.92e9_m_to_AU, 0, 0), convert_speed * glm::dvec3(0, 24.13e3, 0)),
		SimuCore::Structures::Planet("Jupiter", 1.8982e27,  69911, 2000, 2332435, glm::dvec3( 778.57e9_m_to_AU, 0, 0), convert_speed * glm::dvec3(0, 13.07e3, 0)),
		SimuCore::Structures::Planet("Saturne", 5.6834e26,  58232, 3500, 2354450, glm::dvec3(1433.53e9_m_to_AU, 0, 0), convert_speed * glm::dvec3(0,  9.68e3, 0)),
		SimuCore::Structures::Planet("Uranus",  8.6810e25,  25362, 6000, 1875176, glm::dvec3(2872.46e9_m_to_AU, 0, 0), convert_speed * glm::dvec3(0,  6.80e3, 0)),
		SimuCore::Structures::Planet("Neptune", 1.0241e26,  24622, 3000, 3197294, glm::dvec3(4495.06e9_m_to_AU, 0, 0), convert_speed * glm::dvec3(0,  5.43e3, 0))
	};

	PlanetInfo AdaptedSystem::s_start_planet_info = PlanetInfo();
	PlanetInfo AdaptedSystem::s_final_planet_info = PlanetInfo();
	std::vector<AdaptedSystem::vecteur> AdaptedSystem::s_startPlanet_positions = std::vector<AdaptedSystem::vecteur>(); // UA
	std::vector<AdaptedSystem::vecteur> AdaptedSystem::s_finalPlanet_positions = std::vector<AdaptedSystem::vecteur>(); // UA




	AdaptedSystem::AdaptedSystem() :
		m_rocket(Rocket(365.25, std::vector<std::pair<SimuCore::Structures::Impulsion, double>>(), 200._ton_to_kg, 4.5)),
		m_time(0)
	{
	}

	AdaptedSystem::AdaptedSystem(
		PlanetsName start_planet,
		PlanetsName final_planet,
		Rocket rocket,
		double max_duration,
		double dt_seconds
	) :
		m_rocket(rocket),
		m_time(0)
	{
		if (s_initialized && (s_start_planet != start_planet || s_final_planet != final_planet)) {
			std::cerr << "Error creating a system !" << std::endl;
			std::abort();
		}

		s_start_planet = start_planet;
		s_final_planet = final_planet;

		s_MaxTime = max_duration;
		s_deltaTime = dt_seconds;

		InitPlanet(true, start_planet);
		InitPlanet(false, final_planet);

	} // AdaptedSystem

	AdaptedSystem::AdaptedSystem(const AdaptedSystem& sys) :
		m_rocket(sys.m_rocket), // TODO vérifier constructeur par copie
		m_time(sys.m_time)
	{
	}

	AdaptedSystem::RocketState AdaptedSystem::Run(std::function<RocketState()> state) {	
		// on simule jusqu'à m_MaxTime ou la mort de la fusée ou son succès
		const size_t MAX_ITERATIONS = static_cast<const size_t>(daysInSeconds(s_MaxTime) / s_deltaTime);
		//RocketState current_state = state();
		RocketState current_state = Rocket_state();
		size_t iteration = 0;
		size_t start_planet_index = m_start_planet_start_indice;
		size_t final_planet_index = m_final_planet_start_indice;

		while (current_state == RocketState::NEUTRAL && iteration < MAX_ITERATIONS) {
			// on calcul les forces gravitationnelles agissant sur la fusée

			m_rocket.forces = glm::dvec3(0);
			m_rocket.forces += ComputeAttractionForce(
				s_startPlanet_positions[start_planet_index],
				s_start_planet_info.muPlanet,
				m_rocket.position,
				m_rocket.mass
			);

			m_rocket.forces += ComputeAttractionForce(
				s_finalPlanet_positions[final_planet_index],
				s_final_planet_info.muPlanet,
				m_rocket.position,
				m_rocket.mass
			);

			m_rocket.forces += ComputeAttractionForce(glm::dvec3(0), getSun().getMu(), m_rocket.position, m_rocket.mass);

			m_rocket.UpdateFirstPart(s_deltaTime);

			// recalcul des forces pour l'intégrateur
			{
				m_rocket.forces += ComputeAttractionForce(
					s_startPlanet_positions[start_planet_index],
					s_start_planet_info.muPlanet,
					m_rocket.position,
					m_rocket.mass
				);

				m_rocket.forces += ComputeAttractionForce(
					s_finalPlanet_positions[final_planet_index],
					s_final_planet_info.muPlanet,
					m_rocket.position,
					m_rocket.mass
				);

				m_rocket.forces += ComputeAttractionForce(glm::dvec3(0), getSun().getMu(), m_rocket.position, m_rocket.mass);
			}

			glm::dvec3 velocity_before = m_rocket.velocity;
			m_rocket.ApplyImpulsions(m_time, s_deltaTime);
			if (glm::length(m_rocket.velocity - velocity_before) / s_deltaTime >= m_max_acceleration) {
				current_state = RocketState::DEAD_ACCELERATION_TOO_HIGH;
				break;
			}

			m_rocket.UpdateSecondPart(s_deltaTime);

			//current_state = state();
			current_state = Rocket_state();
			m_time += seconds_to_days(s_deltaTime);
			iteration++;
			start_planet_index = (start_planet_index + 1) % s_start_planet_info.nb_iterations_orbit;
			final_planet_index = (final_planet_index + 1) % s_final_planet_info.nb_iterations_orbit;
		}

		return current_state;
	} // Run

	void AdaptedSystem::Reset(double t) {
		if (!s_initialized) {
			std::abort();
		}

		m_time = 0;

		t = days_to_seconds(t);
		m_start_planet_start_indice = static_cast<size_t> (t / s_deltaTime) % s_start_planet_info.nb_iterations_orbit;
		m_final_planet_start_indice = static_cast<size_t>(t / s_deltaTime) % s_final_planet_info.nb_iterations_orbit;

		// besoin de reset la rocket ???  TODO


	} // Reset

	AdaptedSystem::Real AdaptedSystem::Score(const std::vector<std::vector<Real>>& genome) {
		auto [rocket, gen_state] = IndividualToRocket(genome, *this);

		if (gen_state != GenerationState::VALID) {
			return HandleScoreInvalidGenerationState(gen_state);
		}

		m_rocket = rocket;



		// on lance une simulation et on récupère :
		//	• la position finale
		//	• si l'accélération a été trop forte (--> c++ 23 pour std::expected) --> pas besoin, gestion de l'état de la fusée par une énumération RocketState
		//	• si le temps est trop long (--> c++ 23 pour std::expected) --> pas besoin, on arrête la simulation avant + gestion de l'état de la fusée par une énumération RocketState
		//	• la position finale de l'objectif
		//	• le coût énergétique
		//	• le temps de vol
		//

		auto rocket_state = [this]() -> RocketState {
			return Rocket_state();
		};

		RocketState final_state = Run(rocket_state);

		switch (final_state) {
		case RocketState::DEAD_TOUCH_START_PLANET_HIGH_SPEED:
			// on veut que la borne inf soit 0
			// on veut que la borne sup soit 1/m_CstScore
			return (1 / (m_CstScore + glm::length(m_rocket.velocity))); // TODO : ajuster la formule
			break;
		case RocketState::DEAD_TOUCH_START_PLANET_LOW_SPEED:
			// on veut que la borne inf soit 1/m_CstScore
			// on veut que la borne sup soit 2/m_CstScore
			return (1 / (m_CstScore + glm::length(m_rocket.velocity))) + (1 / m_CstScore); // TODO : ajuster la formule
			break;
		case RocketState::DEAD_TOUCH_SUN_HIGH_SPEED:
			// on veut que la borne inf soit 2/m_CstScore
			// on veut que la borne sup soit 3/m_CstScore
			return (1 / (m_CstScore + glm::length(m_rocket.velocity))) + (2 / m_CstScore); // TODO : ajuster la formule
			break;
		case RocketState::DEAD_TOUCH_SUN_LOW_SPEED:
			// on veut que la borne inf soit 3/m_CstScore
			// on veut que la borne sup soit 4/m_CstScore
			return (1 / (m_CstScore + glm::length(m_rocket.velocity))) + (3 / m_CstScore); // TODO : ajuster la formule
			break;
		case RocketState::DEAD_ACCELERATION_TOO_HIGH:
			// on veut que la borne inf soit 4/m_CstScore
			// on veut que la borne sup soit 5/m_CstScore
			return (1 / (m_CstScore + m_rocket.acceleration)) + (4 / m_CstScore); // TODO : ajuster la formule
			break;
		case RocketState::DEAD_GET_TOO_FAR:
			// on veut que la borne inf soit 5/m_CstScore
			// on veut que la borne sup soit 6/m_CstScore
			return (1 / (m_CstScore + (glm::length(m_rocket.position) * glm::length(m_rocket.velocity)))) + (5 / m_CstScore);
			break;

		case RocketState::DEAD_TOUCH_FINAL_PLANET_HIGH_SPEED:
			// on veut que la borne inf soit 6/m_CstScore
			// on veut que la borne sup soit 7/m_CstScore
			return (1 / (m_CstScore + glm::length(m_rocket.velocity))) + (6 / m_CstScore); // TODO : ajuster la formule
			break;
		case RocketState::DEAD_TOUCH_FINAL_PLANET_LOW_SPEED:
			// on veut que la borne inf soit 7/m_CstScore
			// on veut que la borne sup soit 8/m_CstScore
			return (1 / (m_CstScore + glm::length(m_rocket.velocity))) + (7 / m_CstScore); // TODO : ajuster la formule
			break;
		case RocketState::NEUTRAL:
			// on veut que la borne inf soit 8/m_CstScore
			// on veut que la borne sup soit 9/m_CstScore
			return HandleScoreNeutralState(); // TODO : ajuster la formule pour que le score soit dans l'intervalle souhaité
			break;

		case RocketState::VALID:
			// on veut que la borne inf soit 9/m_CstScore
			// on veut que la borne sup soit 10/m_CstScore
			return HandleScoreValidState(); // TODO : ajuster la formule pour que le score soit dans l'intervalle souhaité
			break;
		
		default:
			throw std::runtime_error("Unknown RocketState encountered in Score calculation.");
			break;
		}
	} // Score

	
	AdaptedSystem::Real AdaptedSystem::RingSize_meter() const noexcept {
		const auto& start_planet = getPlanetFromName(s_start_planet);
		return start_planet.maxOrbitRadius() - start_planet.minOrbitRadius();
	} // RingSize_meter

	void AdaptedSystem::Initialize()
	{
		if(s_initialized) {
			std::runtime_error("AdaptedSystem::Initialize() called more than once. This is not allowed because it would cause a significant increase in initialization time due to the recalculation of planet trajectories.");
			std::abort();
		}

		InitPlanet(true, s_start_planet);
		InitPlanet(false, s_final_planet);
		
		s_initialized = true;
	} // Initialize

	void AdaptedSystem::InitPlanet(bool is_start_planet, PlanetsName name) {

		const Planet& planet = getPlanetFromName(name);

		double distance = glm::length(planet.position); // en UA
		double omega = planet.getAngularVelocity(distance, getSun().getMu()) / 3600; // en rad/s


		PlanetInfo& planet_info = is_start_planet ? s_start_planet_info : s_final_planet_info;
		planet_info.distance_to_sun = distance;
		planet_info.angular_velocity = omega;
		planet_info.muPlanet = planet.getMu();
		planet_info.nb_iterations_orbit = static_cast<size_t>((2 * constants::PI) / (omega * s_deltaTime));
		
		Calculate_planet_trajectory(
			omega,
			distance,
			0,
			0,
			is_start_planet ? &s_startPlanet_positions : &s_finalPlanet_positions,
			is_start_planet ? s_start_planet_info : s_final_planet_info
		);
	} // InitPlanet


	void AdaptedSystem::Calculate_planet_trajectory(
		double w,
		double d,
		double roll,
		double pitch,
		std::vector<vecteur>* positions,
		const PlanetInfo& planet_info) const
	{
		/* on construit une base (du plan de la trajectoire) :
		*	aller voir le cahier TIPE pour les détails de la construction de la base
		*/
		const glm::dvec3 u = d * glm::dvec3(std::cos(pitch), 0, std::sin(pitch));
		const glm::dvec3 v = d * glm::dvec3(0, std::cos(roll), std::sin(roll));

		// on calcul le nombre de positions à calculer en fonction de la vitesse angulaire w et 
		// du pas de temps m_deltaTime
		const size_t nb_positions = planet_info.nb_iterations_orbit;
		const double dtheta = (2 * constants::PI) / nb_positions; // pas angulaire entre chaque position
		double theta = 0; // angle initial

		positions->resize(nb_positions); // éviter les reallocations

		/* on calcule la position de la planète à chaque instant */
		for (size_t n = 0; n < nb_positions; n++) {
			(*positions)[n] = (u * std::cos(theta) + v * std::sin(theta));
			theta += dtheta;
		}
	} // Calculate_planet_trajectory

	AdaptedSystem::Real AdaptedSystem::HandleScoreValidState() const
	{
		// On sait que la fusée est en orbite stable autour de la planète cible :
		// On peut donc renvoyer un score parfait ou presque (en fonction du coût énergétique et du temps)
		// Il reste à déterminer ce qu'est un score parfait. Il nous faut un meilleur score que dans le cas neutre. C'est à dire un majorant du score neutre.
		
		constexpr Real Majorant_etat_neutre = 9 / m_CstScore;


		Real cout_energetique = m_rocket.getDeltaM();
		Real tof = m_time;
		return 1/(tof*cout_energetique + m_CstScore) + Majorant_etat_neutre;
	} // HandleScoreValidState

	AdaptedSystem::Real AdaptedSystem::HandleScoreNeutralState() const
	{
		/*Real cout_energetique = m_rocket.getDeltaM();
		Real tof = m_time;

		Real influence_temps_energie = 0.5/(cout_energetique * tof + m_CstScore);
		
		--> On a supprimé l'influence du temps et de l'énergie du score neutre, car dans l'état neutre, 
		la fusée n'est pas encore en orbite stable autour de la planète cible. Donc on veut d'abord encourager
		la fusée à se rapprocher de la planète cible, avant de prendre en compte le temps et l'énergie consommée.
		
		*/



		Real influence_position = 0;

		Planet final_planet = getPlanetFromName(s_final_planet);

		const glm::dvec3 final_planet_position = s_finalPlanet_positions[getFinalPlanetPositionIndice()]; // en UA

		glm::dvec3 projection_on_ring = final_planet.orbitRadius() * glm::normalize(m_rocket.position - final_planet_position); // en km
		Real distance_to_ring = glm::length((m_rocket.position - final_planet_position) - (kilometers_to_AU(1) * projection_on_ring)); // en km
		distance_to_ring = AU_to_kilometers(distance_to_ring); // conversion en km

		
		bool position_ok = (final_planet.minOrbitRadius() <= distance_to_ring && distance_to_ring <= final_planet.maxOrbitRadius());
		if (!position_ok)
		{
			// passage en UA pour avoir de meilleurs nombres
			distance_to_ring = kilometers_to_AU(distance_to_ring); // conversion en UA

			// std::cbrt -> racine cubique

			influence_position = 0.75 / (std::cbrt(distance_to_ring) + m_CstScore); // On veut que l'influence diminue quand on s'éloigne de la distance cible
		}
		else {
			influence_position = 0.75 / m_CstScore; // On veut que l'influence soit maximale quand on est dans la zone cible, mais qu'elle reste inférieure à la borne inf du score neutre pour les positions hors de la zone cible

			// si on est dans ce cas, la fusée se trouve dans l'anneau, mais n'est pas en orbite stable.
			// On peut donc ajouter la vitesse finale dans l'influence. Mais on souhaite toujours que 
			// le score lorsqu'on est dans l'anneau soit supérieur au score lorsqu'on est hors de l'anneau.
			// on fait donc en sorte que l'influence de la vitesse soit au maximum de 0.25/m_CstScore, pour 
			// faire en sorte qu'on balaye des score dont l'étendue en 1/m_CstScore
			influence_position += 0.25 / (std::sqrt(glm::length(m_rocket.velocity)) + m_CstScore);
			// on prend la racine de la vitesse pour que l'influence de la vitesse soit grandit d'autant plus
			// que le vitesse finale est faible.
		}
		// on veut que la borne inf soit 8/epsilon, donc on ajoute 8 / epsilon
		// donc la borne sup est 1/epsilon + 8/epsilon = 9/epsilon

		// TODO : vérifier avec des assertions que les influences sont bien dans les bornes attendues.
		return influence_position + 8/m_CstScore;
	} // HandleScoreNeutralState

	AdaptedSystem::Real AdaptedSystem::HandleScoreInvalidGenerationState(GenerationState gen_state) const {
		// TODO : améliorer en fonction des différents états invalides, pour pénaliser plus ou moins sévèrement
		// afin d'orienter la génération vers des individus valides.

		switch (gen_state)
		{
		case SimuCore::GenerationState::VALID:
			std::cerr << "Warning: HandleScoreInvalidGenerationState called with VALID state." << std::endl;
			std::abort();
			break;
		case SimuCore::GenerationState::INVALID_GENES_OUT_OF_BOUNDS_POSITION:
			return m_LowestScore * 0.001; // Pour pouvoir différencier des autres états invalides
			break;
		case SimuCore::GenerationState::INVALID_GENES_OUT_OF_BOUNDS_VELOCITY:
			return m_LowestScore * 0.01; // Pour pouvoir différencier des autres états invalides
			break;
		case SimuCore::GenerationState::INVALID_GENES_OUT_OF_BOUNDS_IMPULSION:
			return m_LowestScore * 0.1; // Pour pouvoir différencier des autres états invalides
			break;
		default:
			break;
		}
	} // HandleScoreInvalidGenerationState

	glm::dvec3 AdaptedSystem::ComputeAttractionForce(glm::dvec3 attractor_pos, double attractor_mu, glm::dvec3 object_pos, double object_mass) const
	{
		// F = G * m1 * m2 / r^2 * direction

		glm::dvec3 direction = (attractor_pos - object_pos); // dirigé vers la planète attractrice (en UA)
		direction *= 1._AU_to_m; // conversion en mètres
		double distance_carre = glm::dot(direction, direction); // en m²
		direction = glm::normalize(direction); // direction unitaire (sans unité)


		return ((attractor_mu * object_mass / distance_carre) * direction) * 1e-3; // conversion en kN (kg*km/s²).
	}

	bool AdaptedSystem::rocket_collide_with(ObjectName name) const {
		glm::dvec3 object_pos; // en UA
		double object_radius; // en km

		switch (name) {
		case ObjectName::SUN :
			object_pos = glm::dvec3(0);
			object_radius = m_planets[0].getRadius();
			break;

		case ObjectName::START :
			object_pos = getPlanetFromName(s_start_planet).position;
			object_radius = getPlanetFromName(s_start_planet).getRadius();
			break;

		case ObjectName::FINAL :
			object_pos = getPlanetFromName(s_final_planet).position;
			object_radius = getPlanetFromName(s_final_planet).getRadius();
			break;

		default:
			throw std::runtime_error("Unknown ObjectName encountered in rocket_collide_with.");
			break;
		}

		double distance = glm::length(m_rocket.position - object_pos); // en UA

		return AU_to_kilometers(distance) < object_radius + constants::epsilon;
	}
	// rocket_collide_with

	SimuCore::Structures::Planet getPlanetFromName(PlanetsName name)
	{
		return SimuCore::Structures::Planet(Systems::AdaptedSystem::m_planets[static_cast<size_t>(name)]);
	} // getPlanetFromName

	AdaptedSystem::RocketState GetRocketState(const Structures::Rocket& rocket, AdaptedSystem& system) {

		// on sauvegarde l'état actuel du système
		Structures::Rocket saved_rocket = system.m_rocket;
		double saved_time = system.m_time;

		// on remplace la fusée par celle passée en paramètre et on calcule l'état
		system.m_rocket = rocket;
		AdaptedSystem::RocketState state = system.Rocket_state(); // Rocket_state est une fonction const donc il n'y avait pas besoin de sauvegarder/restaurer l'état du système, mais on le fait quand même pour être sûr

		// on restaure l'état du système
		system.m_rocket = saved_rocket;
		system.m_time = saved_time;

		return state;
	}

	AdaptedSystem::RocketState AdaptedSystem::Rocket_state() const {
		// on check la postion
		if (glm::length(m_rocket.position) > m_SolarSystemBound) {
			return RocketState::DEAD_GET_TOO_FAR;
		}


		// on check l'acceleration
		Real acceleration = m_rocket.acceleration * 1.0_km_to_m; // en m/s²
		if (acceleration > m_max_acceleration) {
			return RocketState::DEAD_ACCELERATION_TOO_HIGH;
		}

		// check les collisions --> brute force car seulement 3 comparaisons
		// TODO : Corriger les erreurs : des fonctions membre utilisent m_rocket au lieu de rocket passé en

		Planet start_planet = getPlanetFromName(s_start_planet);
		Planet final_planet = getPlanetFromName(s_final_planet);

		bool collided = rocket_collide_with(ObjectName::SUN) || rocket_collide_with(ObjectName::START) || rocket_collide_with(ObjectName::FINAL);

		if (collided) {
			// déterminer la vitesse au moment de la collision
			Real speed = 0; // en km/s

			// Il faut déterminer si la vitesse est "faible" ou "élevée"
			// On considère qu'une vitesse inférieure à la vitesse de libération de l'astre qu'on a touché est une "faible" vitesse
			Real vitesse_de_liberation;
			if (rocket_collide_with(ObjectName::SUN)) {
				vitesse_de_liberation = getSun().extractionVelocity(getSun().getRadius()); // en km/s
				speed = glm::length(m_rocket.velocity); // en km/s, vitesse relative au soleil

				if (speed < vitesse_de_liberation) {
					return RocketState::DEAD_TOUCH_SUN_LOW_SPEED;
				}
				else {
					return RocketState::DEAD_TOUCH_SUN_HIGH_SPEED;
				}

			}
			else if (rocket_collide_with(ObjectName::START)) {
				vitesse_de_liberation = start_planet.extractionVelocity(start_planet.getRadius()); // en km/s
				speed = glm::length(m_rocket.velocity - start_planet.velocity); // en km/s, vitesse relative à la planète de départ


				if (speed < vitesse_de_liberation) {
					return RocketState::DEAD_TOUCH_START_PLANET_LOW_SPEED;
				}
				else {
					return RocketState::DEAD_TOUCH_START_PLANET_HIGH_SPEED;
				}
			}
			else { // FINAL
				vitesse_de_liberation = final_planet.extractionVelocity(final_planet.getRadius()); // en km/s
				speed = glm::length(m_rocket.velocity - final_planet.velocity); // en km/s, vitesse relative à la planète d'arrivée
			
				if (speed < vitesse_de_liberation) {
					return RocketState::DEAD_TOUCH_FINAL_PLANET_LOW_SPEED;
				}
				else {
					return RocketState::DEAD_TOUCH_FINAL_PLANET_HIGH_SPEED;
				}
			}
		}

		// on check la position finale, et si la fusée est en orbite autour de la planète finale
		glm::dvec3 position_dans_referentiel_planete = m_rocket.position - final_planet.position; // en UA, position de la fusée dans le référentiel de la planète finale
		Real distance = glm::length(position_dans_referentiel_planete); // en UA, distance entre la fusée et la planète finale

		distance = AU_to_kilometers(distance); // conversion en km pour le test de distance_ok
		bool distance_ok = (final_planet.minOrbitRadius() <= distance && distance <= final_planet.maxOrbitRadius());

		if (distance_ok) {
			glm::dvec3 vecteur_vitesse_referentiel_planete_finale = m_rocket.velocity - final_planet.velocity; // en km/s, vitesse de la fusée dans le référentiel de la planète finale
			Real norme_vitesse = glm::length(vecteur_vitesse_referentiel_planete_finale); // en km/s

			/* CONVERSION EN USI */
			distance = kilometers_to_meters(distance); // en m
			norme_vitesse = kilometers_per_seconds_to_meters_per_seconds(norme_vitesse); // en m/s

			position_dans_referentiel_planete *= AU_to_meters(1); // en m
			vecteur_vitesse_referentiel_planete_finale *= kilometers_per_seconds_to_meters_per_seconds(1); // en m/s
			/* fin de la conversion */

			// Toute la suite des calculs se fait en USI

			// TODO check les formules
			double energie_orbitale = 0.5 * norme_vitesse * norme_vitesse - final_planet.getMu() / distance;

			if (energie_orbitale >= 0 - constants::epsilon) return RocketState::NEUTRAL; // vérifier qu'on enleve pas trop de trajectoires


			glm::dvec3 presque_moment_cinetique = glm::cross(position_dans_referentiel_planete, vecteur_vitesse_referentiel_planete_finale);
			double constante_des_aires = glm::length(presque_moment_cinetique);

			double demi_grand_axe = -0.5 * final_planet.getMu() / energie_orbitale;
			double excentricite = std::sqrt(1 + (2 * energie_orbitale * (constante_des_aires / final_planet.getMu()) * (constante_des_aires / final_planet.getMu())));

			double perige = demi_grand_axe * (1 - excentricite);
			double apogee = demi_grand_axe * (1 + excentricite);

			// on reste dans l'anneau si le perigé est plus grand que l'orbite min et que l'apogée est plus petit que l'orbite min
			// car si on dépasse, le soleil a une forte influences
			bool vitesse_ok = (perige >= final_planet.minOrbitRadius()) && (apogee <= final_planet.maxOrbitRadius());

			return vitesse_ok ? RocketState::VALID : RocketState::NEUTRAL;
		}

		else {
			return RocketState::NEUTRAL;
		}
	}
	const char* AdaptedSystem::TypeOfTrajectory(double score) const
	{
		if (score > - constants::epsilon) {
			//constexpr double cste = 1 / m_CstScore;
			score *= m_CstScore; // <=> /= cste

			const char* kinds[10] = {
				"Dead : Collision start panet with high speed",
				"Dead : Collision start panet with low speed",
				"Dead : Collision sun with high speed",
				"Dead : Collision sun with low speed",
				"Dead : Acceleration too high",
				"Dead : Rocket get too far in the system",
				"Dead : Collision final panet with high speed",
				"Dead : Collision final panet with low speed",
				"Alive : Neutral",
				"Alive : Valid",
			};

			return kinds[static_cast<int>(score)];
		}
		else {
			return "Undefined";
		}
	}
	const char* AdaptedSystem::TypeOfTrajectory(RocketState state) const
	{
		constexpr const char* kinds[10] = {
			"Dead : Collision start panet with high speed",
			"Dead : Collision start panet with low speed",
			"Dead : Collision sun with high speed",
			"Dead : Collision sun with low speed",
			"Dead : Acceleration too high",
			"Dead : Rocket get too far in the system",
			"Dead : Collision final panet with high speed",
			"Dead : Collision final panet with low speed",
			"Alive : Neutral",
			"Alive : Valid",
		};

		return kinds[static_cast<uint8_t>(state)];
	}

	void AdaptedSystem::GetRocketTrajectory(std::vector<glm::dvec3>& trajectory)
	{
		// on simule jusqu'à m_MaxTime ou la mort de la fusée ou son succès
		const size_t MAX_ITERATIONS = static_cast<const size_t>(daysInSeconds(s_MaxTime) / s_deltaTime);
		RocketState current_state = Rocket_state();
		size_t iteration = 0;
		size_t start_planet_index = m_start_planet_start_indice;
		size_t final_planet_index = m_final_planet_start_indice;

		trajectory.reserve(MAX_ITERATIONS + 1);

		while (current_state == RocketState::NEUTRAL && iteration < MAX_ITERATIONS) {
			// on calcul les forces gravitationnelles agissant sur la fusée
			trajectory.push_back(m_rocket.position);

			m_rocket.forces = glm::dvec3(0);
			m_rocket.forces += ComputeAttractionForce(
				s_startPlanet_positions[start_planet_index],
				s_start_planet_info.muPlanet,
				m_rocket.position,
				m_rocket.mass
			);

			m_rocket.forces += ComputeAttractionForce(
				s_finalPlanet_positions[final_planet_index],
				s_final_planet_info.muPlanet,
				m_rocket.position,
				m_rocket.mass
			);

			m_rocket.forces += ComputeAttractionForce(glm::dvec3(0), getSun().getMu(), m_rocket.position, m_rocket.mass);

			m_rocket.UpdateFirstPart(s_deltaTime);

			// recalcul des forces pour l'intégrateur
			{
				m_rocket.forces += ComputeAttractionForce(
					s_startPlanet_positions[start_planet_index],
					s_start_planet_info.muPlanet,
					m_rocket.position,
					m_rocket.mass
				);

				m_rocket.forces += ComputeAttractionForce(
					s_finalPlanet_positions[final_planet_index],
					s_final_planet_info.muPlanet,
					m_rocket.position,
					m_rocket.mass
				);

				m_rocket.forces += ComputeAttractionForce(glm::dvec3(0), getSun().getMu(), m_rocket.position, m_rocket.mass);
			}

			glm::dvec3 velocity_before = m_rocket.velocity;
			m_rocket.ApplyImpulsions(m_time, s_deltaTime);
			if (glm::length(m_rocket.velocity - velocity_before) / s_deltaTime >= m_max_acceleration) {
				current_state = RocketState::DEAD_ACCELERATION_TOO_HIGH;
				break;
			}

			m_rocket.UpdateSecondPart(s_deltaTime);

			current_state = Rocket_state();
			iteration++;
			m_time += seconds_to_days(s_deltaTime);
			start_planet_index = (start_planet_index + 1) % s_start_planet_info.nb_iterations_orbit;
			final_planet_index = (final_planet_index + 1) % s_final_planet_info.nb_iterations_orbit;

		}

		trajectory.push_back(m_rocket.position);
	}
}; // namespace SimuCore::Systems
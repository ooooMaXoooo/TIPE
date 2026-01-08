#include <pch.h>
#include <SimuCore/structures/System.h>

#include <SimuCore/integrator/integrator.h>
#include <SimuCore/constants.h>

namespace SimuCore::Systems {

	std::array<SimuCore::Structures::Planet, AdaptedSystem::m_NbPlanets> const AdaptedSystem::m_planets = {
		// (Nom, Masse (kg), Rayon physique (km), Exobase (km), Alt. max (anneau)(km), Position initiale (m), Vitesse initiale (m/s))

		// Le Soleil n'est pas une cible d'orbite finale. Rayon/Altitude = 0.
		SimuCore::Structures::Planet("Soleil",  1.989e30, 696340, 0, -696340),

		// Planètes avec Rayon Physique et Altitude Cible de 200 km (200e3 m)
		SimuCore::Structures::Planet("Mercure", 3.3011e23, 2439.7,    0,     -75, glm::dvec3(  57.91e9, 0, 0), glm::dvec3(0, 47.87e3, 0)),
		SimuCore::Structures::Planet("Venus",   4.8675e24, 6051.8,  175,   10921, glm::dvec3( 108.21e9, 0, 0), glm::dvec3(0, 35.02e3, 0)),
		SimuCore::Structures::Planet("Terre",   5.9722e24, 6371.0,  450,   19611, glm::dvec3(  149.6e9, 0, 0), glm::dvec3(0, 29.78e3, 0)),
		SimuCore::Structures::Planet("Mars",    6.4171e23, 3389.5,  200,    9556, glm::dvec3( 227.92e9, 0, 0), glm::dvec3(0, 24.13e3, 0)),
		SimuCore::Structures::Planet("Jupiter", 1.8982e27,  69911, 2000, 2332435, glm::dvec3( 778.57e9, 0, 0), glm::dvec3(0, 13.07e3, 0)),
		SimuCore::Structures::Planet("Saturne", 5.6834e26,  58232, 3500, 2354450, glm::dvec3(1433.53e9, 0, 0), glm::dvec3(0,  9.68e3, 0)),
		SimuCore::Structures::Planet("Uranus",  8.6810e25,  25362, 6000, 1875176, glm::dvec3(2872.46e9, 0, 0), glm::dvec3(0,  6.80e3, 0)),
		SimuCore::Structures::Planet("Neptune", 1.0241e26,  24622, 3000, 3197294, glm::dvec3(4495.06e9, 0, 0), glm::dvec3(0,  5.43e3, 0))
	};


	AdaptedSystem::AdaptedSystem() :
		m_startPlanet(m_planets[3]),
		m_finalPlanet(m_planets[4]),
		m_sun(m_planets[0]),
		m_rocket(Rocket(daysInSeconds(365.25), std::vector<std::pair<SimuCore::Structures::Impulsion, double>>(), 200e3, 4.5)),
		m_time(0),
		m_startAngle(0),
		m_finalAngle(0),
		m_MaxTime(daysInSeconds(365.25)),
		m_deltaTime(1000)
	{
	}

	AdaptedSystem::AdaptedSystem(PlanetsName start_planet, PlanetsName final_planet, double start_angle, double final_angle, Rocket rocket, double max_duration, double dt_seconds) :
		m_startPlanet(m_planets[3]),
		m_finalPlanet(m_planets[4]),
		m_sun(m_planets[0]),
		m_rocket(rocket),
		m_time(0),
		m_startAngle(start_angle),
		m_finalAngle(final_angle),
		m_MaxTime(daysInSeconds(max_duration)),
		m_deltaTime(dt_seconds)
	{
		assert(start_planet != final_planet && "Les planètes de départ et d'arrivée doivent être différentes.");

		InitPlanet(true, start_planet);
		InitPlanet(false, final_planet);

		SetAnglePlanet(m_startPlanet, start_angle);
		SetAnglePlanet(m_finalPlanet, final_angle);
	}


	AdaptedSystem::RocketState AdaptedSystem::Run(std::function<RocketState(const Rocket&) > state) {
		std::vector<Entity*> entities; // vérifier la mémoire durant l'éxécution
		entities.reserve(4);
		entities.push_back(&m_rocket);
		entities.push_back(&m_sun);
		entities.push_back(&m_startPlanet);
		entities.push_back(&m_finalPlanet);


		// on simule jusqu'à m_MaxTime ou la mort de la fusée ou son succès
		const size_t MAX_ITERATIONS = static_cast<const size_t>(m_MaxTime / m_deltaTime);
		RocketState current_state = RocketState::NEUTRAL;
		size_t iteration = 0;

		while (current_state == RocketState::NEUTRAL && iteration < MAX_ITERATIONS) {
			Integrator::IntegrateStep(entities, m_deltaTime, m_time); // ici à 15h31 jeudi 27/11/2025
			m_rocket = *dynamic_cast<Rocket*>(entities[0]);

			current_state = state(m_rocket);
			m_time += m_deltaTime;
			iteration++;
		}

		return current_state;
	}

	void AdaptedSystem::Reset(double t) {
		m_time = 0;

		m_sun.position = glm::dvec3(0, 0, 0);
		m_sun.velocity = glm::dvec3(0, 0, 0);
		m_sun.forces = glm::dvec3(0, 0, 0);

		// on calcule les angles des planètes à l'instant t
		// Pour cela, on calcule la vitesse angulaire de chaque planète, et on multiplie par t pour obtenir l'angle
		{
			m_startAngle = m_startPlanet.getAngularVelocity(m_sun) * t;
			m_finalAngle = m_finalPlanet.getAngularVelocity(m_sun) * t;

			SetAnglePlanet(m_startPlanet, m_startAngle);
			SetAnglePlanet(m_finalPlanet, m_finalAngle);

			SetPlanetSpeed(m_startPlanet, m_startAngle);
			SetPlanetSpeed(m_finalPlanet, m_finalAngle);
		}

		m_startPlanet.forces = glm::dvec3(0, 0, 0);
		m_finalPlanet.forces = glm::dvec3(0, 0, 0);

		// besoin de reset la rocket ???  TODO
	}

	AdaptedSystem::Real AdaptedSystem::Score(const std::vector<std::vector<Real>>& individu) {
		auto [rocket, gen_state] = IndividualToRocket(individu, *this);
		m_rocket = rocket;

		if (gen_state != GenerationState::VALID) {
			return m_LowScore;
		}



		// on lance une simulation et on récupère :
		//	• la position finale
		//	• si l'accélération a été trop forte (--> c++ 23 pour std::expected) --> pas besoin, gestion de l'état de la fusée par une énumération RocketState
		//	• si le temps est trop long (--> c++ 23 pour std::expected) --> pas besoin, on arrête la simulation avant + gestion de l'état de la fusée par une énumération RocketState
		//	• la position finale de l'objectif
		//	• le coût énergétique
		//	• le temps de vol
		//



		auto rocket_state = [this](const Rocket& rocket) -> RocketState {
			// on check l'acceleration
			constexpr Real acceleration_maximale = constants::g * 5; // 5g - m/s² - Accélération maximale tolérée pour les humains dans la fusée
			Real acceleration = rocket.acceleration;
			if (acceleration > acceleration_maximale) {
				return RocketState::DEAD_ACCELERATION_TOO_HIGH;
			}

			// check les collisions --> brute force car seulement 3 comparaisons
			bool collided = rocket_collide_with(ObjectName::SUN) || rocket_collide_with(ObjectName::START) || rocket_collide_with(ObjectName::FINAL);

			if (collided) {
				// déterminer la vitesse au moment de la collision
				Real speed = glm::length(rocket.velocity);

				// Il faut déterminer si la vitesse est "faible" ou "élevée"
				// On considère qu'une vitesse inférieure à la vitesse de libération de l'astre qu'on a touché est une "faible" vitesse
				Real vitesse_de_liberation;
				if (rocket_collide_with(ObjectName::SUN)) {
					vitesse_de_liberation = m_sun.extractionVelocity(m_sun.getRadius());
				}
				else if (rocket_collide_with(ObjectName::START)) {
					vitesse_de_liberation = m_startPlanet.extractionVelocity(m_startPlanet.getRadius());
				}
				else { // FINAL
					vitesse_de_liberation = m_finalPlanet.extractionVelocity(m_finalPlanet.getRadius());
				}

				if (speed < vitesse_de_liberation) {
					return RocketState::DEAD_TOUCH_PLANET_LOW_SPEED;
				}
				else {
					return RocketState::DEAD_TOUCH_PLANET_HIGH_SPEED;
				}
			}

			// on check la position finale, et si la fusée est en orbite autour de la planète finale
			glm::dvec3 position_dans_referentiel_planete = rocket.position - m_finalPlanet.position;
			Real distance = glm::length(position_dans_referentiel_planete);
			bool distance_ok = (m_finalPlanet.minOrbitRadius() <= distance && distance <= m_finalPlanet.maxOrbitRadius());

			if (distance_ok) {
				glm::dvec3 vecteur_vitesse_referentiel_planete_finale = rocket.velocity - m_finalPlanet.velocity;
				Real norme_vitesse = glm::length(vecteur_vitesse_referentiel_planete_finale);

				// TODO check les formules
				double energie_orbitale = 0.5 * norme_vitesse * norme_vitesse - m_finalPlanet.getMu() / distance;

				if (energie_orbitale >= 0 - constants::epsilon) return RocketState::NEUTRAL; // vérifier qu'on enleve pas trop de trajectoires

				
				glm::dvec3 presque_moment_cinetique = glm::cross(position_dans_referentiel_planete, vecteur_vitesse_referentiel_planete_finale);
				double constante_des_aires = glm::length(presque_moment_cinetique);

				double demi_grand_axe = -0.5 * m_finalPlanet.getMu() / energie_orbitale;
				double excentricite = std::sqrt(1 + (2 * energie_orbitale * (constante_des_aires/ m_finalPlanet.getMu()) * (constante_des_aires / m_finalPlanet.getMu())));

				double perige = demi_grand_axe * (1 - excentricite);
				double apogee = demi_grand_axe * (1 + excentricite);

				// on reste dans l'anneau si le perigé est plus grand que l'orbite min et que l'apogée est plus petit que l'orbite min
				// car si on dépasse, le soleil a une forte influences
				bool vitesse_ok = (perige >= m_finalPlanet.minOrbitRadius()) && (apogee <= m_finalPlanet.maxOrbitRadius());

				return vitesse_ok ? RocketState::VALID : RocketState::NEUTRAL;
			}

			else {
				return RocketState::NEUTRAL;
			}

		};

		RocketState final_state = Run(rocket_state);
		switch (final_state)
		{

		case RocketState::DEAD_TOUCH_PLANET_HIGH_SPEED:
			// on veut que la borne inf soit 1/epsilon
			// on veut que la borne sup soit 2/epsilon
			return (1 / (constants::epsilon + glm::length(m_rocket.velocity)) ) + 1/constants::epsilon; // TODO : ajuster la formule
			break;

		case RocketState::DEAD_ACCELERATION_TOO_HIGH:
			// on veut que la borne inf soit 2/epsilon
			// on veut que la borne sup soit 3/epsilon
			return (1 / (constants::epsilon + m_rocket.acceleration)) + 2 / constants::epsilon; // TODO : ajuster la formule
			break;

		case RocketState::DEAD_TOUCH_PLANET_LOW_SPEED:
			// on veut que la borne inf soit 3/epsilon
			// on veut que la borne sup soit 4/epsilon
			return (1 / (constants::epsilon + glm::length(m_rocket.velocity))) + 3 / constants::epsilon; // TODO : ajuster la formule
			break;

		case RocketState::NEUTRAL:
			return HandleScoreNeutralState();
			break;

		case RocketState::VALID:
			return HandleScoreValidState();
			break;

		default:
			throw std::runtime_error("Unknown RocketState encountered in Score calculation.");
			break;
		}	
	}

	AdaptedSystem::Real AdaptedSystem::RingSize_meter() const {
		return m_startPlanet.maxOrbitRadius() - m_startPlanet.minOrbitRadius();
	}

	void AdaptedSystem::InitPlanet(bool is_start_planet, PlanetsName name) {
		switch (name) {
		case PlanetsName::Mercure:
			if (is_start_planet)	m_startPlanet = m_planets[1];
			else					m_finalPlanet = m_planets[1];
			break;

		case PlanetsName::Venus:
			if (is_start_planet)	m_startPlanet = m_planets[2];
			else					m_finalPlanet = m_planets[2];
			break;

		case PlanetsName::Terre:
			if (is_start_planet)	m_startPlanet = m_planets[3];
			else					m_finalPlanet = m_planets[3];
			break;

		case PlanetsName::Mars:
			if (is_start_planet)	m_startPlanet = m_planets[4];
			else					m_finalPlanet = m_planets[4];
			break;

		case PlanetsName::Jupiter:
			if (is_start_planet)	m_startPlanet = m_planets[5];
			else					m_finalPlanet = m_planets[5];
			break;

		case PlanetsName::Saturne:
			if (is_start_planet)	m_startPlanet = m_planets[6];
			else					m_finalPlanet = m_planets[6];
			break;

		case PlanetsName::Uranus:
			if (is_start_planet)	m_startPlanet = m_planets[7];
			else					m_finalPlanet = m_planets[7];
			break;

		case PlanetsName::Neptune:
			if (is_start_planet)	m_startPlanet = m_planets[8];
			else					m_finalPlanet = m_planets[8];
			break;
		}
	}

	void AdaptedSystem::SetAnglePlanet(Planet& planet, double theta) {
		double norme = glm::length(planet.position);

		double abscisse = std::cos(theta) * norme;
		double ordonnee = std::sin(theta) * norme;

		planet.position = glm::dvec3(abscisse, ordonnee, 0);
	}

	void AdaptedSystem::SetPlanetSpeed(Planet& planet, double theta) {
		double norme = glm::length(planet.velocity);

		double abscisse = std::cos(theta + constants::PI / 2) * norme;
		double ordonnee = std::sin(theta + constants::PI / 2) * norme;

		planet.velocity = glm::dvec3(abscisse, ordonnee, 0);
	}


	AdaptedSystem::Real AdaptedSystem::HandleScoreValidState() const
	{
		// On sait que la fusée est en orbite stable autour de la planète cible :
		// On peut donc renvoyer un score parfait ou presque (en fonction du coût énergétique et du temps)
		// Il reste à déterminer ce qu'est un score parfait. Il nous faut un meilleur score que dans le cas neutre. C'est à dire un majorant du score neutre.
		
		// Un majorant simple est de prendre le cas où la position et la vitesse sont parfaites (donc l'influence position et vitesse est maximale).
		// Comme on ajoute un epsilon pour éviter la division par zéro, on peut prendre 1 / epsilon comme influence maximale, pour chaque influence.
		constexpr Real Majorant_etat_neutre = 5 / constants::epsilon;


		Real cout_energetique = m_rocket.getDeltaM(); // TODO
		Real tof = m_time;
		return 1/(tof*cout_energetique + constants::epsilon) + Majorant_etat_neutre;
	}

	AdaptedSystem::Real AdaptedSystem::HandleScoreNeutralState() const
	{
		Real cout_energetique = m_rocket.getDeltaM(); // TODO 
		Real tof = m_time;

		Real influence_temps_energie = 0.5/(cout_energetique * tof + constants::epsilon);



		Real influence_position = 0.5 / constants::epsilon; // Trouver un maximum pertinent

		glm::dvec3 projection_on_ring = m_finalPlanet.targetRadius() * glm::normalize(m_rocket.position - m_finalPlanet.position);
		Real distance_to_ring = glm::length(m_rocket.position - projection_on_ring);

		
		bool position_ok = (m_finalPlanet.minOrbitRadius() <= distance_to_ring && distance_to_ring <= m_finalPlanet.maxOrbitRadius());
		if (!position_ok)
		{
			influence_position = 0.5 / (distance_to_ring + constants::epsilon); // On veut que l'influence diminue quand on s'éloigne de la distance cible
		}
		// on veut que la borne inf soit 4/epsilon, donc on ajoute 4 / epsilon
		// donc la borne sup est 0.5/epsilon + 0.5/epsilon + 4/epsilon = 5/epsilon




		// TODO : vérifier avec des assertions que les influences sont bien dans les bornes attendues.
		return influence_temps_energie + influence_position + 4/constants::epsilon;
	}

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
			return m_LowScore * 10; // Pour pouvoir différencier des autres états invalides
			break;
		case SimuCore::GenerationState::INVALID_GENES_OUT_OF_BOUNDS_VELOCITY:
			return m_LowScore * 100; // Pour pouvoir différencier des autres états invalides
			break;
		case SimuCore::GenerationState::INVALID_GENES_OUT_OF_BOUNDS_IMPULSION:
			return m_LowScore * 1000; // Pour pouvoir différencier des autres états invalides
			break;
		default:
			break;
		}
	}

	bool AdaptedSystem::rocket_collide_with(ObjectName name) const {
		glm::dvec3 object_pos;
		double object_radius;

		switch (name) {
		case ObjectName::SUN :
			object_pos = m_sun.position;
			object_radius = m_sun.getRadius();
			break;

		case ObjectName::START :
			object_pos = m_startPlanet.position;
			object_radius = m_startPlanet.getRadius();
			break;

		case ObjectName::FINAL :
			object_pos = m_finalPlanet.position;
			object_radius = m_finalPlanet.getRadius();
			break;

		default:
			throw std::runtime_error("Unknown ObjectName encountered in rocket_collide_with.");
			break;
		}

		double distance = glm::length(m_rocket.position - object_pos);

		return distance < object_radius + constants::epsilon;
	}
};
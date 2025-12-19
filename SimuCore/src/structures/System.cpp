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

	void AdaptedSystem::Reset() {
		m_time = 0;

		SetAnglePlanet(m_startPlanet, m_startAngle);
		SetAnglePlanet(m_finalPlanet, m_finalAngle);

		SetPlanetSpeed(m_startPlanet, m_startAngle);
		SetPlanetSpeed(m_finalPlanet, m_finalAngle);

		m_sun.position = glm::dvec3(0, 0, 0);
		m_sun.velocity = glm::dvec3(0, 0, 0);

		m_sun.forces = glm::dvec3(0, 0, 0);
		m_startPlanet.forces = glm::dvec3(0, 0, 0);
		m_finalPlanet.forces = glm::dvec3(0, 0, 0);

		// besoin de reset la rocket ???  TODO
	}

	AdaptedSystem::Real AdaptedSystem::Score(const std::vector<std::vector<Real>>& individu) {
		// On initialise la fusée avec les positions et vitesses initiales.
		const double max_reels_position = RingSize_meter();

		auto convert_to = [max_reels_position](double x, double min, double max) {
			return convertIntervals(-max_reels_position, max_reels_position, -max, max, x);
			};

		//////// gestion de la position initiale qui va conditionner la transformation des vecteurs en impulsions
		// Pour mettre la fusée sur un anneau, la position initiale donnée doit être ajoutée avec un vecteur colinéaire à lui même, de norme
		//		égale au rayon minimal de l'anneau.
		glm::dvec3 position_initiale = glm::dvec3(
			individu[0][0],
			individu[0][1],
			//individu[0][2]
			0
		);

		// si la distance de position_initiale à l'origine n'est pas dans le cerlce de rayon <jsp quoi> alors t'es mort
		if (glm::length(position_initiale) > RingSize_meter()) {
			return m_LowScore; // score pas fou
		}

		// transformation 1er vecteur en position initiale
		{
			// on récupère un vecteur colinéaire au vecteur position_initiale
			glm::dvec3 direction = glm::normalize(position_initiale);
			position_initiale += m_startPlanet.minOrbitRadius() * direction; // anneau centré sur le soleil (en considérant chaque position_initiale possible)
		}



		//////// gestion des impulsions

		//// gestion de la 1ère impulsion (avec son temps t1)
		// sa norme doit être plus grande que la vitesse de la libération et sa direction doit être dans une sphère
		// On réinitialise le temps à t1 ("t1 devient 0").

		// Calculons le nouvel état initial, il suffit de reset les 2 planètes
		// pour cela, on met le temps t1 entre 0 et <le temps qu'il faut aux deux planètes pour revenir à la même configuration (on reste en 2D pour le moment)>

		double temps_1ere_impulsion = individu[2][0];
		{
			double dist_init_planet_to_sun = glm::length(m_startPlanet.position);
			double dist_final_planet_to_sun = glm::length(m_finalPlanet.position);
			double omega_planete_initiale = std::sqrt(m_planets[0].getMu() / (dist_init_planet_to_sun * dist_init_planet_to_sun * dist_init_planet_to_sun));
			double omega_planete_finale = std::sqrt(m_planets[0].getMu() / (dist_final_planet_to_sun * dist_final_planet_to_sun * dist_final_planet_to_sun));
			double temps_1ere_impulsion_max = std::abs(2 * constants::PI / (omega_planete_initiale - omega_planete_finale));

			// on veut remplacer l'intervalle [min_real, max_real] en [0, temps_1ere_impulsion_max]
			/*// Or min_real = - RingSize_meter() = - max_real
			temps_1ere_impulsion += RingSize_meter();			// dans [0, min_real + max_real]
			temps_1ere_impulsion /= 2 * RingSize_meter();		// dans [0, 1]
			temps_1ere_impulsion *= temps_1ere_impulsion_max;	// on a gagné*/

			temps_1ere_impulsion = convert_to(temps_1ere_impulsion, 0, temps_1ere_impulsion_max);


			// modifier les position et vitesses initiales des deux planètes
			// calculons les angles
			m_startAngle = omega_planete_initiale * temps_1ere_impulsion;
			m_startAngle = omega_planete_finale * temps_1ere_impulsion;
			
			SetAnglePlanet(m_startPlanet, m_startAngle);
			SetAnglePlanet(m_finalPlanet, m_finalAngle);

			SetPlanetSpeed(m_startPlanet, m_startAngle);
			SetPlanetSpeed(m_finalPlanet, m_finalAngle);
		}

		// on coupe les vecteur trop grand (pas dans la sphère)
		glm::dvec3 vitesse_initiale = glm::dvec3(
			individu[1][0],
			individu[1][1],
			//individu[1][2]
			0
		);

		if (glm::length(vitesse_initiale) > max_reels_position) {
			return m_LowScore; // score pas fou
		}

		{
			glm::dvec3 direction = glm::normalize(vitesse_initiale);
			vitesse_initiale += m_startPlanet.extractionVelocity( glm::length(position_initiale) ) * direction;
		}


		//// gestion des autres impulsions
		// Les individus sont composés de vecteurs générés aléatoirement dans un pavé droit, où chaque composante est comprise entre [min_real, max_real]
		// Or on souhaite que les impulsions de la fusée soient de "direction aléatoire uniforme" (sphère != pavé droit)

		size_t nombre_impulsions = (individu.size() - 1) / 2;
		std::vector<std::pair<Structures::Impulsion, double>> toutes_les_impulsions;
		toutes_les_impulsions.reserve(nombre_impulsions);
		toutes_les_impulsions.emplace_back(Structures::Impulsion(vitesse_initiale), temps_1ere_impulsion);

		{
			double somme_des_temps = 0;

			for (size_t i = 1; i < nombre_impulsions; i++) {
				// individu[2 * i + 2][0] correspond à l'écart de temps entre l'impulsion i et i-1
				// on remet instant entre 0 et <jsp encore> : TODO (on peut aussi forcer la somme des ecart de temps à être plus petit que m_MaxTime)
				// Afin de faciliter la recherche des solutions, on met les temps entre 0 et m_Max_Time.
				double instant = convert_to(individu[2 * i + 2][0], 0, m_MaxTime);

				glm::dvec3 impulsion = glm::dvec3{
					individu[2 * i + 1][0],
					individu[2 * i + 1][1],
					//convert_to(individu[2 * i + 1][2], -max_reels_position, max_reels_position)
					0
				};

				const auto& [impul_precedente, temps_precedent] = toutes_les_impulsions[i - 1];

				somme_des_temps += instant;
				instant += temps_precedent;

				// on s'assure que l'impulsion impulsion est généré dans la sphère :
				if (glm::length(impulsion) > RingSize_meter()) {
					return m_LowScore;
				}

				toutes_les_impulsions.emplace_back(Structures::Impulsion(impulsion), instant);
			}

			if (somme_des_temps > m_MaxTime) {
				return m_LowScore;
			}
		}


		m_rocket.position = position_initiale + m_startPlanet.position; // anneau centré sur la planète initiale
		m_rocket.velocity = m_startPlanet.velocity; // impulsion pour moduler ce terme.
		m_rocket.acceleration = 0; // une acceleration valide pour que le RocketState initial soit sûr d'être valide. (Pas utiliser dans les calculs)
		m_rocket.setImpulsions(toutes_les_impulsions);
















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
				return RocketState::DEAD;
			}

			// check les collisions --> brute force car seulement 3 comparaisons
			bool collided = rocket_collide_with(ObjectName::SUN) || rocket_collide_with(ObjectName::START) || rocket_collide_with(ObjectName::FINAL);

			if (collided) {
				return RocketState::DEAD;
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
		case RocketState::DEAD:
			return m_LowScore;
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
		constexpr Real Majorant_etat_neutre = 2 / constants::epsilon;


		Real cout_energetique = m_rocket.getDeltaM(); // TODO
		Real tof = m_time;
		return 1/(tof*cout_energetique + constants::epsilon) + Majorant_etat_neutre;
	}

	AdaptedSystem::Real AdaptedSystem::HandleScoreNeutralState() const
	{
		Real cout_energetique = m_rocket.getDeltaM(); // TODO 
		Real tof = m_time;

		Real influence_temps_energie = 1/(cout_energetique * tof + constants::epsilon);



		Real influence_position = 1 / constants::epsilon; // Trouver un maximum pertinent

		glm::dvec3 projection_on_ring = m_finalPlanet.targetRadius() * glm::normalize(m_rocket.position - m_finalPlanet.position);
		Real distance_to_ring = glm::length(m_rocket.position - projection_on_ring);

		
		bool position_ok = (m_finalPlanet.minOrbitRadius() <= distance_to_ring && distance_to_ring <= m_finalPlanet.maxOrbitRadius());
		if (!position_ok)
		{
			influence_position = 1 / (distance_to_ring + constants::epsilon); // On veut que l'influence diminue quand on s'éloigne de la distance cible
		}


		// TODO : vérifier avec des assertions que les influences sont bien dans les bornes attendues.
		return influence_temps_energie + influence_position;
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
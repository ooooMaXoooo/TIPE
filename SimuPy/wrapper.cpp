#include "pch.h"
#include <SimuCore/theory/formula.h>
#include <SimuCore/constants.h>
#include <SimuCore/integrator/integrator.h>
#include <SimuCore.h>

#pragma comment(lib, "SimuCore.lib") 

namespace py = pybind11;

#ifdef _DEBUG
#define MODULE_NAME TIPE_SimuOrbit_d
#else
#define MODULE_NAME TIPE_SimuOrbit
#endif

PYBIND11_MODULE(MODULE_NAME, m) {
    // Ajout du sous-module 'stumpff'
    pybind11::module_ stumpff_mod = m.def_submodule("stumpff", "Fonctions de Stumpff");
    stumpff_mod.def("C", &stumpff::C, R"pbdoc(
        Fonction Stumpff C(z)
        Calcule C(z) selon la valeur de z (positive, négative ou proche de zéro)
        )pbdoc"
    );

    stumpff_mod.def("S", &stumpff::S, R"pbdoc(
        Fonction Stumpff S(z)
        Calcule S(z) avec traitement des cas limites
        )pbdoc"
    );







    py::module_ orbit_mod = m.def_submodule("orbit", "Fonctions orbitales");
    orbit_mod.def("lambert_universal", &orbit::lambert_universal,
        py::arg("r1"), py::arg("r2"), py::arg("tof"), py::arg("mu") = SimuCore::constants::mu,
        "Résout le problème de Lambert universel");

    orbit_mod.def("lambert_batch", &orbit::lambert_batch,
        py::arg("r1"), py::arg("r2_list"), py::arg("tof_list"), py::arg("mu") = SimuCore::constants::mu,
        "Calcule DeltaV (||v1||+||v2||) pour une liste de (r2, tof) en parallèle (OpenMP)");

    orbit_mod.def("lambert_batch_numpy",
        [](const std::array<double, 3>& r1,
            const std::vector<std::array<double, 3>>& r2_list,
            const std::vector<double>& tof_list,
            double mu) -> py::array_t<double>
        {
            size_t N = r2_list.size();
            if (N != tof_list.size()) {
                throw std::invalid_argument("r2_list et tof_list doivent avoir la même taille");
            }

            // Crée directement un numpy.ndarray de taille N
            py::array_t<double> result(N);
            auto r = result.mutable_unchecked<1>();  // accès direct (1D)

            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(N); i++) {
                try {
                    auto res = orbit::lambert_universal(r1, r2_list[i], tof_list[i], mu);
                    const auto& v1 = res.first;
                    const auto& v2 = res.second;

                    if (std::isfinite(v1[0]) && std::isfinite(v2[0])) {

                        auto norme = [](const std::array<double, 3>& u) -> double {
                            return std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
                            };

                        auto substract = [](const std::array<double, 3>& u, const std::array<double, 3>& v) {
                            return std::array<double, 3>{u[0] - v[0], u[1] - v[1], u[2] - v[2] };
                            };

                        auto calcul_vecteur_cercle = [norme, mu](const std::array<double, 3>& r1) {
                            std::array<double, 3> v = { -r1[1], r1[0], r1[2] };
                            auto norme_voulu = std::sqrt(mu / norme(r1)) / norme(r1);
                            for (int i = 0; i < 3; i++) { v[i] *= norme_voulu; }
                            return v;
                            };

                        std::array<double, 3> vecteur_vitesse_orbite_depart = calcul_vecteur_cercle(r1);
                        std::array<double, 3> vecteur_vitesse_orbite_fin = calcul_vecteur_cercle(r2_list[i]);


                        r(i) = norme(substract(v1, vecteur_vitesse_orbite_depart)) + norme(substract(vecteur_vitesse_orbite_fin, v2));
                    }
                    else {
                        r(i) = NAN;
                    }
                }
                catch (...) {
                    r(i) = NAN;
                }
            }
            return result;
        },
        py::arg("r1"), py::arg("r2_list"), py::arg("tof_list"), py::arg("mu") = SimuCore::constants::mu,
        "Version optimisée: retourne directement un numpy.ndarray (sans copie)"
    );


    


	py::module_ genetic_mod = m.def_submodule("genetic", "Modules donnant accès à l'algorithme génétique");


    genetic_mod.def(
		"GetIndividualTrajectory__configD32_2_2", // D pour double, 32 pour uint32_t, 2 pour impulsions, 2 la dimension de l'espace
        [](double dt_seconds,
            double max_time_days,
			double simulation_duration_days,
            double mass_kg,
            uint8_t startPlanet,
            uint8_t finalPlanet,
            std::array<std::array<uint32_t, 2>, 5> genome // 5 chromosomes, 2 gènes par chromosome
            ) -> std::array<std::vector<double>, 2>
        {
            double lifetime = max_time_days; // durée de simulation en jours

            SimuCore::Systems::AdaptedSystem sy{
                static_cast<SimuCore::Systems::PlanetsName>(startPlanet),   // planète de départ
                static_cast<SimuCore::Systems::PlanetsName>(finalPlanet),  // planète d'arrivée
                SimuCore::Structures::Rocket(
                    lifetime, // -> durée de vie de la fusée en jours
                    std::vector<std::pair<SimuCore::Structures::Impulsion, double>>(),
                    mass_kg,
                    2.22),                                     // -> vitesse d'éjection des gaz en km/s  (2.22 pour terre mars, 3 pour terre jupiter)
                lifetime,                                   // durée de simulation en jours
                dt_seconds};

			sy.Initialize(); // initialise les positions des planètes, etc. (doit être appelé avant de pouvoir convertir un individu en fusée)
            std::cout << "[C++] RingSize_meter=" << sy.RingSize_meter()
                << " MaxTime=" << sy.getMaxTime() << std::endl;

            using ConfigType = genetic::Config<double, uint32_t, 5, 2>;
            ConfigType config;

			// initialisation des paramètres de l'algorithme génétique (la pluplart ne sont pas important ici
            {
                config.dimension = 2;							        // on travaille dans l'espace
                config.number_of_vectors = 5;		                    // un vecteur = un gène

                config.enable_saving = false;							// désactive la sauvegarde dans les fichiers car ce n'est pas encore implémenté

                config.population_size = 4;
                config.max_generations = 10;
                config.print_interval = 15;

                config.initial_mutation_probability = 0.8;				// paramètre à ajuster
                config.initial_self_adaptation_probability = 0.8;		// paramètre à ajuster

                config.custom_mutation_proba = -1;					// on prend la proba uniforme


                {
                    // générer des réels tels que |x| < epaisseur de l'anneau
                    config.max_real = sy.RingSize_meter();
                    config.min_real = -config.max_real;
                }

                config.enable_elitism = false;
                config.enable_auto_adaptation = false;
            }
        
            genetic::Individu<ConfigType> ind{};

			ind.set_config(&config); // check, si nécéssaire, l'intervalle [min_real, max_real]
			ind.set_genome(genome);

			std::vector<std::vector<double>> vecs = ind.to_real_vectors();
            std::cout << "[C++] avant IndividualToRocket" << std::endl;
			auto [rocket, gen_status] = SimuCore::IndividualToRocket<double>(vecs, sy);

            if (gen_status != SimuCore::GenerationState::VALID) {
                // génome invalide, retourner un résultat vide ou lever une exception
                return std::array<std::vector<double>, 2>{};
            }


            std::vector<glm::dvec3> trajectory;
			const size_t NbIterations = static_cast<size_t>(max_time_days * 24 * 3600 / dt_seconds);
			trajectory.reserve(NbIterations);
            
            /*std::cout << "[C++] avant GetRocketTrajectory" << std::endl;
            sy.GetRocketTrajectory(trajectory);
            std::cout << "[C++] apres GetRocketTrajectory, taille=" << trajectory.size() << std::endl;*/

            auto position_callback = [&](const glm::dvec3& pos) {
                trajectory.push_back(pos);
                };

            
            std::cout << "[C++] m_time=" << sy.getCurrentTime()
                << " start_idx=" << sy.getStartPlanetPositionIndice()
                << " rocket_pos=(" << rocket.position.x
                << "," << rocket.position.y << ")\n";

            sy.Score(rocket, gen_status, position_callback, false, simulation_duration_days);

			// Convertir la trajectoire en format numpy
			std::array<std::vector<double>, 2> result;

			result[0].reserve(trajectory.size());
			result[1].reserve(trajectory.size());
			for (const auto& p : trajectory) {
				result[0].emplace_back(p.x);
				result[1].emplace_back(p.y);
			}
			return result;
        },
		py::arg("dt_seconds"), py::arg("max_time_days"),
        py::arg("simulation_duration_days"), py::arg("mass_kg"),
		py::arg("startPlanet"), py::arg("finalPlanet"),
		py::arg("genome"),
        "Simule la trajectoire d'une sonde spatiale à partir d'un génome donné (5 chromosomes, 2 gènes par chromosome)"
    );





	py::module_ integrator_mod = m.def_submodule("integrator", "Intégrateurs numériques");

    integrator_mod.def("Simulate",
        [](double dt_seconds, double max_time_days, double mass_kg,
            const std::array<double, 3>& initial_pos_UA,
            const std::array<double, 3>& initial_velocity_km_s,
            size_t final_planet_pos_indice,
            const std::vector<std::array<double, 3>>& impulsions,
            const std::vector<double>& dates_days,
            uint8_t startPlanet,
            uint8_t finalPlanet) -> std::array<std::vector<double>, 2>
        {
            if (impulsions.size() != dates_days.size())
                throw std::invalid_argument("nb_impulsions, impulsions et dates doivent avoir la meme taille");

            glm::dvec3 pos(initial_pos_UA[0], initial_pos_UA[1], initial_pos_UA[2]);
			glm::dvec3 vel(initial_velocity_km_s[0], initial_velocity_km_s[1], initial_velocity_km_s[2]);

            std::vector<glm::dvec3> glm_impulsions;
            glm_impulsions.reserve(impulsions.size());
            for (const auto& imp : impulsions)
                glm_impulsions.push_back(glm::dvec3(imp[0], imp[1], imp[2]));

            std::vector<glm::dvec2> positions;

            SimuCore::Integrator::Simulate(dt_seconds, max_time_days, mass_kg, pos, vel,
                final_planet_pos_indice, glm_impulsions.size(),
                glm_impulsions, dates_days,
                static_cast<SimuCore::Systems::PlanetsName>(startPlanet),
                static_cast<SimuCore::Systems::PlanetsName>(finalPlanet),
                false, &positions);

            std::array<std::vector<double>, 2> result;
            result[0].reserve(positions.size());
            result[1].reserve(positions.size());

            for (const auto& p : positions) {
                result[0].emplace_back(p.x);
                result[1].emplace_back(p.y);
            }

            return result;
        },
        py::arg("dt_seconds"), py::arg("max_time_days"), py::arg("mass_kg"),
        py::arg("initial_pos_UA"), py::arg("initial_velocity_km_s"),
        py::arg("final_planet_pos_indice"),
        py::arg("impulsions_km_s"), py::arg("dates_days"),
        py::arg("startPlanet"), py::arg("finalPlanet"),
        "Simule la trajectoire d'une sonde spatiale avec des impulsions données"
    );
}

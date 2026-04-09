#pragma once

#include <pch.h>
#include "../structures/Rocket.h"
#include "../structures/System.h"
#include "../utility.h"
#include <SimuCore\integrator\integrator.h>

#include <Galib/genetic.hpp>

#include <DataExport/AsyncDataExporter.h>
#include <DataExport/TrajectoryData.h>
#include <DataExport/GeneticData.h>
#include <DataExport/PhysicsData.h>

namespace SimuCore {
	namespace Optimization {

		void writeTrajectory(std::ofstream& filestream, const std::vector<glm::dvec3>& trajectory) {
			if (!filestream.is_open()) {
				std::cerr << "\nErreur lors de l'ouverture d'un fichier" << std::endl;
				std::abort();
			}

			for (const glm::dvec3& vec : trajectory) {
				filestream << vec.x << ';' << vec.y << '\n';
			}
		}

		inline std::string generate_snapshot_directory() {
			// Récupérer le temps actuel
			auto now = std::chrono::system_clock::now();
			std::time_t t_now = std::chrono::system_clock::to_time_t(now);
			std::tm local_tm;

#if defined(_WIN32) || defined(_WIN64)
			localtime_s(&local_tm, &t_now);  // Windows
#else
			localtime_r(&t_now, &local_tm);  // Linux / macOS
#endif

			// Construire le nom sous la forme simu_jour_mois_annee_heure_minute_seconde
			std::ostringstream oss;
			oss << "simu_"
				<< std::setfill('0') << std::setw(2) << local_tm.tm_mday << "_"
				<< std::setfill('0') << std::setw(2) << local_tm.tm_mon + 1 << "_"
				<< local_tm.tm_year + 1900 << "_"
				<< std::setfill('0') << std::setw(2) << local_tm.tm_hour << "_"
				<< std::setfill('0') << std::setw(2) << local_tm.tm_min << "_"
				<< std::setfill('0') << std::setw(2) << local_tm.tm_sec;

			return oss.str();
		}

		inline std::string generate_snapshot_filename(const std::string& extension = "txt", int id=0, const char* suffix = "") {
			std::ostringstream oss;
			oss << "gen_" << id
				<< "_" << suffix
				<< "." << extension;

			return oss.str();
		}

		/// <summary>
		/// Calcul une fusée initialisée de façon à avoir la trajectoire optimale dans le système donné.
		/// </summary>
		/// <typeparam name="NbImpulsions"></typeparam>
		/// <param name="filename"> le chemin d'accès à un fichier de sauvergarde. </param>
		/// <param name="system"> le système dans lequel évolue une fusée. </param>
		/// <returns> </returns>
		template <size_t NbImpulsions>
		SimuCore::Structures::Rocket getBestRocket(const SimuCore::Systems::AdaptedSystem system,
			genetic::CrossoverType cross_type, bool elitism, bool auto_adapt,
			size_t population_size = 100, size_t max_generation = 5000,
			size_t print_interval=100, bool verbose=false,
			size_t snapshot_interval=10, bool saving_in_file=true)
		{
			using ConfigType = genetic::Config<double, uint32_t, 2*NbImpulsions + 1, 2>;
			ConfigType config;

			config.dimension = 2;									// on travaille dans l'espace
			config.number_of_vectors = 2 * NbImpulsions + 1;		// un vecteur = un gène

			config.enable_saving = false;							// désactive la sauvegarde dans les fichiers car ce n'est pas encore implémenté

			config.population_size = population_size;				// paramètre à ajuster
			config.max_generations = max_generation;				// paramètre à ajuster
			config.print_interval  = print_interval;

			config.initial_mutation_probability = 0.8;				// paramètre à ajuster
			config.initial_self_adaptation_probability = 0.8;		// paramètre à ajuster

			config.custom_mutation_proba = -1;					// on prend la proba uniforme



			{
				// générer des réels tels que |x| < epaisseur de l'anneau
				config.max_real = system.RingSize_meter();
				config.min_real = -config.max_real;
			}

			config.enable_elitism = elitism;
			config.enable_auto_adaptation = auto_adapt;
			config.crossover_method = cross_type;

			using Real = ConfigType::real_type;
			using Integer = ConfigType::integer_type;

			// associe à un génome un score
			auto fitness = [&](const std::vector<std::vector<Real>>& vecs) -> Real {
				// vec[0] --> un vecteur à 3 dimension, équivalent à la position initiale
				// vec[1] --> un vecteur à 3 dimension, équivalent à l'impulsion initiale
				// vec[2] --> un vecteur à 3 dimension, dont on n'utilise que la première composante, équivalent à l'instant de la 1ere impulsion
				// etc...

				SimuCore::Systems::AdaptedSystem local_system = system; // on copie le système (problème de concurrence)

				//local_system.Reset();		// On réinitialise le système pour que les individus commençent tous dans les mêmes confitions initiales.
				// On n'a pas besoin de réinitialiser car la fonction Score le fait déjà.

				return local_system.Score(vecs);		// A VERIFIER !!!!! OK ?
			};

			genetic::GeneticAlgorithm<ConfigType> ga(config, fitness); // création d'un algorithme génétique


			using Ind = typename genetic::Individu<ConfigType>;
			using Pop = typename std::vector<Ind>;


			std::filesystem::path file_directory =
				std::filesystem::absolute(std::filesystem::current_path()) /
				"simulation_data" /
				generate_snapshot_directory();

			if (saving_in_file) {
				std::filesystem::create_directories(file_directory);
			}

			AsyncDataExporter generationalExporter;

			auto callback =
				[
				max_generation,
				&system,
				print_interval,
				snapshot_interval,
				saving_in_file,
				verbose,
				file_directory,
				&generationalExporter
				]
				(
					size_t gen,
					double best_fit, const genetic::Individu<ConfigType>& best_ind,
					double worst_fit, const genetic::Individu<ConfigType>& /*worst_ind*/,
					const std::vector<genetic::Individu<ConfigType>>& population
					)
				{
					if (gen % print_interval == 0 || gen == max_generation - 1) {
						///////// petites stats
						{
							size_t count_defined = 0;
							size_t count_valid = 0;
							for (auto& ind : population) {
								double fitness = ind.get_fitness();
								if (fitness > -1) { // -1 pour les erreurs de flottants
									count_valid++;
									count_defined++;
								}
								else if (fitness > SimuCore::Systems::AdaptedSystem::m_LowestScore) {
									count_defined++;
								}
								else {
									std::abort();
								}
							}

							std::cout << "\tValid individuals: " << count_valid << "/" << population.size()
								<< " (" << (100.0 * static_cast<double>(count_valid) / population.size()) << "%)\n";
							std::cout << "\tDefined individuals: " << count_defined << "/" << population.size()
								<< " (" << (100.0 * static_cast<double>(count_defined) / population.size()) << "%)\n";
						}

						///////// Distance au milieu de l'anneau
						{
							constexpr double cste = 1.0 / SimuCore::Systems::AdaptedSystem::m_CstScore;
							if (
								best_fit < 8.75 * cste
								&&
								best_fit > 8 * cste
								)
							{

								double distance_to_final_planet_best = (best_fit - 8 * cste);
								distance_to_final_planet_best /= 0.75;
								distance_to_final_planet_best = 1 / distance_to_final_planet_best;
								distance_to_final_planet_best -= SimuCore::Systems::AdaptedSystem::m_CstScore;
								distance_to_final_planet_best = distance_to_final_planet_best * distance_to_final_planet_best * distance_to_final_planet_best;


								std::cout << "\tBest distance to target: "
									<< AU_to_kilometers(distance_to_final_planet_best) << " (km) = "
									<< distance_to_final_planet_best << " (AU)\n";
							} 
							else if (best_fit < 9 * cste && best_fit > 8.75 * cste) {
								std::cout << "\tBest distance to target: " << 0 << '\n';
								
								long double mecanic_energy = best_fit - 8.75 * cste; // km/s
								mecanic_energy *= 4;
								mecanic_energy = 1 / mecanic_energy;
								mecanic_energy -= SimuCore::Systems::AdaptedSystem::m_CstScore;
								mecanic_energy *= 1e5;
								mecanic_energy = mecanic_energy * mecanic_energy;

								std::cout << "\tBest energy at target: " << mecanic_energy << " (J)\n";
							}
							else if (
								best_fit < 10 * cste
								&&
								best_fit > 9 * cste
								) {
								std::cout << "\tBest distance to target: " << 0 << '\n';
							}
						}

						///////// Type de trajectoire
						std::cout << "\tKind of trajectory (best) : " << system.TypeOfTrajectory(best_fit) << '\n';

						///////// Info physique du meilleur individu
						{
							SimuCore::Systems::AdaptedSystem copy_system = system;
							auto [rocket, state] = SimuCore::IndividualToRocket(best_ind.to_real_vectors(), copy_system);


							auto print_vec = [](glm::dvec3 vec) -> void {
								std::cout << '(' << vec.x << ", " << vec.y << ", " << vec.z << ')';
								};

							///////// Position
							{
								std::cout << "\tBest Rocket physic information :\n";
								std::cout << "\t\tPosition : "; print_vec(rocket.position); std::cout << " (unit : AU)\n";
							}

							///////// Impulsions
							{
								const std::vector<std::pair<SimuCore::Structures::Impulsion, double>>& impuls = rocket.getImpulsions();

								for (int i = 0; i < impuls.size(); i++) {
									auto& [delta_v, date] = impuls[i]; // delta_v en km/s et date en jours

									std::cout << "\t\tImpulsion " << i + 1<< " : ";
									print_vec(delta_v.GetDeltaV_vec());
									std::cout << " (unit : km/s)\n";
									std::cout << "\t\tDate " << i + 1 << " : " << date << " (jours)\n";

								}
							}

							///////// Position planete depart
							{
								std::cout << "\tStart planet :\n";
								std::cout << "\t\tPosition : ";
								print_vec(copy_system.getStartPlanetPositions()[
									copy_system.getStartPlanetStartIndice()
								]);
								std::cout << " (unit : AU)\n";
							}

							///////// Position planete arrive
							{
								std::cout << "\tFinal planet :\n";
								std::cout << "\t\tPosition : ";
								print_vec(copy_system.getFinalPlanetPositions()[
									copy_system.getFinalPlanetStartIndice()
								]);
								std::cout << " (unit : AU)\n";
							}
						}						
					}

					if (saving_in_file) {
						if (gen == 0) {
							std::filesystem::path filepath_start = file_directory / "start_planet.txt";
							std::filesystem::path filepath_final = file_directory / "final_planet.txt";

							// on ecrit les données des trajectoires des planètes dans un fichier
							//std::string filename_start = filepath_start.string();
							//std::string filename_final = filepath_final.string();

							//std::ofstream file_start(filename_start);
							//std::ofstream file_final(filename_final);

							AsyncDataExporter planetsDataExporter;
							TrajectoryData traj_data_start(system.getStartPlanetPositions(), 2);
							TrajectoryData traj_data_final(system.getFinalPlanetPositions(), 2);

							try
							{
								planetsDataExporter.enqueue(std::make_shared<TrajectoryData>(traj_data_start), filepath_start);
								planetsDataExporter.enqueue(std::make_shared<TrajectoryData>(traj_data_final), filepath_final);
							}
							catch (const std::exception& e)
							{
								std::cerr << "\n\n" << e.what() << std::endl;
								exit(EXIT_FAILURE);
							}


							//writeTrajectory(file_start, system.getStartPlanetPositions());
							//writeTrajectory(file_final, system.getFinalPlanetPositions());

							//file_start.close();
							//file_final.close();
						}

						if (gen % snapshot_interval == 0 || gen == max_generation - 1) {
							SimuCore::Systems::AdaptedSystem copy_system = system;
							///////// Envoi des données
							{
								std::filesystem::path filepath =
									file_directory /
									generate_snapshot_filename("txt", gen + 1, "rocket");

								auto [rocket, state] = SimuCore::IndividualToRocket(best_ind.to_real_vectors(), copy_system);
								copy_system.SetRocket(rocket);
								std::vector<glm::dvec3> trajectory;
								copy_system.GetRocketTrajectory(trajectory);

								TrajectoryData trajectoryData{ trajectory, 2 };
								try {
									generationalExporter.enqueue(std::make_shared<TrajectoryData>(trajectoryData), filepath);
								}
								catch (const std::exception& e) {
									std::cerr << "\n\nTrajectory writing error ||\t" << e.what() << std::endl;
									exit(EXIT_FAILURE);
								}



								// envoi des données physiques

								double tof = copy_system.getMaxTime();	// en jours
								double dt = copy_system.getDeltaTime(); // en secondes

								uint8_t startPlanetIndex = SimuCore::Systems::AdaptedSystem::GetStartPlanetID();
								uint8_t finalPlanetIndex = SimuCore::Systems::AdaptedSystem::GetFinalPlanetID();

								glm::dvec3 startPlanet_initialPosition		= copy_system.GetStartPlanet_StartPosition();
								glm::dvec3 startPlanet_finalPosition		= copy_system.GetStartPlanet_CurrentPosition();
								glm::dvec3 finalPlanet_initialPosition		= copy_system.GetFinalPlanet_StartPosition();
								glm::dvec3 finalPlanet_finalPosition		= copy_system.GetFinalPlanet_CurrentPosition();

								glm::dvec3 finalVelocity_rocket = copy_system.GetRocketVelocity(); // en km/s


								std::vector<std::pair<SimuCore::Structures::Impulsion, double>> impulsions = rocket.getImpulsions();

								std::vector<double> impulsions_times;
								std::vector<glm::dvec3> impulsions_vectors;

								impulsions_times.reserve(impulsions.size());
								impulsions_vectors.reserve(impulsions.size());

								for (auto& impuls : impulsions) {
									auto& [vec, instant] = impuls;

									impulsions_times.emplace_back(instant);
									impulsions_vectors.emplace_back(vec.GetDeltaV_vec());
								}

								constexpr uint8_t dimension = 2;

								PhysicsData phys_data{ tof, dt, startPlanetIndex, finalPlanetIndex,
									startPlanet_initialPosition,
									finalPlanet_initialPosition,
									startPlanet_finalPosition,
									finalPlanet_finalPosition,
									finalVelocity_rocket,
									impulsions_times,
									impulsions_vectors,
									dimension
								};

								std::filesystem::path filepath_physics_data =
									file_directory /
									generate_snapshot_filename("txt", gen + 1, "physics");

								try {
									generationalExporter.enqueue(std::make_shared<PhysicsData>(phys_data), filepath_physics_data);
								}
								catch (const std::exception& e) {
									std::cerr << "\n\nTrajectory writing error ||\t" << e.what() << std::endl;
									exit(EXIT_FAILURE);
								}
							}
						}
					}
				};

			ga.reset(config); // initialisation de l'algo génétique
			genetic::Individu<ConfigType> best = ga.run(verbose, callback);	// on lance l'algorithme en affichant des logs dans la console



			// *** ----> pas besoin de la suite pour le moment. Juste un exemple d'utilisation futur.
			/*
			/// changement de la configuration de l'algo génétique.
			config.enable_elitism = false;
			ga.reset(config);				// on met à jour la config pour l'algo génétique.
			ga.run(true, nullptr);			// on lance l'algo de nouveau.
			*/

			SimuCore::Systems::AdaptedSystem copy = system;
			auto [rocket, gen_state] = SimuCore::IndividualToRocket(best.to_real_vectors(), copy);
			return rocket;
		}
				
	}
}
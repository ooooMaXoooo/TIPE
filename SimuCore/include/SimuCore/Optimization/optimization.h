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
#include <DataExport/StatisticsData.h>
#include <DataExport/RocketData.h>

#include <SimuCore/statistics/HAC.h>


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
			const size_t population_size = 100, size_t max_generation = 5000,
			size_t print_interval=100, bool verbose=false,
			size_t snapshot_interval=10, bool saving_in_file=true,
			bool calculate_statistics=true)
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

			config.custom_mutation_proba = 1 * config.number_of_vectors * config.dimension * 0.003125;					// on prend la proba uniforme (-1 pour uniforme)



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
			using Ind = typename genetic::Individu<ConfigType>;
			using Pop = typename std::vector<Ind>;


			std::vector<std::unique_ptr<RocketData>> rockets_data;
			rockets_data.resize(config.population_size);
			// associe à un génome un score
			auto fitness = [&](const Ind& individu, size_t indice, bool last_evaluation) -> Real {
				// vec[0] --> un vecteur à 3 dimension, équivalent à la position initiale
				// vec[1] --> un vecteur à 3 dimension, équivalent à l'impulsion initiale
				// vec[2] --> un vecteur à 3 dimension, dont on n'utilise que la première composante, équivalent à l'instant de la 1ere impulsion
				// etc...
				std::vector<std::vector<Real>> vecs = individu.to_real_vectors();

				SimuCore::Systems::AdaptedSystem local_system = system; // on copie le système (problème de concurrence)
				double score = SimuCore::Systems::AdaptedSystem::m_LowestScore; // score par défaut pour les individus invalides
				auto [rocket, gen_state] = IndividualToRocket(vecs, local_system);

				if (last_evaluation) { // TODO : ajouter calculate_statistics en condition pour limiter les calculs lorsque c'est possible

					//local_system.Reset();		// On réinitialise le système pour que les individus commençent tous dans les mêmes confitions initiales.
					// On n'a pas besoin de réinitialiser car la fonction Score le fait déjà.

					rockets_data[indice] = std::make_unique<RocketData>(
						indice,
						-1.0,	// tof
						system.getMaxTime(), // maxTime
						rocket.GetInitialImpulsion(), // impulsion initiale
						rocket.position, // position initiale
						rocket.velocity, // vitesse initiale (toujours nulle)
						local_system.getFinalPlanetStartIndice(),
						static_cast<uint8_t>(local_system.getStartPlanetName()),
						static_cast<uint8_t>(local_system.getFinalPlanetName()),
						NbImpulsions,
						individu.to_integer_vectors()
					); 

					/// TODO : ajouter le stockage de l'individu génétique en forme binaire pour pouvoir le simuler et avoir sa trajectoire en python

					RocketData& current_rocket_data = *rockets_data[indice];


					auto position_callback = [&](const glm::dvec3& pos) {
						current_rocket_data.RegisterNewCartesianPosition(pos);
						};


					score = local_system.Score(rocket, gen_state, position_callback);
					double tof = local_system.getCurrentTime();
					current_rocket_data.SetTof(tof);

					const std::vector<std::pair<SimuCore::Structures::Impulsion, double>>& impulsions = rocket.getImpulsions();

					for (size_t i = 1; i < impulsions.size(); i++) {
						auto& [impulsion, time] = impulsions[i];

						if (time < tof)
							current_rocket_data.RegisterNewImpulsion(impulsion.GetDeltaV_vec(), time);
					}
				}
				else {
					score = local_system.Score(rocket, gen_state);
				}


				return score;
			};

			genetic::GeneticAlgorithm<ConfigType> ga(config, fitness); // création d'un algorithme génétique


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
				calculate_statistics,
				verbose,
				file_directory,
				&generationalExporter,
				&rockets_data
				]
				(
					size_t gen,
					double best_fit, int best_i,
					double worst_fit, int /*worst_i*/,
					const std::vector<genetic::Individu<ConfigType>>& population
					)
				{
					if (gen % print_interval == 0 || gen == max_generation - 1) {

						// best genes 
						{
							std::vector<std::vector<uint32_t>> genes = population[best_i].to_integer_vectors();
							
							std::cout << "Best genes :\n";
							for (size_t j = 0; j < genes.size(); ++j) {
								for (size_t k = 0; k < genes[j].size(); ++k) {
									std::cout << '\t' << genes[j][k] << ' ';
								}
								std::cout << '\n';
							}
						}

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
								best_fit < 8.33 * cste
								&&
								best_fit > 8 * cste
								)
							{

								double distance_to_final_planet_best = (best_fit - 8 * cste);
								distance_to_final_planet_best /= 0.33;
								distance_to_final_planet_best = 1 / distance_to_final_planet_best;
								distance_to_final_planet_best -= SimuCore::Systems::AdaptedSystem::m_CstScore;
								distance_to_final_planet_best = distance_to_final_planet_best * distance_to_final_planet_best * distance_to_final_planet_best;


								std::cout << "\tBest distance to target: "
									<< AU_to_kilometers(distance_to_final_planet_best) << " (km) = "
									<< distance_to_final_planet_best << " (AU)\n";
							} 
							else if (best_fit < 8.66 * cste && best_fit > 8.33 * cste) {
								std::cout << "\tBest distance to target: " << 0 << '\n';
								
								long double mecanic_energy = best_fit - 8.33 * cste; // km/s
								mecanic_energy /= 0.33;
								mecanic_energy = 1 / mecanic_energy;
								mecanic_energy -= SimuCore::Systems::AdaptedSystem::m_CstScore;
								mecanic_energy *= 1e24;
								mecanic_energy = std::sqrt(mecanic_energy);

								// calcul du minimum de l'énergie potentielle effective en valeur absolue, qu'on enlève à mecanic_energy

								std::cout << "\tEtat de diffusion\n";
								std::cout << "\tBest energy at target: " << mecanic_energy << " (J)\n";
							}
							else if (best_fit < 9 * cste && best_fit > 8.66 * cste) {
								std::cout << "\tBest distance to target: " << 0 << '\n';

								long double ratio = best_fit - 8.66 * cste;
								ratio /= (cste * 0.34);

								std::cout << "\tEtat lie\n";
								std::cout << std::fixed << std::setprecision(30);
								std::cout << "\tBest ratio power of tkt (energie meca / min energie pot effective) at target: " << ratio * 100 << "%\n";
								std::cout << std::setprecision(3) << std::defaultfloat;
							}
						}

						///////// Type de trajectoire
						std::cout << "\tKind of trajectory (best) : " << system.TypeOfTrajectory(best_fit) << '\n';

						///////// Info physique du meilleur individu
						{
							SimuCore::Systems::AdaptedSystem copy_system = system;
							auto [rocket, state] = SimuCore::IndividualToRocket(population[best_i].to_real_vectors(), copy_system);


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
						}

						if (gen % snapshot_interval == 0 || gen == max_generation - 1) {
							SimuCore::Systems::AdaptedSystem copy_system = system;
							///////// Envoi des données
							{
								// envoi des trajectoires et des impulsions du meilleur individu et des planètes de la génération dans un fichier
								{
									std::filesystem::path filepath =
										file_directory /
										generate_snapshot_filename("txt", gen + 1, "rocket");



									auto [rocket, state] = SimuCore::IndividualToRocket(population[best_i].to_real_vectors(), copy_system);
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

									// en jours
									double max_time = copy_system.getMaxTime();
									// en secondes
									double dt = copy_system.getDeltaTime();

									uint8_t startPlanetIndex = SimuCore::Systems::AdaptedSystem::GetStartPlanetID();
									uint8_t finalPlanetIndex = SimuCore::Systems::AdaptedSystem::GetFinalPlanetID();

									glm::dvec3 startPlanet_initialPosition = copy_system.GetStartPlanet_StartPosition();
									glm::dvec3 startPlanet_finalPosition = copy_system.GetStartPlanet_CurrentPosition();
									glm::dvec3 finalPlanet_initialPosition = copy_system.GetFinalPlanet_StartPosition();
									glm::dvec3 finalPlanet_finalPosition = copy_system.GetFinalPlanet_CurrentPosition();

									// en km/s
									glm::dvec3 finalVelocity_rocket = copy_system.GetRocketVelocity();


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

									bool etat_lie = false;

									// en km
									auto [r_min, r_max] = copy_system.GetApsidesAroundFinalPlanet(&etat_lie);


									r_min = meters_to_kilometers(r_min);
									r_max = meters_to_kilometers(r_max);

									constexpr uint8_t dimension = 2;


									double tof = trajectory.size() * dt / 86400.0; // en jours
									double delta_v = rocket.getDeltaV(tof);

									std::cout << "\tTime of flight : " << tof << " (jours)\n";
									std::cout << "\tDelta V : " << delta_v << " (km/s)\n";

									PhysicsData phys_data{ max_time, dt, startPlanetIndex, finalPlanetIndex,
										startPlanet_initialPosition,
										finalPlanet_initialPosition,
										startPlanet_finalPosition,
										finalPlanet_finalPosition,
										finalVelocity_rocket,
										impulsions_times,
										impulsions_vectors,
										tof, delta_v,
										etat_lie, r_min, r_max,
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


								///////////////////////////////////////////
								//                 Stats                 //
								///////////////////////////////////////////

								if (calculate_statistics) {
									/** On veut :
									* - le nombre de patates (clusters) de la population
									* - meilleur score de la population
									* - pire score de la population
									* - score moyen de la population
									* - le nombre de fusée ayant un même état, pour chaque état, i.e la taille de chaque état
									*/

									/**
									*
									*	 0 ~~> Collision start panet with high speed
									*	 1 ~~> Collision start panet with low speed
									*	 2 ~~> Collision sun with high speed
									*	 3 ~~> Collision sun with low speed
									*	 4 ~~> Acceleration too high
									*	 5 ~~> Rocket get too far in the system
									*	 6 ~~> Collision final panet with high speed
									*	 7 ~~> Collision final panet with low speed
									*	 8 ~~> Neutral : not in ring
									*	 9 ~~> Neutral : in ring but not in orbit
									*	10 ~~> Neutral : in orbit but too large
									*	11 ~~> Valid
									* 
									*/
									int kinds[12] = { 0 };

									auto type_of_trajectory = [&](double fitness) -> int {
										if (fitness < 0) return -1;

										constexpr double cste = 1.0 / SimuCore::Systems::AdaptedSystem::m_CstScore;

										if (fitness < 8*cste) {
											return static_cast<int>(fitness / cste);
										}
										else if (fitness < 8.33*cste) {
											return 8;
										}
										else if (fitness < 8.66*cste) {
											return 9;
										}
										else if (fitness < 9*cste) {
											return 10;
										}
										else {
											return 11;
										}
									};

									double score_mean = 0;

									for (const auto& ind : population) {

										if (ind.have_been_evaluated()) {
											double fitness = ind.get_fitness();
											score_mean += fitness;

											kinds[type_of_trajectory(fitness)]++;
										}
										else {
											assert(false && "Un individu n'a pas ete evalue, alors que tous les individus devraient l'etre à ce stade de l'algorithme genetique.");
											std::abort();
										}
									} // boucle for sur la population

									score_mean /= population.size();
									
									std::vector<SimuCore::Statistics::Cluster> clusters = SimuCore::Statistics::HAC(rockets_data);

									StatisticsData stats_data{clusters.size(), best_fit, worst_fit, score_mean, kinds};


									std::filesystem::path filepath_stats_data =
										file_directory /
										generate_snapshot_filename("txt", gen + 1, "stats");

									try {
										generationalExporter.enqueue(std::make_shared<StatisticsData>(stats_data), filepath_stats_data);
									}
									catch (const std::exception& e) {
										std::cerr << "\n\nTrajectory writing error ||\t" << e.what() << std::endl;
										exit(EXIT_FAILURE);
									}

								} // if calculate_statistics
							} // bloc envoi donnéees




							// envoi rocketsData
							{
								// trouvons le meilleur individu de la population
								size_t idx = 0;
								while (idx < rockets_data.size() && rockets_data[idx]->GetId() != best_i) {
									idx++;
								}


								std::filesystem::path filepath_rockets_data =
									file_directory /
									generate_snapshot_filename("txt", gen + 1, "rockets_data");
								try {
									generationalExporter.enqueue(std::make_shared<RocketData>(*rockets_data[idx]), filepath_rockets_data);
								}
								catch (const std::exception& e) {
									std::cerr << "\n\nRocketData writing error ||\t" << e.what() << std::endl;
									exit(EXIT_FAILURE);

								}
							}

						} // if snapshot
					}
					rockets_data.clear();
					rockets_data.resize(population.size());
				}; // callback

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
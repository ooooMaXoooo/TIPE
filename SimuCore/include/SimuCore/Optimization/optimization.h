#pragma once

#include <pch.h>
#include "../structures/Rocket.h"
#include "../structures/System.h"

#include <Galib/genetic.hpp>

#include <DataExport/AsyncDataExporter.h>

#include <DataExport/RocketData.h>
#include <DataExport/ClustersData.h>

#include <SimuCore/statistics/HAC.h>


namespace SimuCore {
	namespace Optimization {

		void writeTrajectory(std::ofstream& filestream, const std::vector<glm::dvec3>& trajectory);

		std::string generate_snapshot_directory();

		/// <summary>
		/// create a filename
		/// </summary>
		/// <returns> prefix + id + suffix + .extension </returns>
		std::string generate_snapshot_filename(const char* prefix = "gen_", size_t id = 0, const char* suffix = "", const std::string& extension = "txt");


		void sendPlanetsTrajectory(const std::filesystem::path& filepath_start, const std::filesystem::path& filepath_final, SimuCore::Systems::AdaptedSystem system);
		void sendIndividualTrajectory(const SimuCore::Structures::Rocket& rocket, GenerationState gen_state, const std::filesystem::path& filepath, SimuCore::Systems::AdaptedSystem* system, AsyncDataExporter& exporter);
		void sendRocketPhysics(const SimuCore::Structures::Rocket& rocket, const std::filesystem::path& filepath, SimuCore::Systems::AdaptedSystem* system, AsyncDataExporter& exporter);
		void sendStatistics(double best_fit, double worst_fit, double mean_score, int kinds[12], const std::vector<SimuCore::Statistics::Cluster>& clusters, const std::filesystem::path& filepath, AsyncDataExporter& exporter);



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

			config.custom_mutation_proba = 0.5 * 0.003125;					// on prend la proba uniforme (-1 pour uniforme)
			//config.custom_mutation_proba = 1 * config.number_of_vectors * config.dimension * 0.003125;					// on prend la proba uniforme (-1 pour uniforme)
			//config.custom_mutation_proba = -1;					// on prend la proba uniforme (-1 pour uniforme)



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
			auto fitness = [&](const Ind& individu, size_t indice, bool last_evaluation, int gen) -> Real {
				// vec[0] --> un vecteur à 3 dimension, équivalent à la position initiale
				// vec[1] --> un vecteur à 3 dimension, équivalent à l'impulsion initiale
				// vec[2] --> un vecteur à 3 dimension, dont on n'utilise que la première composante, équivalent à l'instant de la 1ere impulsion
				// etc...
				std::vector<std::vector<Real>> vecs = individu.to_real_vectors();

				SimuCore::Systems::AdaptedSystem local_system = system; // on copie le système (problème de concurrence)
				double score = SimuCore::Systems::AdaptedSystem::m_LowestScore; // score par défaut pour les individus invalides
				auto [rocket, gen_state] = IndividualToRocket(vecs, local_system);

				if (last_evaluation && (gen % snapshot_interval == 0)) { // TODO : ajouter calculate_statistics en condition pour limiter les calculs lorsque c'est possible

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
						constexpr double cste = 1.0 / SimuCore::Systems::AdaptedSystem::m_CstScore;

						///////// petites stats
						{
							size_t count_defined = 0;
							size_t count_valid = 0;
							for (auto& ind : population) {
								double fitness = ind.get_fitness();
								if (fitness > 9 * cste) {
									count_valid++;
									count_defined++;
								}
								else if (fitness > SimuCore::Systems::AdaptedSystem::m_LowestScore) {
									count_defined++;
								}
								else {
									throw std::runtime_error("un score trop bas existe");
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

						/*

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
					
						*/
					
					}

					if (saving_in_file) {
						// trajectoires des planètes de départ et d'arrivée
						if (gen == 0) {
							std::filesystem::path filepath_start = file_directory / "start_planet.txt";
							std::filesystem::path filepath_final = file_directory / "final_planet.txt";

							sendPlanetsTrajectory(filepath_start, filepath_final, system);
						}

						if (gen % snapshot_interval == 0 || gen == max_generation - 1) {
							// trajectoire meilleur individu
							{
								SimuCore::Systems::AdaptedSystem copy_system = system;
								auto [rocket, state] = SimuCore::IndividualToRocket(population[best_i].to_real_vectors(), copy_system);

								std::filesystem::path dir_best = file_directory / "Bests";
								std::filesystem::create_directories(dir_best);

								{
									std::filesystem::create_directories(dir_best / "Trajectories");
									std::filesystem::path filepath =
										dir_best / "Trajectories" /
										generate_snapshot_filename("gen_", gen+1, "", "traj");


									sendIndividualTrajectory(rocket, state, filepath, &copy_system, generationalExporter);
								}

								// envoi des données physiques
								{
									std::filesystem::create_directories(dir_best / "Physics");
									std::filesystem::path filepath_physics_data =
										dir_best / "Physics" /
										generate_snapshot_filename("gen_", gen+1, "", "phys");

									sendRocketPhysics(rocket, filepath_physics_data, &copy_system, generationalExporter);
								}

								// envoi rocketsData best
								{
									// trouvons le meilleur individu de la population
									size_t idx = 0;
									while (idx < rockets_data.size() && rockets_data[idx]->GetId() != best_i) {
										idx++;
									}


									std::filesystem::create_directories(dir_best / "RocketsData");
									std::filesystem::path filepath_rockets_data =
										dir_best / "RocketsData" /
										generate_snapshot_filename("gen_" , gen+1, "", "rck");

									try {
										generationalExporter.enqueue(std::make_shared<RocketData>(*rockets_data[idx]), filepath_rockets_data);
									}
									catch (const std::exception& e) {
										std::cerr << "\n\nRocketData writing error ||\t" << e.what() << std::endl;
										exit(EXIT_FAILURE);

									}
								} // rocketsData best
							}
								


							///////////////////////////////////////////
							//                 Stats                 //
							///////////////////////////////////////////

							if (calculate_statistics) {

								std::vector<SimuCore::Statistics::Cluster> clusters = SimuCore::Statistics::HAC(rockets_data);

								std::filesystem::path stats_dir = file_directory / "Stats";
								std::filesystem::create_directories(stats_dir);

								// generation's statistics
								{
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
									int kinds[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

									auto type_of_trajectory = [&](double fitness) -> int {
										if (fitness < 0) return -1;

										constexpr double cste = 1.0 / SimuCore::Systems::AdaptedSystem::m_CstScore;

										if (fitness < 8 * cste) {
											return static_cast<int>(fitness / cste);
										}
										else if (fitness < 8.33 * cste) {
											return 8;
										}
										else if (fitness < 8.66 * cste) {
											return 9;
										}
										else if (fitness < 9 * cste) {
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

											if (fitness > 0) {
												score_mean += fitness;

												int kind = type_of_trajectory(fitness);
												kinds[kind]++;

												if (kind < -0.5 || kinds[kind] < -0.5) {
													throw std::runtime_error("pas cool");
												}
											}
										}
										else {
											throw std::runtime_error("Un individu n'a pas ete evalue, alors que tous les individus devraient l'etre à ce stade de l'algorithme genetique.");
											std::abort();
										}
									} // boucle for sur la population

									score_mean /= population.size();

									std::filesystem::create_directories(stats_dir / "Simple");
									std::filesystem::path filepath_stats_data =
										stats_dir / "Simple" /
										generate_snapshot_filename("gen_", gen + 1, "", "stats");

									sendStatistics(best_fit, worst_fit, score_mean, kinds, clusters, filepath_stats_data, generationalExporter);
								} // generation's statistics


								// envoi Clusters
								{
									std::filesystem::path HAC_dir =
										stats_dir / "HAC";

									std::filesystem::create_directories(HAC_dir);

									std::filesystem::path filepath_cluster =
										HAC_dir /
										generate_snapshot_filename("gen_", gen + 1, "", "cluster");

									try {
										generationalExporter.enqueue(std::make_shared<ClustersData>(clusters), filepath_cluster);
									}
									catch (const std::exception& e) {
										std::cerr << "\n\nClustersData writing error ||\t" << e.what() << std::endl;
										exit(EXIT_FAILURE);
									}
								} // Clusters

							} // if calculate_statistics
							


							// envoi rocketsData all
							{
								std::filesystem::path file_directory_all_rockets =
									file_directory / "RocketsData" /  ("gen_" + std::to_string(gen + 1));

								std::filesystem::create_directories(file_directory_all_rockets);

								for (size_t idx = 0; idx < rockets_data.size(); idx++) {
									std::filesystem::path filepath_rocket_data =
										file_directory_all_rockets /
										generate_snapshot_filename("ind_", idx, "", "rck");

									try {
										generationalExporter.enqueue(std::make_shared<RocketData>(*rockets_data[idx]), filepath_rocket_data);
									}
									catch (const std::exception& e) {
										std::cerr << "\n\nRocketData writing error ||\t" << e.what() << std::endl;
										exit(EXIT_FAILURE);
									}
								}
							} // rocketsData all
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
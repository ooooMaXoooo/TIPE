#pragma once

#include <pch.h>
#include "../structures/Rocket.h"
#include "../structures/System.h"
#include "../utility.h"
#include <SimuCore\integrator\integrator.h>

#include <Galib/genetic.hpp>
#include <DataExport\Snapshot.h>
#include <DataExport\HDF5WriterAsync.h>
#include <DataExport\GenerationStats.h>
#include <DataExport\TrajectorySoA.h>
#include <DataExport\Accumulators.h>

namespace SimuCore {
	namespace Optimization {

		inline std::string generate_snapshot_filename(const std::string& extension = "h5") {
			// Récupérer le temps actuel
			auto now = std::chrono::system_clock::now();
			std::time_t t_now = std::chrono::system_clock::to_time_t(now);
			std::tm local_tm;

#if defined(_WIN32) || defined(_WIN64)
			localtime_s(&local_tm, &t_now);  // Windows
#else
			localtime_r(&t_now, &local_tm);  // Linux / macOS
#endif

			// Construire le nom sous la forme simu_jour_mois_annee_heure_minute_seconde.extension
			std::ostringstream oss;
			oss << "simu_"
				<< std::setfill('0') << std::setw(2) << local_tm.tm_mday << "_"
				<< std::setfill('0') << std::setw(2) << local_tm.tm_mon + 1 << "_"
				<< local_tm.tm_year + 1900 << "_"
				<< std::setfill('0') << std::setw(2) << local_tm.tm_hour << "_"
				<< std::setfill('0') << std::setw(2) << local_tm.tm_min << "_"
				<< std::setfill('0') << std::setw(2) << local_tm.tm_sec
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
		SimuCore::Structures::Rocket getBestRocket(const char* filename, const SimuCore::Systems::AdaptedSystem system,
			genetic::CrossoverType cross_type, bool elitism, bool auto_adapt,
			size_t population_size = 100, size_t max_generation = 5000,
			size_t print_interval=100, bool verbose=false, size_t snapshot_interval=10)
		{
			using ConfigType = genetic::Config<double, uint16_t, 2*NbImpulsions + 1, 3>;
			ConfigType config;

			config.dimension = 3;									// on travaille dans l'espace
			config.number_of_vectors = 2 * NbImpulsions + 1;		// un vecteur = un gène

			config.enable_saving = false;							// désactive la sauvegarde dans les fichiers car ce n'est pas encore implémenté

			config.population_size = population_size;				// paramètre à ajuster
			config.max_generations = max_generation;				// paramètre à ajuster
			config.print_interval  = print_interval;

			config.initial_mutation_probability = 0.2;				// paramètre à ajuster
			config.initial_self_adaptation_probability = 0.8;		// paramètre à ajuster



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


			std::string filename_save = generate_snapshot_filename("h5");
			HDF5WriterAsync snapshot_sink(filename_save); // nom de fichier adapté
			GenerationAccumulator accumulator(NbImpulsions);
			double last_best_score = -std::numeric_limits<double>::infinity();

			auto callback =
				[
					&accumulator,
					&system,
					snapshot_interval,
					&last_best_score,
					&snapshot_sink
				]
				(
					size_t gen,
					double best_fit, const auto& best_ind,
					double worst_fit, const auto& /*worst_ind*/,
					const auto& population
					)
				{
					if (gen % snapshot_interval != 0) return;

					// Détecter si on a un nouveau meilleur score
					bool is_new_best = std::abs(best_fit - last_best_score) > SimuCore::constants::epsilon;
					last_best_score = best_fit;

					// --- Accumuler population complète (avec validité) ---
					accumulator.push_population(population, system);

					// --- Finalisation du snapshot ---
					dataExport::Snapshot snapshot = accumulator.finalize(
						gen,
						dataExport::BestIndividualUpdate{} // sera mis à jour si nécessaire
					);

					// Mettre à jour les scores globaux
					snapshot.stats.best_score = best_fit;
					snapshot.stats.worst_score = worst_fit;

					// --- Mettre à jour la trajectoire du meilleur si nouveau meilleur ---
					if (is_new_best) {
						snapshot.best_update.is_new_best = true;
						snapshot.best_update.trajectory = calculate_trajectory(best_ind, system);
					}

					// --- Envoi vers le writer asynchrone ---
					snapshot_sink.enqueue(std::make_shared<dataExport::Snapshot>(std::move(snapshot)));
				};



			auto callback_2 =
				[print_interval,
				max_generation
				]
				(
					size_t gen,
					double best_fit, const auto& best_ind,
					double worst_fit, const auto& /*worst_ind*/,
					const auto& population
					)
				{
					if (gen % print_interval == 0 || gen == max_generation - 1) {
						size_t count_defined = 0;
						size_t count_valid = 0;
						for (auto& ind : population) {
							double fitness = ind.get_fitness();
							if (fitness > 0) {
								count_valid++;
								count_defined++;
							}
							else if (fitness > SimuCore::Systems::AdaptedSystem::m_LowestScore) {
								count_defined++;
							}
						}

						std::cout << "\tValid individuals: " << count_valid << "/" << population.size()
							<< " (" << (100.0 * static_cast<double>(count_valid) / population.size()) << "%)\n";
						std::cout << "\tDefined individuals: " << count_defined << "/" << population.size()
							<< " (" << (100.0 * static_cast<double>(count_defined) / population.size()) << "%)\n";

						constexpr double cste = 1.0 / SimuCore::Systems::AdaptedSystem::m_CstScore;
						if (
							best_fit < 7.751 * cste
							&& 
							best_fit > 7 * cste
							) 
						{

							double distance_to_final_planet_best = (best_fit - 7*cste);
							distance_to_final_planet_best /= 0.75;
							distance_to_final_planet_best = 1 / distance_to_final_planet_best;
							distance_to_final_planet_best -= SimuCore::Systems::AdaptedSystem::m_CstScore;

							std::cout << "\tBest distance to target: "
								<< AU_to_kilometers(distance_to_final_planet_best) << " (km) = "
								<< distance_to_final_planet_best << " (AU)\n";
						}
					}
				};

			ga.reset(config); // initialisation de l'algo génétique
			//genetic::Individu<ConfigType> best = ga.run(verbose, callback);	// on lance l'algorithme en affichant des logs dans la console
			genetic::Individu<ConfigType> best = ga.run(verbose, callback_2);	// on lance l'algorithme en affichant des logs dans la console



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
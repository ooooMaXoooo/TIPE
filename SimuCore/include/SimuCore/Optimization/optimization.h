#pragma once

#include <pch.h>
#include "../structures/Rocket.h"
#include "../structures/System.h"
#include "../utility.h"

#include <Galib/genetic.hpp>

namespace SimuCore {
	namespace Optimization {

		/// <summary>
		/// Calcul une fusée initialisée de façon à avoir la trajectoire optimale dans le système donné.
		/// </summary>
		/// <typeparam name="NbImpulsions"></typeparam>
		/// <param name="filename"> le chemin d'accès à un fichier de sauvergarde. </param>
		/// <param name="system"> le système dans lequel évolue une fusée. </param>
		/// <returns> </returns>
		template <size_t NbImpulsions>
		SimuCore::Structures::Rocket getBestRocket(const char* filename, SimuCore::Systems::AdaptedSystem system,
			genetic::CrossoverType cross_type, bool elitism, bool auto_adapt,
			size_t population_size = 100, size_t max_generation = 5000,
			size_t print_interval=100)
		{
			using ConfigType = genetic::Config<double, uint64_t, 2*NbImpulsions + 1, 3>;
			ConfigType config;

			config.dimension = 3;									// on travaille dans l'espace.
			config.number_of_vectors = 2 * NbImpulsions + 1;		// p0, v0

			config.enable_saving = false;							// désactive la sauvegarde dans les fichiers car ce n'est pas encore implémenté.

			config.population_size = population_size;				// paramètre à ajuster
			config.max_generations = max_generation;				// paramètre à ajuster
			config.print_interval  = print_interval;



			{
				// générer des réels tels que |x| < epaisseur de l'anneau  --> dans [-5, 5] pour éviter les problèmes de l'algo génétique

				config.max_real = system.RingSize_meter();
				config.min_real = -config.max_real;
			}

			config.enable_elitism = elitism;
			config.enable_auto_adaptation = auto_adapt;
			config.crossover_method = cross_type;

			using Real = ConfigType::real_type;
			using Integer = ConfigType::integer_type;
			auto fitness = [&](const std::vector<std::vector<Real>>& vecs) -> Real {
				// vec[0] --> un vecteur à 3 dimension, équivalent à la position initiale
				// vec[1] --> un vecteur à 3 dimension, équivalent à la vitesse initiale

				thread_local SimuCore::Systems::AdaptedSystem local_system = system;

				//local_system.Reset();					// On réinitialise le système pour que les individus commençent tous dans les mêmes confitions initiales.
				// On n'a pas besoin de réinitialiser car la fonction Score le fait déjà.
				return local_system.Score(vecs);		// A VERIFIER !!!!! OK ?
				};

			genetic::GeneticAlgorithm<ConfigType> ga(config, fitness); // création d'un algorithme génétique

			ga.reset(config);
			ga.run(true, nullptr);	// on lance l'algorithme en affichant des logs dans la console. Il n'y a pas encore de fonction de callback pour le moment.



			// *** ----> pas besoin de la suite pour le moment. Juste un exemple d'utilisation futur.
			/*
			/// changement de la configuration de l'algo génétique.
			config.enable_elitism = false;
			ga.reset(config);				// on met à jour la config pour l'algo génétique.
			ga.run(true, nullptr);			// on lance l'algo de nouveau.
			*/

			return SimuCore::Structures::Rocket(0, std::vector<std::pair<SimuCore::Structures::Impulsion, double>>(), 500, 5.4);
		}
				
	}
}
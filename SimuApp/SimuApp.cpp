#ifdef _OPENMP
#pragma message("OpenMP activé")
#else
#pragma message("OpenMP NON activé")
#endif

#include <omp.h>
#include <SimuCore.h>
#pragma comment(lib, "SimuCore.lib")

#include <SimuCore/Optimization/optimization.h>



static void rosenbrock_tests() {
    // Configuration pour Rosenbrock (minimum global en (1, 1))
    using ConfigType = genetic::Config<double, uint32_t, 1, 2>;
    ConfigType config;
    config.population_size = 500;    // NOLINT (cppcoreguidlines-avoid-magic-numbers)
    config.max_generations = 10000;  // NOLINT (cppcoreguidlines-avoid-magic-numbers)
    config.number_of_vectors = 1;   // Un seul vecteur
    config.dimension = 2;           // 2D: (x, y)
    config.min_real = -1e10;         // NOLINT (cppcoreguidlines-avoid-magic-numbers)
    config.max_real = 1e10;          // NOLINT (cppcoreguidlines-avoid-magic-numbers)
    config.integer_bits = 32;       // NOLINT (cppcoreguidlines-avoid-magic-numbers)
    config.print_interval = 1000;    // NOLINT (cppcoreguidlines-avoid-magic-numbers)
    config.tournament_size = 2;

    // Fonction de Rosenbrock: f(x,y) = (1-x)² + 100(y-x²)²
    // On veut minimiser, donc maximiser -f
    using Real = ConfigType::real_type;
    auto rosenbrock = [](const std::vector<std::vector<Real>>& vecs) -> Real {
        double x = vecs[0][0];
        double y = vecs[0][1];

        double term1 = (1.0 - x) * (1.0 - x);
        double term2 = 100.0 * (y - x * x) * (y - x * x);

        return -(term1 + term2);  // Négatif pour maximisation
        };
    genetic::GeneticAlgorithm<ConfigType> ga(config, rosenbrock);

    auto print_info_end_algo = [&ga]() -> void {
        // Le minimum global devrait être proche de (1, 1)
        auto best_vecs = ga.get_best_individual().to_real_vectors();
        std::cout << "Best solution found: (" << best_vecs[0][0] << ", " << best_vecs[0][1] << ")\n";
        std::cout << "Expected: (1.0, 1.0)\n";
        std::cout << "\n";
        };

    /// =========================== Cas 1 =========================== ///
    // Elitisme         : true
    // Auto-adaptation  : true
    // Cross over       : single point
    /// ============================================================= ///
    std::cout << "=== Cas 1 : elitism - auto adaptation - single point ===\n";
    config.enable_elitism = true;
    config.enable_auto_adaptation = true;
    config.crossover_method = genetic::CrossoverType::SINGLE_POINT_BIT_LEVEL;

    ga.reset(config);
    ga.run();
    print_info_end_algo();

    /// =========================== Cas 2 =========================== ///
    // Elitisme         : false
    // Auto-adaptation  : true
    // Cross over       : single point
    /// ============================================================= ///
    std::cout << "=== Cas 2 : no elitism - auto adaptation - single point ===\n";
    config.enable_elitism = false;
    config.enable_auto_adaptation = true;
    config.crossover_method = genetic::CrossoverType::SINGLE_POINT_BIT_LEVEL;

    ga.reset(config);
    ga.run();
    print_info_end_algo();

    /// =========================== Cas 3 =========================== ///
    // Elitisme         : true
    // Auto-adaptation  : false
    // Cross over       : single point
    /// ============================================================= ///
    std::cout << "=== Cas 3 : elitism - no auto adaptation - single point ===\n";
    config.enable_elitism = true;
    config.enable_auto_adaptation = false;
    config.crossover_method = genetic::CrossoverType::SINGLE_POINT_BIT_LEVEL;

    ga.reset(config);
    ga.run();
    print_info_end_algo();

    /// =========================== Cas 4 =========================== ///
    // Elitisme         : false
    // Auto-adaptation  : false
    // Cross over       : single point
    /// ============================================================= ///
    std::cout << "=== Cas 4 : no elitism - no auto adaptation - single point ===\n";
    config.enable_elitism = false;
    config.enable_auto_adaptation = false;
    config.crossover_method = genetic::CrossoverType::SINGLE_POINT_BIT_LEVEL;

    ga.reset(config);
    ga.run();
    print_info_end_algo();

    /// =========================== Cas 5 =========================== ///
    // Elitisme         : true
    // Auto-adaptation  : true
    // Cross over       : uniform
    /// ============================================================= ///
    std::cout << "=== Cas 5 : elitism - auto adaptation - uniform ===\n";
    config.enable_elitism = true;
    config.enable_auto_adaptation = true;
    config.crossover_method = genetic::CrossoverType::UNIFORM_BIT_LEVEL;

    ga.reset(config);
    ga.run();
    print_info_end_algo();

    /// =========================== Cas 6 =========================== ///
    // Elitisme         : false
    // Auto-adaptation  : true
    // Cross over       : uniform
    /// ============================================================= ///
    std::cout << "=== Cas 6 : no elitism - auto adaptation - uniform ===\n";
    config.enable_elitism = false;
    config.enable_auto_adaptation = true;
    config.crossover_method = genetic::CrossoverType::UNIFORM_BIT_LEVEL;

    ga.reset(config);
    ga.run();
    print_info_end_algo();

    // ============================ Cas 7 =========================== ///
    // Elitisme         : true
    // Auto-adaptation  : false
    // Cross over       : uniform
    /// ============================================================= ///
    std::cout << "=== Cas 7 : elitism - no auto adaptation - uniform ===\n";
    config.enable_elitism = true;
    config.enable_auto_adaptation = false;
    config.crossover_method = genetic::CrossoverType::UNIFORM_BIT_LEVEL;

    ga.reset(config);
    ga.run();
    print_info_end_algo();

    // ============================ Cas 8 =========================== ///
    // Elitisme         : false
    // Auto-adaptation  : false
    // Cross over       : uniform
    /// ============================================================= ///
    std::cout << "=== Cas 8 : no elitism - no auto adaptation - uniform ===\n";
    config.enable_elitism = false;
    config.enable_auto_adaptation = false;
    config.crossover_method = genetic::CrossoverType::UNIFORM_BIT_LEVEL;

    ga.reset(config);
    ga.run();
    print_info_end_algo();

    // ============================ Cas 9 =========================== ///
    // Elitisme         : true
    // Auto-adaptation  : false
    //      • taux de mutation fixe : 0.1
    // Cross over       : uniform
    /// ============================================================= ///
    std::cout << "=== Cas 9 : elitism - no auto adaptation - uniform - mutation rate (0.1) ===\n";
    config.enable_elitism = false;
    config.enable_auto_adaptation = false;
    config.custom_mutation_proba = 0.1;  // NOLINT (cppcoreguidlines-avoid-magic-numbers)
    config.crossover_method = genetic::CrossoverType::UNIFORM_BIT_LEVEL;

    ga.reset(config);
    ga.run();
    print_info_end_algo();

    // ============================ Cas 10 =========================== ///
    // Elitisme         : false
    // Auto-adaptation  : false
    //      • taux de mutation fixe : 0.1
    // Cross over       : uniform
    /// ============================================================= ///
    std::cout << "=== Cas 10 : no elitism - no auto adaptation - uniform - mutation rate (0.1) ===\n";
    config.enable_elitism = false;
    config.enable_auto_adaptation = false;
    config.custom_mutation_proba = 0.1;  // NOLINT (cppcoreguidlines-avoid-magic-numbers)
    config.crossover_method = genetic::CrossoverType::UNIFORM_BIT_LEVEL;

    ga.reset(config);
    ga.run();
    print_info_end_algo();

    /// =========================== Cas 11 =========================== ///
    // Elitisme         : true
    // Auto-adaptation  : true
    // Cross over       : uniform
    /// ============================================================= ///
    std::cout << "=== Cas 11 : elitism - auto adaptation - uniform ===\n";
    config.enable_elitism = true;
    config.enable_auto_adaptation = true;
    config.crossover_method = genetic::CrossoverType::UNIFORM_BIT_LEVEL;

    config.initial_mutation_probability = 0.01;         // NOLINT (cppcoreguidlines-avoid-magic-numbers)
    config.initial_self_adaptation_probability = 0.05;  // NOLINT (cppcoreguidlines-avoid-magic-numbers)

    ga.reset(config);
    ga.run();
    print_info_end_algo();
}





static void test_algo_bornes() {
    // Configuration pour Rosenbrock (minimum global en (1, 1))
    using ConfigType = genetic::Config<double, uint64_t, 1, 2>;
    ConfigType config;
    config.population_size = 60;    // NOLINT (cppcoreguidlines-avoid-magic-numbers)
    config.max_generations = 50000;  // NOLINT (cppcoreguidlines-avoid-magic-numbers)
    config.number_of_vectors = 1;   // Un seul vecteur
    config.dimension = 2;           // 2D: (x, y)
    config.min_real = -5;         // NOLINT (cppcoreguidlines-avoid-magic-numbers)
    config.max_real = 5;          // NOLINT (cppcoreguidlines-avoid-magic-numbers)
    config.integer_bits = 32;       // NOLINT (cppcoreguidlines-avoid-magic-numbers)
    config.print_interval = 1000;    // NOLINT (cppcoreguidlines-avoid-magic-numbers)
    config.tournament_size = 2;
	config.initial_self_adaptation_probability = 0.6;  // NOLINT (cppcoreguidlines-avoid-magic-numbers)

    // Fonction de Rosenbrock: f(x,y) = (1-x)² + 100(y-x²)²
    // On veut minimiser, donc maximiser -f
    using Real = ConfigType::real_type;
    auto rosenbrock = [](const std::vector<std::vector<Real>>& vecs) -> Real {
        double x = convertIntervals(-5, 5, -1e5, 1e5, vecs[0][0]);
        double y = convertIntervals(-5, 5, -1e5, 1e5, vecs[0][1]);

        double term1 = (1.0 - x) * (1.0 - x);
        double term2 = 100.0 * (y - x * x) * (y - x * x);

        return -(term1 + term2);  // Négatif pour maximisation
        };
    genetic::GeneticAlgorithm<ConfigType> ga(config, rosenbrock);

    auto print_info_end_algo = [&ga, &rosenbrock]() -> void {
        // Le minimum global devrait être proche de (1, 1)
        auto best_vecs = ga.get_best_individual().to_real_vectors();
        std::cout << "Best solution found: (" << convertIntervals(-5, 5, -1e5, 1e5, best_vecs[0][0]) << ", " <<
            convertIntervals(-5, 5, -1e5, 1e5, best_vecs[0][1]) << ")\n";
		std::cout << "\t~~> Score: " << rosenbrock(best_vecs) << '\n';
        std::cout << "Expected: (1.0, 1.0)\n";
        std::cout << "\n";
        };

    /// =========================== Cas 1 =========================== ///
    // Elitisme         : true
    // Auto-adaptation  : true
    // Cross over       : uniform
    /// ============================================================= ///
    std::cout << "=== Cas 1 : elitism - auto adaptation - uniform ===\n";
    config.enable_elitism = true;
    config.enable_auto_adaptation = true;
    config.crossover_method = genetic::CrossoverType::UNIFORM_BIT_LEVEL;

    ga.reset(config);
    ga.run();
       
    std::fflush(stdout);

    print_info_end_algo();

    /// =========================== Cas 2 =========================== ///
    // Elitisme         : true
    // Auto-adaptation  : false
    // Cross over       : uniform
    /// ============================================================= ///
    /*std::cout << "=== Cas 2 : elitism - no auto adaptation - uniform ===\n";
    config.enable_elitism = true;
    config.enable_auto_adaptation = false;
    config.crossover_method = genetic::CrossoverType::UNIFORM_BIT_LEVEL;

    ga.reset(config);
    ga.run();
    print_info_end_algo();*/
}




int main(int argc, char** argv) {
    
    /*genetic::print_info();
    std::cout << '\n';

    try {
        //rosenbrock_tests();
        test_algo_bornes();
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }

    std::cout << "Appuyez sur une touche pour continuer ...";
    std::cin.get();
    return 0;*/
    

    
	double lifetime = 500;
	SimuCore::Systems::AdaptedSystem sy(SimuCore::Systems::PlanetsName::Terre, SimuCore::Systems::PlanetsName::Mars, 0, 0, SimuCore::Structures::Rocket(lifetime, std::vector<std::pair<SimuCore::Structures::Impulsion, double>>(), 700000, 4.5), lifetime, 3600);

    genetic::CrossoverType cross_type = genetic::CrossoverType::UNIFORM_BIT_LEVEL;
    bool elitism = false;
    bool auto_adapt = true;
    size_t population_size = 100;
    size_t max_generation  =  10000;
    size_t print_interval  =   10;

	SimuCore::Optimization::getBestRocket<2>("", sy,
        cross_type, elitism, auto_adapt,
        population_size, max_generation, print_interval);
    

    //std::cout << convertIntervals(2, 3, -5, 5, 2.75) << '\n';
}
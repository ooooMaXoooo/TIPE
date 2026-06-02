#ifdef _OPENMP
#pragma message("OpenMP activé")
#else
#pragma message("OpenMP NON activé")
#endif

#include <omp.h>
#include <SimuCore.h>
#pragma comment(lib, "SimuCore.lib")

#include <SimuCore/Optimization/optimization.h>
#include <SimuCore/Units/all.h>


int main(int argc, char** argv) {
    
	double lifetime = 250; // durée de simulation en jours
	SimuCore::Systems::AdaptedSystem sy(
        SimuCore::Systems::PlanetsName::Terre,   // planète de départ
        SimuCore::Systems::PlanetsName::Mars,  // planète d'arrivée
        SimuCore::Structures::Rocket(
			lifetime, // -> durée de vie de la fusée en jours
            std::vector<std::pair<SimuCore::Structures::Impulsion, double>>(),
            1._ton_to_kg,
			2.22),                                     // -> vitesse d'éjection des gaz en km/s  (2.22 pour terre mars, 3 pour terre jupiter)
		lifetime,                                   // durée de simulation en jours
		3600);                                      // -> pas de temps en secondes

	sy.Initialize();

    genetic::CrossoverType cross_type = genetic::CrossoverType::UNIFORM_BIT_LEVEL; // ce paramètre ne change rien, on a implémenter en dur un UCLC
    bool elitism = true;                       // diminituion de la vitesse de perte de diversité ?
    bool auto_adapt = false;                    // a tester
    size_t population_size =  500;
    size_t max_generation  =  1000;
    size_t print_interval  =  1;
    bool verbose = true;
    size_t snapshot_interval = 1;
    bool save_in_file = true;
	bool calculate_statistics = true;

    constexpr size_t nombre_d_impulsions = 2;

	SimuCore::Optimization::getBestRocket<nombre_d_impulsions>(sy,
        cross_type, elitism, auto_adapt,
        population_size, max_generation,
        print_interval, verbose,
        snapshot_interval, save_in_file,
        calculate_statistics);
}
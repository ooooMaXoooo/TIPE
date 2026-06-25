#ifdef _OPENMP
#pragma message("OpenMP activé")
#else
#pragma message("OpenMP NON activé")
#endif

#include <omp.h>
#include <SimuCore.h>
#include <GaLib/genetic.hpp>
#pragma comment(lib, "SimuCore.lib")

#include <SimuCore/Optimization/optimization.h>
#include <SimuCore/Units/all.h>

void lambert_by_parts(bool is_jupiter) {
    constexpr double mu_terre = 3.986004418 * 1e14;
    constexpr double mu_mars = 4.282837 * 1e13;
    constexpr double mu_jupiter = 1.26686534 * 1e17;
    constexpr double mu_sun = 1.32712440042 * 1e20;

    const double distFinalSun = is_jupiter ? (5.2 * SimuCore::constants::AU) : (1.52 * SimuCore::constants::AU);


    constexpr double Nt1            = 2;
    constexpr double Nt2            = 2;
    constexpr double Nt3            = 2;
    constexpr double Nr1            = 2;
    constexpr double Nth1           = 2;
    constexpr double Nr2            = 2;
    constexpr double Nth2           = 2;
    constexpr double NthStart       = 2;
    constexpr double NthFinal       = 2;
    constexpr double NthFinalPlanet = 2;


    constexpr double r1_min = (6378 + 450) * 1e3;
    constexpr double r1_max = (6378 + 253111) * 1e3;

    constexpr double r2_min_mars = (3400 + 200) * 1e3;
    constexpr double r2_max_mars = (3400 + 126093) * 1e3;

    constexpr double r2_min_jupiter = (71500 + 2000) * 1e3;
    constexpr double r2_max_jupiter = (71500 + 23319047) * 1e3;

    constexpr double RStart = r1_max;

    constexpr double RFinal_mars = r2_max_mars;
    constexpr double RFinal_jupiter = r2_max_jupiter;

    const double mu_final = is_jupiter ? mu_jupiter : mu_mars;
    const double r2_min = is_jupiter ? r2_min_jupiter : r2_min_mars;
    const double r2_max = is_jupiter ? r2_max_jupiter : r2_max_mars;
    const double RFinal = is_jupiter ? RFinal_jupiter : RFinal_mars;

    constexpr double pi2 = 2 * SimuCore::constants::PI;

    using ConfigType = genetic::Config<long double, uint64_t, 8, 1>;
    ConfigType config;
    config.population_size = 1000;          // NOLINT (cppcoreguidlines-avoid-magic-numbers)
    config.max_generations = 100;          // NOLINT (cppcoreguidlines-avoid-magic-numbers)
    config.number_of_vectors = 8;           
    config.dimension = 1;                   // 2D: (x, y)
    config.min_real = 0;                    // NOLINT (cppcoreguidlines-avoid-magic-numbers)
    config.max_real = daysInSeconds(300);   // NOLINT (cppcoreguidlines-avoid-magic-numbers)
    //config.integer_bits = 64;               // NOLINT (cppcoreguidlines-avoid-magic-numbers)
    config.print_interval = 1;             // NOLINT (cppcoreguidlines-avoid-magic-numbers)
    config.tournament_size = 2;

    // On veut minimiser, donc maximiser -f
    using Real = ConfigType::real_type;
    using Integer = ConfigType::integer_type;
    using Ind = typename genetic::Individu<ConfigType>;
    using Pop = typename std::vector<Ind>;

    auto fitness_vector = [&](std::vector<std::vector<Real>> vecs) -> std::array<double, 17> {
        double tof1 = daysInSeconds(convertIntervals(config.min_real, config.max_real, 1, 15, vecs[0][0]));
        double tof2 = daysInSeconds(convertIntervals(config.min_real, config.max_real, 200, 300, vecs[1][0]));
        double tof3 = daysInSeconds(convertIntervals(config.min_real, config.max_real, 1, 15, vecs[2][0]));
        double theta1 = convertIntervals(config.min_real, config.max_real, 0.005, pi2 - 0.005, vecs[3][0]);
        double theta2 = convertIntervals(config.min_real, config.max_real, 0.005, pi2 - 0.005, vecs[4][0]);
        double thetaStart = convertIntervals(config.min_real, config.max_real, 0.005, pi2 - 0.005, vecs[5][0]);
        double thetaFinal = convertIntervals(config.min_real, config.max_real, 0.005, pi2 - 0.005, vecs[6][0]);
        double thetaFinalPlanet = convertIntervals(config.min_real, config.max_real, 0.005, pi2 - 0.005, vecs[7][0]);


        double tof1_min = tof1 - daysInSeconds(0.5);
        double tof1_max = tof1 + daysInSeconds(0.5);

        double tof2_min = tof2 - daysInSeconds(0.5);
        double tof2_max = tof2 + daysInSeconds(0.5);

        double tof3_min = tof3 - daysInSeconds(0.5);
        double tof3_max = tof3 + daysInSeconds(0.5);

        double theta1_min = theta1 - 0.005;
        double theta1_max = theta1 + 0.005;

        double theta2_min = theta2 - 0.005;
        double theta2_max = theta2 + 0.005;

        double thetaStart_min = thetaStart - 0.005;
        double thetaStart_max = thetaStart + 0.005;

        double thetaFinal_min = thetaFinal - 0.005;
        double thetaFinal_max = thetaFinal + 0.005;

        double thetaFinalPlanet_min = thetaFinalPlanet - 0.005;
        double thetaFinalPlanet_max = thetaFinalPlanet + 0.005;

        return orbit::lambert_by_parts(
            mu_terre, mu_final, mu_sun,
            distFinalSun,
            tof1_min, tof1_max,
            tof2_min, tof2_max,
            tof3_min, tof3_max,
            r1_min, r1_max,
            theta1_min, theta1_max,
            r2_min, r2_max,
            theta2_min, theta2_max,
            RStart,
            thetaStart_min, thetaStart_max,
            RFinal,
            thetaFinal_min, thetaFinal_max,
            thetaFinalPlanet_min, thetaFinalPlanet_max,
            Nt1, Nt2, Nt3,
            Nr1,
            Nth1,
            Nr2,
            Nth2,
            NthStart,
            NthFinal,
            NthFinalPlanet);
    };

    auto fitness = [&](const Ind& individu, size_t indice, bool last_evaluation, int gen) -> Real {
        std::vector<std::vector<Real>> vecs = individu.to_real_vectors();
        return -fitness_vector(vecs)[0];  // Négatif pour maximisation
    };
    genetic::GeneticAlgorithm<ConfigType> ga(config, fitness);

    config.enable_elitism = true;
    config.enable_auto_adaptation = true;
    config.crossover_method = genetic::CrossoverType::SINGLE_POINT_BIT_LEVEL;

    ga.reset(config);
    ga.run();


    std::vector<std::vector<Real>> vecs = ga.get_best_individual().to_real_vectors();

    std::array<double, 17> res = fitness_vector(vecs);

    std::cout << "delta_v : " << res[0] << " (m/s)\n";
    std::cout << "tof 1 : " << seconds_to_days(res[2]) << " days\n";
    std::cout << "tof 2 : " << seconds_to_days(res[3]) << " days\n";
    std::cout << "tof 3 : " << seconds_to_days(res[4]) << " days\n";
    std::cout << "tof : " << seconds_to_days(res[2] + res[3] + res[4]) << " days\n";

    std::cout << "thetaFinalPlanet : " << res[1] << " rad\n";
    std::cout << "thetaStart : " << res[5] << " rad\n";
    std::cout << "thetaFinal : " << res[6] << " rad\n";
    std::cout << "r1 : " << res[7] << " m\n";
    std::cout << "theta1 : " << res[8] << " rad\n";
    std::cout << "r2 : " << res[9] << " m\n";
    std::cout << "theta2 : " << res[10] << " rad\n";
    std::cout << "traj1_start : (" << res[11] << ", " << res[12] << ") m\n";
    std::cout << "traj2_start : (" << res[13] << ", " << res[14] << ") m\n";
    std::cout << "traj3_start : (" << res[15] << ", " << res[16] << ") m\n";

    std::cout << "RStart : " << RStart << " m\n";
    std::cout << "RFinal : " << RFinal << " m\n";
    std::cout << "distFinalSun : " << distFinalSun << " m\n";
    std::cout << "mu_final : " << mu_final << " USI\n";
}


int main(int argc, char** argv) {

    constexpr bool lambertByPart = true;
    constexpr double jupiter = true;

    if (lambertByPart) {
        lambert_by_parts(jupiter);
    }
    else {

        double lifetime = 258; // durée de simulation en jours
        SimuCore::Systems::AdaptedSystem sy(
            SimuCore::Systems::PlanetsName::Terre,   // planète de départ
            SimuCore::Systems::PlanetsName::Mars,    // planète d'arrivée
            SimuCore::Structures::Rocket(
                lifetime, // -> durée de vie de la fusée en jours
                std::vector<std::pair<SimuCore::Structures::Impulsion, double>>(),
                1._ton_to_kg,
                2.22),                                     // -> vitesse d'éjection des gaz en km/s  (2.22 pour terre mars, 3 pour terre jupiter)
            lifetime,                                   // durée de simulation en jours
            3600);                                      // -> pas de temps en secondes

        sy.Initialize();

        genetic::CrossoverType cross_type = genetic::CrossoverType::UNIFORM_BIT_LEVEL; // ce paramètre ne change rien, on a implémenter en dur un UCLC
        bool elitism = true;                        // diminituion de la vitesse de perte de diversité ?
        bool auto_adapt = false;                    // a tester
        size_t population_size = 10000;
        size_t max_generation = 150;
        size_t print_interval = 1;
        bool verbose = true;
        size_t snapshot_interval = 5;
        bool save_in_file = true;
        bool calculate_statistics = false;

        constexpr size_t nombre_d_impulsions = 2;

        SimuCore::Optimization::getBestRocket<nombre_d_impulsions>(sy,
            cross_type, elitism, auto_adapt,
            population_size, max_generation,
            print_interval, verbose,
            snapshot_interval, save_in_file,
            calculate_statistics);
    }
}
#pragma once

//#include <boost/type_index.hpp>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace genetic {

enum class CrossoverType : uint8_t {  // On spécifie le type pour l'optimisation mémoire
    SINGLE_POINT_BIT_LEVEL,
    UNIFORM_BIT_LEVEL
};

// todo Implémenter différents types de sélection
enum class SelectionType : uint8_t {  // On spécifie le type pour l'optimisation mémoire
    TOURNAMENT,
    ROULETTE_WHEEL
};

// Concept définissant les requis d'une configuration valide
template <typename T>
concept ConfigConcept = requires {
    // Types requis
    typename T::real_type;
    typename T::integer_type;
    requires std::is_floating_point_v<typename T::real_type>;
    requires std::is_unsigned_v<typename T::integer_type>;

    // Constantes template et leurs contraintes
    { T::max_vectors } -> std::convertible_to<size_t>;
    { T::max_dimension } -> std::convertible_to<size_t>;
    requires T::max_vectors > 0;
    requires T::max_dimension > 0;
} && requires(const T& config) {
    // Paramètres de représentation des nombres
    { config.min_real } -> std::convertible_to<typename T::real_type>;
    { config.max_real } -> std::convertible_to<typename T::real_type>;
    { config.integer_bits } -> std::convertible_to<size_t>;
    { config.integer_max } -> std::convertible_to<typename T::integer_type>;

    // Paramètres de l'algorithme génétique
    { config.dimension } -> std::convertible_to<size_t>;
    { config.enable_elitism } -> std::convertible_to<bool>;
    { config.population_size } -> std::convertible_to<size_t>;
    { config.max_generations } -> std::convertible_to<size_t>;
    { config.tournament_size } -> std::convertible_to<size_t>;
    { config.number_of_vectors } -> std::convertible_to<size_t>;
    { config.enable_auto_adaptation } -> std::convertible_to<bool>;
    { config.crossover_method } -> std::convertible_to<CrossoverType>;
    { config.initial_mutation_probability } -> std::convertible_to<typename T::real_type>;
    { config.initial_self_adaptation_probability } -> std::convertible_to<typename T::real_type>;

    // Paramètre spécifique à la méthode de coisement uniforme
    { config.uniform_crossover_probability } -> std::convertible_to<double>;

    // Paramètres de sauvegarde
    { config.enable_saving } -> std::convertible_to<bool>;
    { config.save_interval } -> std::convertible_to<size_t>;
    { config.print_interval } -> std::convertible_to<size_t>;
    { config.save_directory } -> std::convertible_to<const char*>;
    { config.save_extension } -> std::convertible_to<const char*>;
    { config.save_config_extension } -> std::convertible_to<const char*>;

    // Méthodes requises
    { config.validate() } -> std::same_as<void>;
    { config.get_half_population_size() } -> std::convertible_to<size_t>;
    { config.get_real_size() } -> std::convertible_to<typename T::real_type>;
    { config.get_fixed_proba() } -> std::convertible_to<typename T::real_type>;
    { config.get_integer_max() } -> std::convertible_to<typename T::integer_type>;

    // Variables utiles
    { T::proba_distribution } -> std::convertible_to<std::uniform_real_distribution<typename T::real_type>>;
};

/**
 * @brief Configuration structure for the genetic algorithm
 *
 * @tparam Real Floating-point type for real-valued genes (float, double, long double)
 * @tparam Integer Unsigned integer type for integer-valued genes (uint8_t, uint16_t, uint32_t, uint64_t)
 * @tparam MaxVectors Maximum number of vectors per individual
 * @tparam MaxDimension Maximum dimension per vector
 */
template <typename Real, typename Integer, size_t MaxVectors = 2, size_t MaxDimension = 3>
struct Config {
    // Validation des types à la compilation
    static_assert(std::is_floating_point_v<Real>,
                  "Real type must be a floating-point type (float, double, long double)");
    static_assert(std::is_unsigned_v<Integer>,
                  "Integer type must be an unsigned integer type (uint8_t, uint16_t, uint32_t, uint64_t)");

    static_assert(MaxVectors > 0, "There must be at least 1 vector, e.g at least 1 chromosome");
    static_assert(MaxDimension > 0, "The dimension must be at least 1, e.g there is at least 1 gene per chromosome");

    // ========== Alias des types ==========

    using real_type = Real;
    using integer_type = Integer;

    // ========== Paramètres de représentation des nombres ==========

    /// Valeur minimale pour les gènes réels
    Real min_real = Real(-1000);  // NOLINT (cppcoreguidelines-avoid-magic-numbers)

    /// Valeur maximale pour les gènes réels
    Real max_real = Real(1000);  // NOLINT (cppcoreguidelines-avoid-magic-numbers)

    /// Nombre de bits utilisés pour encoder un gène entier (doit être <= sizeof(Integer)*8)
    size_t integer_bits = std::numeric_limits<Integer>::digits;

    /// Valeur maximale pour les gènes entiers (calculée automatiquement si laissée à 0)
    Integer integer_max = 0;

    // ========== Paramètres de l'algorithme génétique ==========

    static constexpr size_t max_vectors = MaxVectors;
    static constexpr size_t max_dimension = MaxDimension;

    /// booléant permettant l'élitisme
    bool enable_elitism = true;

    /// booléant permettant l'élitisme
    bool enable_auto_adaptation = true;

    /// Probabilité de mutation si l'auto-adaptation est désactivée (si négative, on prendra la proba uniforme)
    double custom_mutation_proba = -1;

    /// Taille de la population (doit être un nombre pair)
    size_t population_size = 1000;  // NOLINT (cppcoreguildelines-avoid-magic-numbers)

    /// Nombre de vecteurs (chromosomes) par individu
    size_t number_of_vectors = max_vectors;

    /// Dimension de chaque vecteur (nombre de gènes par chromosome)
    size_t dimension = max_dimension;

    /// Nombre maximum de générations
    size_t max_generations = 1000;  // NOLINT (cppcoreguidelines-avoid-magic-numbers)

    /// Probabilité de mutation initiale (entre 0 et 1)
    Real initial_mutation_probability = Real(0.01);  // NOLINT (cppcoreguidelines-avoid-magic-numbers)

    /// Probabilité de mutation de chaque bit de chaque gène du chromosome de mutation
    Real initial_self_adaptation_probability = 0.05;  // NOLINT (cppcoreguidelines-avoid-magic-numbers)

    /// Taille du tournoi pour la sélection
    size_t tournament_size = 2;

    /// Méthode de crossover à utiliser
    CrossoverType crossover_method = CrossoverType::SINGLE_POINT_BIT_LEVEL;

    // ========== Paramètre spécifique à la méthode de coisement uniforme ==========

    /// probabilité de donner le bit du parent 1 à l'enfant 1
    double uniform_crossover_probability = 0.5;  // NOLINT (cppcoreguidelines-avoid-magic-numbers)

    // ========== Paramètres de sauvegarde ==========

    /// Active la sauvegarde des générations
    bool enable_saving = false;

    /// Intervalle de générations entre deux sauvegardes
    size_t save_interval = 10;  // NOLINT (cppcoreguidelines-avoid-magic-numbers)

    /// Intervalle de générations entre deux affichages
    size_t print_interval = 10;  // NOLINT (cppcoreguidelines-avoid-magic-numbers)

    /// Répertoire de sauvegarde
    const char* save_directory = "./data";

    /// Extension des fichiers de sauvegarde pour les générations
    const char* save_extension = ".gen";

    /// Extension des fichiers de sauvegarde pour la configuration
    const char* save_config_extension = ".config";

    // ========== Variable pour utile pour la place en mémoire ==========
    inline static std::uniform_real_distribution<Real> proba_distribution = std::uniform_real_distribution<Real>(0, 1);

    /**
     * @brief Default constructor with validation
     */
    Config() { validate(); }

    /**
     * @brief Validates the configuration parameters
     * @throws std::invalid_argument if parameters are invalid
     */
    void validate() const {
        // Validation de la plage des réels
        if (min_real >= max_real) {
            throw std::invalid_argument("min_real must be less than max_real");
        }

        // Validation des bits pour les entiers
        if (integer_bits == 0 || integer_bits > std::numeric_limits<Integer>::digits) {
            throw std::invalid_argument("integer_bits must be between 1 and " +
                                        std::to_string(std::numeric_limits<Integer>::digits));
        }

        // Validation des dimensions
        if (number_of_vectors == 0 || number_of_vectors > max_vectors) {
            throw std::invalid_argument("number_of_vectors must lie between 1 and " + std::to_string(max_vectors));
        }

        if (dimension == 0 || dimension > max_dimension) {
            throw std::invalid_argument("dimension must lie between 1 and " + std::to_string(max_dimension));
        }

        // Validation des probabilités de mutation
        if (initial_mutation_probability < 0 || initial_mutation_probability > 1) {
            throw std::invalid_argument("initial_mutation_probability must be between 0 and 1");
        }

        if (initial_self_adaptation_probability < 0 || initial_self_adaptation_probability > 1) {
            throw std::invalid_argument("initial_self_adaptation_probability must be between 0 and 1");
        }

        // Validation du nombre de générations
        if (max_generations == 0) {
            throw std::invalid_argument("max_generations must be non-null");
        }

        // Validation de la taille du tournoi
        if (tournament_size == 0) {
            throw std::invalid_argument("tournament_size must be non-null");
        }

        // validation de la taille de la population
        if (population_size % 2 == 1) {
            throw std::invalid_argument("population size must be an even number (for pairing)");
        }
        if (population_size == 0) {
            throw std::invalid_argument("population size must non null");
        }
    }

    /**
     * @brief Gets the computed real size (max - min)
     */
    [[nodiscard]] constexpr Real get_real_size() const noexcept { return max_real - min_real; }

    /**
     * @brief Gets the computed maximum value for integers based on bits
     */
    [[nodiscard]] Integer get_integer_max() const noexcept {
        if (integer_max > 0) {
            return integer_max;
        }

        if (integer_bits == std::numeric_limits<Integer>::digits) {
            return std::numeric_limits<Integer>::max();
        }

        return (Integer(1) << integer_bits) - Integer(1);
    }

    /**
     * @brief Gets half the population size
     */
    [[nodiscard]] constexpr size_t get_half_population_size() const noexcept { return population_size / 2; }

    [[nodiscard]] Real get_fixed_proba() const {
        if (custom_mutation_proba < 0) {
            return 1 / static_cast<Real>(integer_bits * dimension * number_of_vectors);
        }

        return custom_mutation_proba;
    }

    /**
     * @brief Output operator for debugging
     */
    friend std::ostream& operator<<(std::ostream& os, const Config& config) {
        os << "=============== Configuration settings ===============\n";
        /*os << "Real type: " << boost::typeindex::type_id<Real>().pretty_name() << '\n';
        os << "Integer type: " << boost::typeindex::type_id<Integer>().pretty_name() << '\n';*/

        os << "\n\t=== Numbers representation ===\n";
        os << "Reals = [" << config.min_real << ", " << config.max_real << "]\n";
        os << "Integers are represented on " << config.integer_bits << " bits\n";
        os << "Integer maximum value: " << config.get_integer_max() << '\n';

        os << "\n\t=== Genetic settings ===\n";
        os << "Population size: " << config.population_size << '\n';
        os << "Max generations: " << config.max_generations << '\n';
        os << "Number of vectors: " << config.number_of_vectors << '\n';
        os << "Dimension: " << config.dimension << '\n';
        os << "Max number of vectors: " << Config::max_vectors << '\n';
        os << "Max dimension: " << Config::max_dimension << '\n';
        os << "Tournament size: " << config.tournament_size << '\n';
        os << "Elitism: " << (config.enable_elitism ? "true" : "false") << '\n';
        os << "Auto-adaptation: " << (config.enable_auto_adaptation ? "true" : "false") << '\n';
        if (!config.enable_auto_adaptation) {
            os << "\tProbability across all chromosomes: " << config.get_fixed_proba() << '\n';
        }
        os << "Initial mutation probability: " << config.initial_mutation_probability << '\n';
        os << "Initial self adaptation probability: " << config.initial_self_adaptation_probability << '\n';

        os << "Crossover method: ";
        switch (config.crossover_method) {
            case genetic::CrossoverType::SINGLE_POINT_BIT_LEVEL:
                os << "single point bit level\n";
                break;
            case genetic::CrossoverType::UNIFORM_BIT_LEVEL:
                os << "uniform bit level\n";
                os << "Uniform crossover swap probability: " << config.uniform_crossover_probability << '\n';
                break;
        }

        os << "\n\t=== IO settings ===\n";
        os << "Interval between consecutive print: " << config.print_interval << '\n';
        os << "Is saving enable " << (config.enable_saving ? "true" : "false") << '\n';
        if (!config.enable_saving) {
            return os;
        }
        os << "Interval between consecutive save: " << config.save_interval << '\n';
        os << "Save directory: " << config.save_directory << '\n';
        os << "Generation file extension: " << config.save_extension << '\n';
        os << "Configuration file extension: " << config.save_config_extension << '\n';

        return os;
    }
};

// Vérifications statiques que Config satisfait son propre concept
static_assert(ConfigConcept<Config<double, uint32_t>>, "Config must satisfy ConfigConcept");
// NOLINTNEXTLINE (cppcoreguidlines-avoid-magic-numbers)
static_assert(ConfigConcept<Config<float, uint8_t, 5, 10>>, "Config must satisfy ConfigConcept with custom parameters");

}  // namespace genetic
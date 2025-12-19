#pragma once

#include <array>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

#include "genetic_config.hpp"
#include "genetic_utils.hpp"

namespace genetic {

/**
 * @brief Represents an individual in the genetic algorithm
 *
 * @tparam ConfigType the specific type of a config (for static allocation)
 */
template <ConfigConcept ConfigType>
class Individu {
public:
    using Real = typename ConfigType::real_type;
    using Integer = typename ConfigType::integer_type;

    using gene = Integer;
    using chromosome = std::array<gene, ConfigType::max_dimension>;
    using genome = std::array<chromosome, ConfigType::max_vectors>;
    using proba_array = std::array<Integer, ConfigType::max_vectors + 1>;

private:
    genome m_genome;
    proba_array m_mutation_probas;
    const ConfigType* m_config;

    // Dimensions actuelles (peuvent être < Max)
    size_t m_num_vectors;
    size_t m_dimension;

    Real m_fitness;
    bool m_have_been_evaluated = false;

public:
    /**
     * @brief Constructs an individual with random genes
     * @tparam RNG The type of the random number generator
     * @param config Configuration reference
     * @param rng The random number generator instance
     */
    template <typename RNG>
    Individu(const ConfigType& config, RNG& rng)
        : m_config(&config),
          m_num_vectors(config.number_of_vectors),
          m_dimension(config.dimension),
          m_fitness(std::numeric_limits<Real>::lowest()) {
        if (m_num_vectors > ConfigType::max_vectors) {
            throw std::invalid_argument("number_of_vectors exceeds max_vectors");
        }
        if (m_dimension > ConfigType::max_dimension) {
            throw std::invalid_argument("dimension exceeds max_dimension");
        }

        std::uniform_int_distribution<Integer> gene_dist(0, config.get_integer_max());

        // Initialiser les gènes aléatoirement
        for (size_t i = 0; i < m_num_vectors; i++) {
            for (size_t j = 0; j < m_dimension; j++) {
                m_genome[i][j] = gene_dist(rng);
            }
        }

        // Initialiser les probabilités de mutation
        Integer proba_bin =
            utils::proba_to_bin<Real, Integer>(config.initial_mutation_probability, config.integer_bits);

        for (size_t i = 0; i < m_num_vectors; i++) {
            m_mutation_probas[i] = proba_bin;
        }

        // Initialiser la probabilité d'AUTO-ADAPTATION
        Integer self_proba_bin =
            utils::proba_to_bin<Real, Integer>(config.initial_self_adaptation_probability, config.integer_bits);
        m_mutation_probas[m_num_vectors] = self_proba_bin;
    }

    /**
     * @brief Default constructor for container compatibility
     */
	Individu() : m_config(nullptr), m_num_vectors(0), m_dimension(0), m_fitness(std::numeric_limits<Real>::lowest()) {}

    // ========== Getters ==========

    [[nodiscard]] const genome& get_genome() const { return m_genome; }
    [[nodiscard]] size_t get_num_vectors() const { return m_num_vectors; }
    [[nodiscard]] size_t get_dimension() const { return m_dimension; }

    [[nodiscard]] gene get_gene(size_t chromo, size_t gene_idx) const { return m_genome.at(chromo).at(gene_idx); }

    [[nodiscard]] const chromosome& get_chromosome(size_t index) const { return m_genome.at(index); }

    [[nodiscard]] Integer get_mutation_proba(size_t index) const { return m_mutation_probas.at(index); }
    [[nodiscard]] proba_array get_mutation_probas() const { return m_mutation_probas; }

    [[nodiscard]] Real get_fitness() const { return m_fitness; }
    [[nodiscard]] bool have_been_evaluated() const { return m_have_been_evaluated; }

    [[nodiscard]] uint8_t get_bit(size_t chromo, size_t gene_idx, size_t bit_idx) const {
        const Integer mask = (Integer(1) << bit_idx);
        return (get_gene(chromo, gene_idx) & mask) >> bit_idx;
    }

    // ========== Setters ==========

    void set_gene(size_t chromo, size_t gene_idx, Integer value) { m_genome[chromo][gene_idx] = value; }

    void set_chromosome(size_t index, const chromosome& c) { m_genome[index] = c; }

    void set_mutation_proba(size_t index, Integer value) { m_mutation_probas[index] = value; }

    void set_config(const ConfigType* config) {
        m_config = config;
        m_num_vectors = config->number_of_vectors;
        m_dimension = config->dimension;
    }

    void set_fitness(Real fitness) {
        m_fitness = fitness;
        m_have_been_evaluated = true;
    }

    // ========== Genetic operations ==========

    /**
     * @brief Performs mutation on this individual
     * @tparam rng Random number generator
     */
    template <typename RNG>
    void mutate(RNG& rng) {
        mutate_vectors(rng);
        mutate_proba(rng);

        invalidate_fitness();
    }

    /**
     * @brief Converts genome to real-valued vectors
     * @return Vector of vectors containing the decoded real values
     */
    std::vector<std::vector<Real>> to_real_vectors() const {
        if (m_config == nullptr) {
            return {};
        }

        std::vector<std::vector<Real>> result(m_num_vectors, std::vector<Real>(m_dimension));

        for (size_t i = 0; i < m_num_vectors; i++) {
            for (size_t j = 0; j < m_dimension; j++) {
                /*result[i][j] = utils::bin_to_real<Real, Integer>(m_genome[i][j], m_config->min_real, m_config->max_real,
                                                                 m_config->integer_bits);*/
                result[i][j] = utils::bin_to_real(m_genome[i][j], m_config->min_real, m_config->max_real,
                    m_config->integer_bits);
            }
        }

        return result;
    }

    /**
     * @brief Output operator for debugging
     */
    friend std::ostream& operator<<(std::ostream& os, const Individu& ind) {
        if (ind.m_config == nullptr) {
            os << "[Uninitialized individual]";
            return os;
        }

        auto real_vecs = ind.to_real_vectors();
        for (size_t i = 0; i < ind.m_num_vectors; i++) {
            os << "[";
            for (size_t j = 0; j < ind.m_dimension; j++) {
                os << real_vecs[i][j];
                if (j < ind.m_dimension - 1) {
                    os << ", ";
                }
            }
            os << "]";
            if (i < ind.m_num_vectors - 1) {
                os << "\n";
            }
        }
        return os;
    }

    /**
     * @brief Invalide la fitness après une modification génétique (mutation/crossover)
     */
    void invalidate_fitness() {
        m_have_been_evaluated = false;
        // réinitialiser la fitness à la pire valeur
        if (m_config) {
            m_fitness = std::numeric_limits<Real>::lowest();
        }
    }

private:
    /**
     * @brief Performs mutation on this individual's mutation chromosome
     * @tparam rng Random number generator
     */
    template <typename RNG>
    void mutate_proba(RNG& rng) {
        // Ce bloc ne s'exécute que si l'auto-adaptation est activée.
        if (!m_config || !m_config->enable_auto_adaptation) {
            return;
        }

        Real self_mutation_proba =
            utils::bin_to_proba<Real, Integer>(m_mutation_probas[m_num_vectors], m_config->integer_bits);
        /*#pragma omp critical
        {
            std::cout << "mutate individu, proba : " << self_mutation_proba << '\n';
        }*/
        mutate_single_array(rng, self_mutation_proba, m_mutation_probas, m_num_vectors + 1);
    }

    /**
     * @brief Performs mutation on this individual's vectors chromosomes
     * @tparam rng Random number generator
     */
    template <typename RNG>
    void mutate_vectors(RNG& rng) {
        if (!m_config) {
            return;
        }

        if (m_config->enable_auto_adaptation) {
            for (size_t i = 0; i < m_num_vectors; i++) {
                // Taux calculé à chaque vecteur i
                Real current_rate = utils::bin_to_proba<Real, Integer>(m_mutation_probas[i], m_config->integer_bits);
                mutate_single_array(rng, current_rate, m_genome[i], m_dimension);
            }
        } else {
            // Calcule le taux fixe une seule fois (ou le lit depuis la config)
            Real fixed_rate = m_config->get_fixed_proba();

            for (size_t i = 0; i < m_num_vectors; i++) {
                // Le taux 'fixed_rate' est constant pour tous les vecteurs
                mutate_single_array(rng, fixed_rate, m_genome[i], m_dimension);
            }
        }
    }

    /**
     * @brief Performs mutation on a single chromosome
     * @tparam rng Random number generator
     * @tparam N
     */
    template <typename RNG, size_t N>
    void mutate_single_array(RNG& rng, Real rate, std::array<Integer, N>& target_array, size_t array_size) {
        for (size_t locus = 0; locus < array_size; locus++) {
            for (size_t bit = 0; bit < m_config->integer_bits; ++bit) {
                if (ConfigType::proba_distribution(rng) <= rate) {
                    target_array[locus] ^= (Integer(1) << bit);
                }
            }
        }
    }
};

}  // namespace genetic
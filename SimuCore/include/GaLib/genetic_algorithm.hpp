#pragma once

#include <algorithm>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include "genetic_config.hpp"
#include "genetic_individu.hpp"
#include "genetic_utils.hpp"


#ifdef _OPENMP
#include <omp.h>
#endif



namespace genetic {

/**
 * @brief Main genetic algorithm class
 *
 * @tparam ConfigType the specific type of a config
 */
template <ConfigConcept ConfigType>
class GeneticAlgorithm {
public:
    using Real = typename ConfigType::real_type;
    using Integer = typename ConfigType::integer_type;

    using Individual = Individu<ConfigType>;
    using Population = std::vector<Individual>;

    /**
     * @brief Fitness function type
     * Takes vectors of real values and returns a fitness score (higher is better)
     */
    using FitnessFunction = std::function<Real(const std::vector<std::vector<Real>>&)>;

private:
    ConfigType m_config;
    FitnessFunction m_fitness_func;
    unsigned int m_seed;
    std::mt19937_64 m_rng;
    Population m_population;

    Real m_best_fitness;
    Individual m_best_individual;

    Real m_worst_fitness;
    Individual m_worst_individual;

    // section pour la fonction selection
    std::uniform_int_distribution<size_t> m_index_dist;
    Population m_selected;

    // section commune aux fonctions suivantes : crossover, bit_level_crossover, bit_level_crossover_probas
    std::uniform_int_distribution<size_t> m_cut_dist;
    std::uniform_int_distribution<size_t> m_cut_dist_proba;

public:
    /**
     * @brief Constructs the genetic algorithm
     *
     * @param config Configuration parameters
     * @param fitness_func Fitness function to maximize
     * @param seed Random seed (0 for random seed)
     */
    GeneticAlgorithm(const ConfigType& config, FitnessFunction fitness_func, unsigned int seed = 0)
        : m_config(config),
          m_fitness_func(fitness_func),
          m_seed(seed),
          m_rng(seed == 0 ? std::random_device{}() : seed),
          m_best_fitness(std::numeric_limits<Real>::lowest()),
          m_worst_fitness(std::numeric_limits<Real>::max()),
          m_index_dist(0, config.population_size - 1),
          m_cut_dist(0, (config.dimension * config.integer_bits) - 1),
          m_cut_dist_proba(0, ((config.number_of_vectors + 1) * config.integer_bits) - 1) {
        m_config.validate();
        initialize_population();
        m_selected.resize(m_config.get_half_population_size());
    }

    /**
     * @brief Runs the complete genetic algorithm
     *
     * @param verbose Print progress information
     * @param callback Optional callback called after each generation
     */
    void run(bool verbose = true, std::function<void(size_t, Real, const Individual&)> callback = nullptr) {
        if (verbose) {
            std::cout << "Starting genetic algorithm..." << '\n';
            std::cout << m_config;
            std::cout << '\n';
        }

        auto affiche_proba_array =
            [&](const Individual& ind)
            {
                std::cout << "[";
                auto probas = ind.get_mutation_probas();
                for (int i = 0; i < probas.size(); i++) {
                    std::cout << genetic::utils::bin_to_proba<Real>(probas[i], m_config.integer_bits);
                    if (i < probas.size() - 1) {
                        std::cout << ", ";
                    }
                }
                std::cout << ']';
            };

        for (size_t gen = 0; gen < m_config.max_generations; gen++) {
            step(gen);

            if (verbose && (gen % m_config.print_interval == 0 || gen == m_config.max_generations - 1)) {
                std::cout << "Generation " << gen + 1 << "/" << m_config.max_generations;

                // affichage best
                {
                    std::cout << "\n\tBest fitness: " << m_best_fitness;
                    if (m_config.enable_auto_adaptation) {
                        std::cout << " ~ Proba Array: ";
                        affiche_proba_array(m_best_individual);
                    }
                }

                // affichage worst
                {
                    std::cout << "\n\tWorst fitness: " << m_worst_fitness;
                    if (m_config.enable_auto_adaptation) {
                        std::cout << " ~ Proba Array: ";
                        affiche_proba_array(m_worst_individual);
                    }
                }

                // affichage populationSize/2
                {
                    std::cout << "\n\tHalf fitness: " << m_population[m_config.get_half_population_size()].get_fitness();
                    if (m_config.enable_auto_adaptation) {
                        std::cout << " ~ Proba Array: ";
                        affiche_proba_array(m_population[m_config.get_half_population_size()]);
                    }
                }

                // affichage random
                {
                    size_t indice = m_index_dist(m_rng);
                    std::cout << "\n\tRandom fitness: " << m_population[indice].get_fitness();
                    if (m_config.enable_auto_adaptation) {
                        std::cout << " ~ Proba Array: ";
                        affiche_proba_array(m_population[indice]);
                    }
                }
                std::cout << '\n';
            }

            if (callback) {
                callback(gen, m_best_fitness, m_best_individual);
            }
        }

        if (verbose) {
            std::cout << "\nFinal best fitness: " << m_best_fitness << std::endl;
            std::cout << "Best individual:\n" << m_best_individual << std::endl;
        }
    }

    void reset(const ConfigType& config) {
        config.validate();

        m_config = config;
        m_best_fitness = std::numeric_limits<Real>::lowest();
        m_worst_fitness = m_best_fitness;
        m_index_dist = std::uniform_int_distribution<size_t>(0, config.population_size - 1);
        m_cut_dist = std::uniform_int_distribution<size_t>(0, (config.dimension * config.integer_bits) - 1);
        m_cut_dist_proba =
            std::uniform_int_distribution<size_t>(0, ((config.number_of_vectors + 1) * config.integer_bits) - 1);

        m_best_individual = Individual(m_config, m_rng);
        initialize_population();

        m_selected.clear();
        m_selected.resize(m_config.get_half_population_size());
    }

    // ========== Getters ==========

    [[nodiscard]] Real get_best_fitness() const { return m_best_fitness; }
    [[nodiscard]] Real get_worst_fitness() const { return m_worst_fitness; }
    [[nodiscard]] const Individual& get_best_individual() const { return m_best_individual; }
    [[nodiscard]] const Individual& get_worst_individual() const { return m_worst_individual; }
    [[nodiscard]] const Population& get_population() const { return m_population; }
    [[nodiscard]] const Config<Real, Integer>& get_config() const { return m_config; }

private:
    /**
     * @brief Initializes the population with random individuals
     */
    void initialize_population() {
        m_population.clear();
        m_population.reserve(m_config.population_size);

        for (size_t i = 0; i < m_config.population_size; i++) {
            m_population.emplace_back(m_config, m_rng);
        }

        update_best();
    }

    void evaluate_population() {
        #pragma omp parallel
        {
            #pragma omp for
            for (int i = 0; i < m_config.population_size; i++) {
                auto real_vecs = m_population[i].to_real_vectors();
                Real eval = m_fitness_func(real_vecs);
                m_population[i].set_fitness(eval);

                /*#pragma omp critical
                {
                    std::cout << "Individu " << i << " fitness: " << eval << '\n';
                }*/
            }
        }
    }


    /**
     * @brief Evaluates the fitness of an individual, Not thread-safe
     */
    Real evaluate(Individual& ind) {
        if (ind.have_been_evaluated()) {
            return ind.get_fitness();
        }

        auto real_vecs = ind.to_real_vectors();
        Real eval = m_fitness_func(real_vecs);
        ind.set_fitness(eval);
        return eval;
    }

    /**
     * @brief Tournament selection
     */
    void selection(int gen) {
        const int half_pop = static_cast<int>(m_config.get_half_population_size());
        const size_t tournament_size = m_config.tournament_size;
        const size_t pop_size = m_population.size();

        #pragma omp parallel
        {
            // RNG thread-local
            std::mt19937 rng(m_seed + omp_get_thread_num() + gen);
            std::uniform_int_distribution<size_t> index_dist(0, pop_size - 1);

            // Chaque thread gère plusieurs tournois
            #pragma omp for
            for (int i = 0; i < half_pop; i++) {
                size_t best_idx = index_dist(rng); // tirage initial
                Real best_eval = m_population[best_idx].get_fitness();

                // Tirage des k-1 autres concurrents du tournoi
                for (size_t t = 1; t < tournament_size; t++) {
                    size_t idx = index_dist(rng);
                    Real eval = m_population[idx].get_fitness();
                    if (eval > best_eval) {
                        best_eval = eval;
                        best_idx = idx;
                    }
                }

                // Stockage du gagnant
                m_selected[i] = m_population[best_idx];
            }
        }
    }


    /**
     * @brief Crossover between two parents at bit level
     */
    std::pair<Individual, Individual> crossover(const Individual& p1, const Individual& p2, std::mt19937& rng) {
        Individual child1 = p1;
        Individual child2 = p2;

        const size_t num_vecs = m_config.number_of_vectors;

        // Crossover des chromosomes de données
        for (size_t chromo = 0; chromo < num_vecs; chromo++) {
            bit_level_crossover(p1, p2, child1, child2, chromo, rng);
        }

        // Crossover des probabilités de mutation
        if (m_config.enable_auto_adaptation) {
            bit_level_crossover_probas(p1, p2, child1, child2, rng);
        }

        child1.invalidate_fitness();
        child2.invalidate_fitness();

        return {child1, child2};
    }

    /**
     * @brief Bit-level crossover for a chromosome
     */
    void bit_level_crossover(const Individual& p1, const Individual& p2, Individual& child1, Individual& child2,
                             size_t chromo, std::mt19937& rng) {
        switch (m_config.crossover_method) {
            case CrossoverType::SINGLE_POINT_BIT_LEVEL:
                bit_level_crossover_single_point_bit_level(p1, p2, child1, child2, chromo, rng);
                break;
            case CrossoverType::UNIFORM_BIT_LEVEL:
                bit_level_crossover_uniform_bit_level(p1, p2, child1, child2, chromo, rng);
                break;
        }
    }

    /**
     * @brief Bit-level crossover for a chromosome using single point bit level
     */
    void bit_level_crossover_single_point_bit_level(const Individual& p1, const Individual& p2, Individual& child1,
                                                    Individual& child2, size_t chromo, std::mt19937& rng) {
        size_t cut_point = m_cut_dist(rng);

        size_t k = cut_point / m_config.integer_bits;
        size_t k_prime = cut_point % m_config.integer_bits;

        // Copier avant le cut
        for (size_t i = 0; i < k; i++) {
            child1.set_gene(chromo, i, p1.get_gene(chromo, i));
            child2.set_gene(chromo, i, p2.get_gene(chromo, i));
        }

        // Gérer le gène coupé
        if (k < m_config.dimension) {
            Integer g1 = p1.get_gene(chromo, k);
            Integer g2 = p2.get_gene(chromo, k);

            Integer mask1 = (Integer(1) << k_prime) - Integer(1);
            Integer mask2 = ~mask1;

            child1.set_gene(chromo, k, (g1 & mask1) | (g2 & mask2));
            child2.set_gene(chromo, k, (g2 & mask1) | (g1 & mask2));
        }

        // Copier après le cut (de l'autre parent)
        for (size_t i = k + 1; i < m_config.dimension; i++) {
            child1.set_gene(chromo, i, p2.get_gene(chromo, i));
            child2.set_gene(chromo, i, p1.get_gene(chromo, i));
        }
    }

    /**
     * @brief Bit-level crossover for a chromosome using uniform bit level
     */
    void bit_level_crossover_uniform_bit_level(const Individual& p1, const Individual& p2, Individual& child1,
                                               Individual& child2, size_t chromo, std::mt19937& rng) {
        // on itère sur chaque gène
        for (size_t k = 0; k < m_config.dimension; k++) {
            Integer g1 = p1.get_gene(chromo, k);
            Integer g2 = p2.get_gene(chromo, k);

            // on itère sur chaque bit
            Integer gene_child_1 = 0;
            Integer gene_child_2 = 0;
            for (size_t i = 0; i < m_config.integer_bits; i++) {
                Integer mask = Integer(1) << i;
                if (ConfigType::proba_distribution(rng) <= m_config.uniform_crossover_probability) {
                    // on met le bit du parent 1 pour l'enfant 1
                    gene_child_1 |= g1 & mask;
                    gene_child_2 |= g2 & mask;
                } else {
                    // on met le bit du parent 1 pour l'enfant 2
                    gene_child_1 |= g2 & mask;
                    gene_child_2 |= g1 & mask;
                }
            }

            child1.set_gene(chromo, k, gene_child_1);
            child2.set_gene(chromo, k, gene_child_2);
        }
    }

    /**
     * @brief Bit-level crossover for mutation probabilities
     */
    void bit_level_crossover_probas(const Individual& p1, const Individual& p2, Individual& child1,
                                    Individual& child2, std::mt19937& rng) {
        switch (m_config.crossover_method) {
            case CrossoverType::SINGLE_POINT_BIT_LEVEL:
                bit_level_crossover_single_point_bit_level_probas(p1, p2, child1, child2, rng);
                break;
            case CrossoverType::UNIFORM_BIT_LEVEL:
                bit_level_crossover_uniform_bit_level_probas(p1, p2, child1, child2, rng);
                break;
        }
    }

    /**
     * @brief Bit-level crossover for mutation probabilities using single point bit level
     */
    void bit_level_crossover_single_point_bit_level_probas(const Individual& p1, const Individual& p2,
                                                           Individual& child1, Individual& child2, std::mt19937& rng) {
        size_t cut_point = m_cut_dist_proba(rng);

        size_t k = cut_point / m_config.integer_bits;
        size_t k_prime = cut_point % m_config.integer_bits;

        for (size_t i = 0; i < k; i++) {
            child1.set_mutation_proba(i, p1.get_mutation_proba(i));
            child2.set_mutation_proba(i, p2.get_mutation_proba(i));
        }

        if (k < m_config.number_of_vectors + 1) {
            Integer pr1 = p1.get_mutation_proba(k);
            Integer pr2 = p2.get_mutation_proba(k);

            Integer mask1 = (Integer(1) << k_prime) - Integer(1);
            Integer mask2 = ~mask1;

            child1.set_mutation_proba(k, (pr1 & mask1) | (pr2 & mask2));
            child2.set_mutation_proba(k, (pr2 & mask1) | (pr1 & mask2));
        }

        for (size_t i = k + 1; i < m_config.number_of_vectors + 1; i++) {
            child1.set_mutation_proba(i, p2.get_mutation_proba(i));
            child2.set_mutation_proba(i, p1.get_mutation_proba(i));
        }
    }

    /**
     * @brief Bit-level crossover for mutation probabilities using uniform bit level
     */
    void bit_level_crossover_uniform_bit_level_probas(const Individual& p1, const Individual& p2, Individual& child1,
                                                      Individual& child2, std::mt19937& rng) {
        // on itère sur chaque probabilités
        for (size_t k = 0; k < m_config.number_of_vectors + 1; k++) {
            Integer g1 = p1.get_mutation_proba(k);
            Integer g2 = p2.get_mutation_proba(k);

            // on itère sur chaque bit
            Integer gene_child_1 = 0;
            Integer gene_child_2 = 0;
            for (size_t i = 0; i < m_config.integer_bits; i++) {
                Integer mask = Integer(1) << i;
                if (ConfigType::proba_distribution(rng) <= m_config.uniform_crossover_probability) {
                    // on met le bit du parent 1 pour l'enfant 1
                    gene_child_1 |= g1 & mask;
                    gene_child_2 |= g2 & mask;
                } else {
                    // on met le bit du parent 1 pour l'enfant 2
                    gene_child_1 |= g2 & mask;
                    gene_child_2 |= g1 & mask;
                }
            }

            child1.set_mutation_proba(k, gene_child_1);
            child2.set_mutation_proba(k, gene_child_2);
        }
    }

    /**
     * @brief Creates offspring from selected parents
     */
    void create_offspring(int gen) {
        // Mélanger les parents
        const size_t half = m_selected.size();
        std::shuffle(m_selected.begin(), m_selected.end(), m_rng);

        // Créer des paires et faire crossover
        #pragma omp parallel
        {
            std::mt19937 rng(m_seed + omp_get_thread_num() + gen);

            #pragma omp for
            for (int i = 0; i < half; i += 2) {
                // Ajouter les enfants du crossover
                auto [child1, child2] = crossover(m_selected[i], m_selected[i + 1], rng);
                m_population[i] = std::move(child1);
                m_population[i + 1] = std::move(child2);
            }
        }

        // on le fait une deuxième fois pour avoir un population complète :
        // Mélanger les parents
        std::shuffle(m_selected.begin(), m_selected.end(), m_rng);

        #pragma omp parallel
        {
            std::mt19937 rng(m_seed + omp_get_thread_num() + gen + 123456);

            // Créer des paires et faire crossover
            #pragma omp for
            for (int i = 0; i < half; i += 2) {
                // Ajouter les enfants du crossover
                auto [child1, child2] = crossover(m_selected[i], m_selected[i + 1], rng);
                m_population[half + i] = std::move(child1);
                m_population[half + i + 1] = std::move(child2);
            }
        }
    }

    /**
     * @brief Applies mutation to entire population
     */
    void mutate_population(int gen) {
        #pragma omp parallel
        {
            std::mt19937 rng(m_seed + omp_get_thread_num() + gen);

            #pragma omp for
            for (int i = 0; i < m_population.size(); i++) {
                m_population[i].mutate(rng);
            }
        }
    }

    /**
     * @brief Runs one generation
     */
    void step(int gen) {
		//std::cout << "Generation " << gen << ":\n";
        // évaluation
        evaluate_population(); // multithreadé

        // Sélection
        selection(gen); // multithreadé

        // sauvegarder le meilleur pour plus tard si il y a élitisme
        if (m_config.enable_elitism) {
            update_best(); // multithreadé
        }

        // Crossover
        create_offspring(gen);

        // Mutation
        mutate_population(gen); // multithreadé (peut-être pas entièrement, i.e. pas au niveau de chaque individu)

        // mettre le meilleur quelque part dans la nouvelle population si il y a élitisme
        if (m_config.enable_elitism) {
            add_best();
        }

		//std::cout << "\nAfter mutation of generation " << gen << ":\n";
        evaluate_population(); // les mutés n'ont pas de score valide
        // Mise à jour du meilleur
        update_best(); // multithreadé // update aussi le pire

        //std::cout << "\n\n\n\n";
    }

    /**
     * @brief Updates the best individual found so far
     * @pre Selection function need to have been executed
     */
    void update_best() {
        Real best_fitness_local = std::numeric_limits<Real>::lowest();
        Real worst_fitness_local = std::numeric_limits<Real>::max();
        Individual best_ind_local;
        Individual worst_ind_local;

        #pragma omp parallel
        {
            Real thread_best = std::numeric_limits<Real>::lowest();
            Individual thread_best_ind;

            Real thread_worst = std::numeric_limits<Real>::lowest();
            Individual thread_worst_ind;

            // Parcours parallèle
            #pragma omp for nowait
            for (int i = 0; i < static_cast<int>(m_population.size()); ++i) {
                Real f = m_population[i].get_fitness();
                if (f > thread_best) {
                    thread_best = f;
                    thread_best_ind = m_population[i];
                }
                else if (f < thread_worst) {
                    thread_worst = f;
                    thread_worst_ind = m_population[i];
                }
            }

            // Mise à jour globale sécurisée
            #pragma omp critical
            {
                if (thread_best > best_fitness_local) {
                    best_fitness_local = thread_best;
                    best_ind_local = thread_best_ind;
                }
                
                if (thread_worst < worst_fitness_local) {
                    worst_fitness_local = thread_worst;
                    worst_ind_local = thread_worst_ind;
                }
            }
        }

        // Mise à jour finale
        m_best_fitness = best_fitness_local;
        m_best_individual = best_ind_local;

        m_worst_fitness = worst_fitness_local;
        m_worst_individual = worst_ind_local;
    }


    void add_best() {
        // on sélectionne aléatoirement un individu, qu'on remplace par le meilleur actuel pré-sauvegardé
        size_t indice = m_index_dist(m_rng);
        m_population[indice] = m_best_individual;
    }
};
}  // namespace genetic
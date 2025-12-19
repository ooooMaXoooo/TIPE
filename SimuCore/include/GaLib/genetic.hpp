/**
 * @file genetic.hpp
 * @brief Main header file for the genetic algorithm library
 *
 * This header-only library provides a flexible and generic implementation
 * of genetic algorithms with the following features:
 *
 * - Template-based design for custom Real and Integer types
 * - Configurable population size, mutation rates, and other parameters
 * - Custom fitness functions via std::function or lambdas
 * - Self-adaptive mutation probabilities
 * - Tournament selection
 * - Bit-level crossover
 *
 * @example example_simple.cpp
 * @example example_Rosenbrock.cpp
 *
 * @version 1.0
 * @author Lemoine Maxence
 */

#pragma once

// Include all components
#include "genetic_algorithm.hpp"  // IWYU pragma: export
#include "genetic_config.hpp"     // IWYU pragma: export
#include "genetic_individu.hpp"   // IWYU pragma: export
#include "genetic_utils.hpp"      // IWYU pragma: export

namespace genetic {

/**
 * @brief Library version information
 */
constexpr const char* VERSION = "1.0.0";

/**
 * @brief Prints library information
 */
inline void print_info() {
    std::cout << "Genetic Algorithm Library v" << VERSION << '\n';
    std::cout << "A flexible header-only library for genetic algorithms" << '\n';
}

}  // namespace genetic
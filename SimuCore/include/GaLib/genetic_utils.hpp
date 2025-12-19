#pragma once

#include <limits>
#include <type_traits>

namespace genetic::utils {

/**
 * @brief Converts a binary integer to a real number in [min_real, max_real]
 *
 * @tparam Real Floating-point type
 * @tparam Integer Unsigned integer type
 * @param bin_value Binary encoded value
 * @param min_real Minimum real value
 * @param max_real Maximum real value
 * @param num_bits Number of bits used for encoding
 * @return Real value in [min_real, max_real]
 */
template <typename Real, typename Integer>
inline Real bin_to_real(Integer bin_value, Real min_real, Real max_real, size_t num_bits) {
    static_assert(std::is_floating_point_v<Real>, "Real must be floating-point");
    static_assert(std::is_unsigned_v<Integer>, "Integer must be unsigned");

    // Calculer la valeur maximale pour num_bits
    Integer max_int_value;
    if (num_bits == std::numeric_limits<Integer>::digits) {
        max_int_value = std::numeric_limits<Integer>::max();
    } else {
        max_int_value = (Integer(1) << num_bits) - Integer(1);
    }

    // Normaliser entre 0 et 1
    Real normalized = static_cast<Real>(bin_value) / static_cast<Real>(max_int_value);

    // Mapper sur [min_real, max_real]
    return min_real + (normalized * (max_real - min_real));
}

/**
 * @brief Converts a real number to binary representation
 *
 * @tparam Real Floating-point type
 * @tparam Integer Unsigned integer type
 * @param real_value Real value in [min_real, max_real]
 * @param min_real Minimum real value
 * @param max_real Maximum real value
 * @param num_bits Number of bits used for encoding
 * @return Binary encoded value
 */
template <typename Real, typename Integer>
inline Integer real_to_bin(Real real_value, Real min_real, Real max_real, size_t num_bits) {
    static_assert(std::is_floating_point_v<Real>, "Real must be floating-point");
    static_assert(std::is_unsigned_v<Integer>, "Integer must be unsigned");

    // Calculer la valeur maximale pour num_bits
    Integer max_int_value;
    if (num_bits == std::numeric_limits<Integer>::digits) {
        max_int_value = std::numeric_limits<Integer>::max();
    } else {
        max_int_value = (Integer(1) << num_bits) - Integer(1);
    }

    // Clamp la valeur dans [min_real, max_real]
    if (real_value < min_real) {
        real_value = min_real;
    }
    if (real_value > max_real) {
        real_value = max_real;
    }

    // Normaliser entre 0 et 1
    Real normalized = (real_value - min_real) / (max_real - min_real);

    // Convertir en entier
    return static_cast<Integer>(normalized * static_cast<Real>(max_int_value));
}

/**
 * @brief Converts a probability [0, 1] to binary representation
 *
 * @tparam Real Floating-point type
 * @tparam Integer Unsigned integer type
 * @param proba Probability in [0, 1]
 * @param num_bits Number of bits used for encoding
 * @return Binary encoded probability
 */
template <typename Real, typename Integer>
inline Integer proba_to_bin(Real proba, size_t num_bits) {
    static_assert(std::is_floating_point_v<Real>, "Real must be floating-point");
    static_assert(std::is_unsigned_v<Integer>, "Integer must be unsigned");

    // Clamp entre 0 et 1
    if (proba < Real(0)) {
        proba = Real(0);
    }
    if (proba > Real(1)) {
        proba = Real(1);
    }

    // Calculer la valeur maximale pour num_bits
    Integer max_int_value;
    if (num_bits == std::numeric_limits<Integer>::digits) {
        max_int_value = std::numeric_limits<Integer>::max();
    } else {
        max_int_value = (Integer(1) << num_bits) - Integer(1);
    }

    return static_cast<Integer>(proba * static_cast<Real>(max_int_value));
}

/**
 * @brief Converts binary representation to probability [0, 1]
 *
 * @tparam Real Floating-point type
 * @tparam Integer Unsigned integer type
 * @param bin_proba Binary encoded probability
 * @param num_bits Number of bits used for encoding
 * @return Probability in [0, 1]
 */
template <typename Real, typename Integer>
inline Real bin_to_proba(Integer bin_proba, size_t num_bits) {
    static_assert(std::is_floating_point_v<Real>, "Real must be floating-point");
    static_assert(std::is_unsigned_v<Integer>, "Integer must be unsigned");

    // Calculer la valeur maximale pour num_bits
    Integer max_int_value;
    if (num_bits == std::numeric_limits<Integer>::digits) {
        max_int_value = std::numeric_limits<Integer>::max();
    } else {
        max_int_value = (Integer(1) << num_bits) - Integer(1);
    }

    return static_cast<Real>(bin_proba) / static_cast<Real>(max_int_value);
}

}  // namespace genetic::utils

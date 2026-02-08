#pragma once

constexpr inline double daysInSeconds(double days) noexcept {
	return 24 * 60 * 60 * days;
}

/// <summary>
/// Convertit linéairement une valeur depuis l'intervalle [min_1, max_1] vers l'intervalle [min_2, max_2].
/// </summary>
/// <param name="min_1">Borne minimale de l'intervalle source.</param>
/// <param name="max_1">Borne maximale de l'intervalle source.</param>
/// <param name="min_2">Borne minimale de l'intervalle cible.</param>
/// <param name="max_2">Borne maximale de l'intervalle cible.</param>
/// <param name="value">Valeur à convertir depuis l'intervalle source.</param>
/// <returns>La valeur mappée dans l'intervalle cible obtenue par interpolation linéaire. Si max_1 == min_1, le comportement est indéfini (division par zéro).</returns>
constexpr inline double convertIntervals(double min_1, double max_1, double min_2, double max_2, double value) {
	// <proportion dans I1> * <Taille de I2> + min_I2
	return ((value-min_1) / (max_1 - min_1)) * (max_2 - min_2) + min_2;
}
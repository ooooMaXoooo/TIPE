#pragma once

constexpr inline double daysInSeconds(double days) noexcept {
	return 24 * 60 * 60 * days;
}

constexpr inline double convertIntervals(double min_1, double max_1, double min_2, double max_2, double value) {
	// <proportion dans I1> * <Taille de I2> + min_I2
	return ((value-min_1) / (max_1 - min_1)) * (max_2 - min_2) + min_2;
}

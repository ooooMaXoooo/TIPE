#pragma once
#include <SimuCore/constants.h>


constexpr long double operator""_s_to_h(long double t) noexcept {
	return t / 3600;
}

constexpr long double operator""_h_to_s(long double t) noexcept {
	return t * 3600;
}

constexpr long double operator""_s_to_d(long double t) noexcept {
	return t / 86400;
}

constexpr long double operator""_d_to_s(long double t) noexcept {
	return t * 86400;
}

constexpr long double operator""_h_to_d(long double t) noexcept {
	return t / 24;
}

constexpr long double operator""_d_to_h(long double t) noexcept {
	return t * 24;
}




constexpr inline double seconds_to_hours(double t_s) noexcept {
	return t_s / 3600;
}

constexpr inline double hours_to_seconds(double t_h) noexcept {
	return t_h * 3600;
}

constexpr inline double seconds_to_days(double t_s) noexcept {
	return t_s / 86400;
}

constexpr inline double days_to_seconds(double t_d) noexcept {
	return t_d * 86400;
}

constexpr inline double hours_to_days(double t_h) noexcept {
	return t_h / 24;
}

constexpr inline double days_to_hours(double t_d) noexcept {
	return t_d * 24;
}

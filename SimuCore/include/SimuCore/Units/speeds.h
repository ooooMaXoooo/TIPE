#pragma once
#include <SimuCore/constants.h>


constexpr long double operator""_m_per_s__to__km_per_h(long double v_m_per_s) {
	return v_m_per_s * 3.6;
}

constexpr long double operator""_km_per_h__to__m_per_s(long double v_km_per_h) {
	return v_km_per_h / 3.6;
}

constexpr inline double meters_per_seconds_to_kilometers_per_hours(double v_m_per_s) noexcept {
	return v_m_per_s * 3.6;
}

constexpr inline double kilometers_per_hours_to_meters_per_seconds(double v_km_per_h) noexcept {
	return v_km_per_h / 3.6;
}

constexpr long double operator""_m_per_s__to__km_per_s(long double v_m_per_s) {
	return v_m_per_s * 1e-3;
}

constexpr long double operator""_km_per_s__to__m_per_s(long double v_km_per_s) {
	return v_km_per_s * 1e+3;
}

constexpr inline double meters_per_seconds_to_kilometers_per_seconds(double v_m_per_s) noexcept {
	return v_m_per_s * 1e-3;
}

constexpr inline double kilometers_per_seconds_to_meters_per_seconds(double v_km_per_s) noexcept {
	return v_km_per_s * 1e+3;
}

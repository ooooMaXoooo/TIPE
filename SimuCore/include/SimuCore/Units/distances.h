#pragma once
#include <SimuCore/constants.h>


constexpr long double operator""_m_to_AU(long double d) noexcept {
	return d / SimuCore::constants::AU;
}

constexpr long double operator""_AU_to_m(long double d) noexcept {
	return d * SimuCore::constants::AU;
}

constexpr long double operator""_m_to_km(long double d) noexcept {
	return d * 1e-3;
}

constexpr long double operator""_km_to_m(long double d) noexcept {
	return d * 1e3;
}

constexpr long double operator""_km_to_AU(long double d) noexcept {
	return d / 1.496e+8;
}

constexpr long double operator""_AU_to_km(long double d) noexcept {
	return d * 1.496e+8;
}


constexpr inline double meters_to_AU(double d_meters) noexcept {
	return d_meters / SimuCore::constants::AU;
}

constexpr inline double AU_to_meters(double d_UA) noexcept {
	return d_UA * SimuCore::constants::AU;
}

constexpr inline double meters_to_kilometers(double d_meters) noexcept {
	return d_meters * 1e-3;
}

constexpr inline double kilometers_to_meters(double d_km) noexcept {
	return d_km * 1e+3;
}

constexpr inline double AU_to_kilometers(double d_AU) noexcept {
	return d_AU * 1.496e+8;
}

constexpr inline double kilometers_to_AU(double d_km) noexcept {
	return d_km / 1.496e+8;
}
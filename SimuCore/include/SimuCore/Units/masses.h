#pragma once

constexpr long double operator""_kg_to_ton(long double m_kg) noexcept {
	return m_kg * 1e-3;
}

constexpr inline double kilograms_to_tons(double m_kg) noexcept {
	return m_kg * 1e-3;
}

constexpr long double operator""_ton_to_kg(long double m_ton) noexcept {
	return m_ton * 1e+3;
}

constexpr inline double tons_to_kilograms(double m_ton) noexcept {
	return m_ton * 1e+3;
}
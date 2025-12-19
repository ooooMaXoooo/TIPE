#pragma once

namespace SimuCore {
	namespace constants {
		constexpr double mu = 1.32712440018e20; // Gravitation du Soleil
		constexpr double G = 6.67428e-11; // --> à passer en km
		constexpr double g = 9.80665; // m/s²

		constexpr double softening = 1e1; // mètre --> à passer en km

		constexpr double PI = 3.14159265358979323846;

		constexpr double AU = 149597870700.0; // mètre --> à passer en km

		constexpr double epsilon = 1e-12;
	}
};
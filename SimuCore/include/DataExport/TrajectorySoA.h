#pragma once
#include <pch.h>

namespace dataExport {

	struct TrajectorySoA { // SoA = Structure of Arrays  (par opposition à Array of Structures, AoS)
		std::vector<double> t;			// temps
		std::vector<double> px, py, pz; // positions
		std::vector<double> vx, vy, vz; // vitesses

		bool is_valid; // validité de la trajectoire
	}; // struct TrajectorySoa

	struct BestIndividualUpdate {
		bool is_new_best = false;

		// Trajectoire uniquement si is_new_best == true
		TrajectorySoA trajectory;
	};

}; // namespace dataExport
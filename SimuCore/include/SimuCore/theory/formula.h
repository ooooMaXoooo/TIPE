#pragma once

#include <pch.h>
#include <SimuCore/constants.h>
#include <SimuCore/structures/structures.h>

namespace stumpff {
    /* --- Stumpff robustes --- */
    double C(double x);
    double S(double x);
};

namespace orbit {
    /// <summary>
	/// Renvoie les vitesses initiale et finale pour une manœuvre de Lambert entre r1 et r2 en un temps tof. Etant donné le paramètre gravitationnel mu de l'astre attracteur.
    /// </summary>
    /// <param name="r1">vecteur position initial</param>
    /// <param name="r2">vecteur position final</param>
    /// <param name="tof">temps total du voyage</param>
    /// <param name="mu">paramètre gravitationnel</param>
    /// <returns></returns>
    std::pair<std::array<double, 3>, std::array<double, 3>> lambert_universal(
        const std::array<double, 3>& r1,
        const std::array<double, 3>& r2,
        double tof,
        double mu = SimuCore::constants::mu
    );

    std::vector<double> lambert_batch(
        const std::array<double, 3>& r1,
        const std::vector<std::array<double, 3>>& r2_list,
        const std::vector<double>& tof_list,
        double mu = SimuCore::constants::mu
    );
};

glm::dvec3 forceAttractionGrav(const SimuCore::Structures::Entity& from, const SimuCore::Structures::Entity& on);
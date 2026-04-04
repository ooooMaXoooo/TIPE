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


/// <summary>
/// 
/// </summary>
/// <param name="distance_to_central_body"> en m </param>
/// <param name="mass_of_system"> en kg </param>
/// <param name="mu_central_body"> en USI, i.e m^3 / s^2 </param>
/// <param name="constante_des_aires"> en USI, i.e m^2 / s </param>
/// <param name="system_energy"> en J, kg * m^2 / s^2 </param>
/// <param name="is_trajectory_elliptic"> un pointeur sur un booléen </param>
/// <returns></returns>
std::pair<double, double> calcul_perige_et_apoge(
    double distance_to_central_body,
    double mass_of_system,
    double mu_central_body,
    double constante_des_aires,
    double system_energy,
    bool* is_trajectory_elliptic
);
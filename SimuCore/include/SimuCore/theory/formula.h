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

    /// <summary>
    /// </summary>
    /// <returns>le delta_v minimal, puis les 10 scalaires qui définissent l'état du système pour l'otpimal</returns>
    std::array<double, 17> lambert_by_parts(
        double mu_start, double mu_final, double mu_sun,
        double distFinalSun,
        double tof1_min, double tof1_max,
        double tof2_min, double tof2_max,
        double tof3_min, double tof3_max,
        double r1_min, double r1_max, double theta1_min, double theta1_max,
        double r2_min, double r2_max, double theta2_min, double theta2_max,
        double RStart, double thetaStart_min, double thetaStart_max,
        double RFinal, double thetaFinal_min, double thetaFinal_max,
        double thetaFinalPlanet_min, double thetaFinalPlanet_max,
        double Nt1, double Nt2, double Nt3,
        double Nr1, double Nth1,
        double Nr2, double Nth2,
        double NthStart, double NthFinal,
        double NthFinalPlanet
    );

    /// <summary>
    /// </summary>
    /// <returns>le delta_v minimal, puis les 10 scalaires qui définissent l'état du système pour l'otpimal</returns>
    std::array<double, 11> lambert_by_parts_OpenMP(
        double mu_start, double mu_final, double mu_sun,
        double distFinalSun,
        double tof1_min, double tof1_max,
        double tof2_min, double tof2_max,
        double tof3_min, double tof3_max,
        double r1_min, double r1_max, double theta1_min, double theta1_max,
        double r2_min, double r2_max, double theta2_min, double theta2_max,
        double RStart, double thetaStart_min, double thetaStart_max,
        double RFinal, double thetaFinal_min, double thetaFinal_max,
        double thetaFinalPlanet_min, double thetaFinalPlanet_max,
        double Nt1, double Nt2, double Nt3,
        double Nr1, double Nth1,
        double Nr2, double Nth2,
        double NthStart, double NthFinal,
        double NthFinalPlanet
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
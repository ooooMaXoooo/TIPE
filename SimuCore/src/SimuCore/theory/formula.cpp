#include "pch.h"
#include "theory/formula.h"
#include "constants.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <Units/all.h>


double stumpff::C(double z) {
    if (z > 1e-8) return (1 - std::cos(std::sqrt(z))) / z;
    if (z < -1e-8) return (std::cosh(std::sqrt(-z)) - 1) / -z;
    return 0.5 - z / 24 + z * z / 720;
}

double stumpff::S(double z) {
    if (z > 1e-8) {
        double s = std::sqrt(z);
        return (s - std::sin(s)) / (s * s * s);
    }
    if (z < -1e-8) {
        double s = std::sqrt(-z);
        return (std::sinh(s) - s) / (s * s * s);
    }
    return 1.0 / 6 - z / 120 + z * z / 5040;
}

std::pair<std::array<double, 3>, std::array<double, 3>> orbit::lambert_universal(
    const std::array<double, 3>& r1,
    const std::array<double, 3>& r2,
    double tof,
    double mu
) {
    auto norm = [](const std::array<double, 3>& v) {
        return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    };

    auto dot = [](const std::array<double, 3>& a, const std::array<double, 3>& b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    };

    constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

    double r1_norm = norm(r1);
    double r2_norm = norm(r2);
    double cos_dtheta = std::max(-1.0, std::min(1.0, dot(r1, r2) / (r1_norm * r2_norm)));
    double dtheta = std::acos(cos_dtheta);

    if (std::abs(dtheta) < 1e-8) return { {NaN, NaN, NaN}, {NaN, NaN, NaN} };

    double A = std::sin(dtheta) * std::sqrt(r1_norm * r2_norm / (1 - cos_dtheta));
    if (A == 0) return { {NaN, NaN, NaN}, {NaN, NaN, NaN} };

    auto tof_of_z = [&](double z) -> double {
        double C = stumpff::C(z);
        double S = stumpff::S(z);
        if (C <= 0) return NaN;
        double y = r1_norm + r2_norm + A * (z * S - 1) / std::sqrt(C);
        if (y <= 0) return NaN;
        double chi = std::sqrt(y / C);
        return (std::pow(chi, 3) * S + A * std::sqrt(y)) / std::sqrt(mu);
        };

    double z = NaN;
    for (int i = 0; i < 1199; ++i) {
        double z1 = -4 * SimuCore::constants::PI * SimuCore::constants::PI + i * (8 * SimuCore::constants::PI * SimuCore::constants::PI / 1200);
        double z2 = z1 + (8 * SimuCore::constants::PI * SimuCore::constants::PI / 1200);
        double f1 = tof_of_z(z1) - tof;
        double f2 = tof_of_z(z2) - tof;
        if (std::isfinite(f1) && std::isfinite(f2) && f1 * f2 < 0) {
            z = (z1 + z2) / 2;
            break;
        }
    }
    if (!std::isfinite(z)) return { {NaN, NaN, NaN}, {NaN, NaN, NaN} };

    double C = stumpff::C(z);
    double S = stumpff::S(z);
    double y = r1_norm + r2_norm + A * (z * S - 1) / std::sqrt(C);
    double chi = std::sqrt(y / C);
    double f = 1 - y / r1_norm;
    double g = A * std::sqrt(y / mu);
    double gdot = 1 - y / r2_norm;

    if (std::abs(g) < 1e-6) return { {NaN, NaN, NaN}, {NaN, NaN, NaN} };

    std::array<double, 3> v1, v2;
    for (int i = 0; i < 3; ++i) {
        v1[i] = (r2[i] - f * r1[i]) / g;
        v2[i] = (gdot * r2[i] - r1[i]) / g;
    }

    return { v1, v2 };
}


std::vector<double> orbit::lambert_batch(
    const std::array<double, 3>& r1,
    const std::vector<std::array<double, 3>>& r2_list,
    const std::vector<double>& tof_list,
    double mu
) {
    size_t N = r2_list.size();
    if (N != tof_list.size()) {
        throw std::invalid_argument("r2_list and tof_list must have same length");
    }

    std::vector<double> deltaV(N, std::numeric_limits<double>::quiet_NaN());

    // Optionnel : fixer le nombre de threads (ou garder OpenMP auto)
    // omp_set_num_threads(8);

#pragma omp parallel for
    for (int idx = 0; idx < static_cast<int>(N); ++idx) {
        try {
            const auto& r2 = r2_list[idx];
            double tof = tof_list[idx];

            auto res = orbit::lambert_universal(r1, r2, tof, mu);
            const auto& v1 = res.first;
            const auto& v2 = res.second;

            // Vérifier la validité
            bool ok = true;
            for (int k = 0; k < 3; ++k) {
                if (!std::isfinite(v1[k]) || !std::isfinite(v2[k])) { ok = false; break; }
            }
            if (!ok) {
                deltaV[idx] = std::numeric_limits<double>::quiet_NaN();
                continue;
            }

            auto norme = [](const std::array<double, 3>& u) -> double {
                return std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
            };

            auto substract = [](const std::array<double, 3>& u, const std::array<double, 3>& v) {
                return std::array<double, 3>{u[0] - v[0], u[1] - v[1], u[2] - v[2] };
            };

            auto calcul_vecteur_cercle = [norme, mu](const std::array<double, 3>& r1) {
                std::array<double, 3> v = { -r1[1], r1[0], r1[2] };
                auto norme_voulu = std::sqrt(mu / norme(r1)) / norme(r1);
                for (int i = 0; i < 3; i++) { v[i] *= norme_voulu; }
                return v;
            };
            
            std::array<double, 3> vecteur_vitesse_orbite_depart = calcul_vecteur_cercle(r1);
            std::array<double, 3> vecteur_vitesse_orbite_fin = calcul_vecteur_cercle(r2);


            deltaV[idx] = norme(substract(v1, vecteur_vitesse_orbite_depart)) + norme(substract(vecteur_vitesse_orbite_fin, v2));
        }
        catch (...) {
            // en cas d'exception, garder NaN
            deltaV[idx] = std::numeric_limits<double>::quiet_NaN();
        }
    }

    return deltaV;
}

std::array<double, 17> orbit::lambert_by_parts(
    double mu_start, double mu_final, double mu_sun,
    double distFinalSun,
    double tof1_min, double tof1_max,
    double tof2_min, double tof2_max,
    double tof3_min, double tof3_max,
    double r1_min, double r1_max,
    double phi1_min, double phi1_max,
    double r2_min, double r2_max,
    double phi2_min, double phi2_max,
    double RStart,
    double phiStart_min, double phiStart_max,
    double RFinal,
    double phiFinal_min, double phiFinal_max,
    double phiFinalPlanet_min, double phiFinalPlanet_max,
    double Nt1, double Nt2, double Nt3,
    double Nr1, double Nphi1,
    double Nr2, double Nphi2,
    double NphiStart, double NphiFinal, double NphiFinalPlanet) {

    std::array<double, 17> res{};

    double delta_v_best = std::numeric_limits<double>::infinity();
    double phiFinalPlanet_best = phiFinalPlanet_min;
    double tof1_best = tof1_min;
    double tof2_best = tof2_min;
    double tof3_best = tof3_min;
    double phiStart_best = phiStart_min;
    double phiFinal_best = phiFinal_min;
    double r1_best = r1_min;
    double phi1_best = phi1_min;
    double r2_best = r2_min;
    double phi2_best = phi2_min;

    // Variables pour l'export Leapfrog Python
    double v1_start_x_best = 0.0;
    double v1_start_y_best = 0.0;
    double v2_start_x_best = 0.0;
    double v2_start_y_best = 0.0;
    double v3_start_x_best = 0.0;
    double v3_start_y_best = 0.0;

    {
        constexpr double pos_terre_x = AU_to_meters(1.);
        constexpr double pos_terre_y = 0;

        double d1 = (phiFinalPlanet_max - phiFinalPlanet_min) / NphiFinalPlanet;
        double d2 = (tof1_max - tof1_min) / Nt1;
        double d3 = (tof2_max - tof2_min) / Nt2;
        double d4 = (tof3_max - tof3_min) / Nt3;
        double d5 = (phiStart_max - phiStart_min) / NphiStart;
        double d6 = (phiFinal_max - phiFinal_min) / NphiFinal;
        double d7 = (r1_max - r1_min) / Nr1;
        double d8 = (phi1_max - phi1_min) / Nphi1;
        double d9 = (r2_max - r2_min) / Nr2;
        double d10 = (phi2_max - phi2_min) / Nphi2;

        for (double phiFinalPlanet = phiFinalPlanet_min; phiFinalPlanet < phiFinalPlanet_max; phiFinalPlanet += d1) {
            std::array<double, 3> pos_finalPlanet;
            {
                double cos_phi = std::cos(phiFinalPlanet);
                double sin_phi = std::sin(phiFinalPlanet);

                pos_finalPlanet[0] = cos_phi * distFinalSun;
                pos_finalPlanet[1] = sin_phi * distFinalSun;
                pos_finalPlanet[2] = 0;
            }

            for (double tof1 = tof1_min; tof1 < tof1_max; tof1 += d2)
                for (double tof2 = tof2_min; tof2 < tof2_max; tof2 += d3)
                    for (double tof3 = tof3_min; tof3 < tof3_max; tof3 += d4)
                        for (double phiStart = phiStart_min; phiStart < phiStart_max; phiStart += d5)
                        {
                            // CORRECTION : Coordonnées LOCALES relatives à la Terre
                            std::array<double, 3> p_intermediate_start_local;
                            p_intermediate_start_local[0] = std::cos(phiStart) * RStart;
                            p_intermediate_start_local[1] = std::sin(phiStart) * RStart;
                            p_intermediate_start_local[2] = 0;

                            // Coordonnées absolues héliocentriques pour l'Arc 2
                            std::array<double, 3> p_intermediate_start_helio;
                            p_intermediate_start_helio[0] = p_intermediate_start_local[0] + pos_terre_x;
                            p_intermediate_start_helio[1] = p_intermediate_start_local[1] + pos_terre_y;
                            p_intermediate_start_helio[2] = 0;

                            for (double r1 = r1_min; r1 < r1_max; r1 += d7)
                                for (double phi1 = phi1_min; phi1 < phi1_max; phi1 += d8)
                                {
                                    // CORRECTION : Coordonnées LOCALES relatives à la Terre
                                    std::array<double, 3> p_initial_local;
                                    p_initial_local[0] = std::cos(phi1) * r1;
                                    p_initial_local[1] = std::sin(phi1) * r1;
                                    p_initial_local[2] = 0;

                                    // Lambert autour de la Terre utilise les vecteurs locaux
                                    auto [traj1_start, traj1_final] = orbit::lambert_universal(p_initial_local, p_intermediate_start_local, tof1, mu_start);

                                    std::array<double, 3> v_circ_init;
                                    {
                                        double speed = std::sqrt(mu_start / r1);
                                        glm::dvec2 u_r = glm::normalize(glm::dvec2(p_initial_local[0], p_initial_local[1]));
                                        v_circ_init[0] = -speed * u_r.y;
                                        v_circ_init[1] = speed * u_r.x;
                                        v_circ_init[2] = 0;
                                    }

                                    double delta_v_start = std::sqrt(
                                        std::pow(traj1_start[0] - v_circ_init[0], 2) +
                                        std::pow(traj1_start[1] - v_circ_init[1], 2) +
                                        std::pow(traj1_start[2] - v_circ_init[2], 2)
                                    );

                                    for (double phiFinal = phiFinal_min; phiFinal < phiFinal_max; phiFinal += d6)
                                    {
                                        // CORRECTION : Coordonnées LOCALES relatives à la planète d'arrivée
                                        std::array<double, 3> p_intermediate_final_local;
                                        p_intermediate_final_local[0] = std::cos(phiFinal) * RFinal;
                                        p_intermediate_final_local[1] = std::sin(phiFinal) * RFinal;
                                        p_intermediate_final_local[2] = 0;

                                        // Coordonnées absolues héliocentriques pour l'Arc 2
                                        std::array<double, 3> p_intermediate_final_helio;
                                        p_intermediate_final_helio[0] = p_intermediate_final_local[0] + pos_finalPlanet[0];
                                        p_intermediate_final_helio[1] = p_intermediate_final_local[1] + pos_finalPlanet[1];
                                        p_intermediate_final_helio[2] = 0;

                                        for (double r2 = r2_min; r2 < r2_max; r2 += d9)
                                            for (double phi2 = phi2_min; phi2 < phi2_max; phi2 += d10)
                                            {
                                                // CORRECTION : Coordonnées LOCALES relatives à la planète d'arrivée
                                                std::array<double, 3> p_final_local;
                                                p_final_local[0] = std::cos(phi2) * r2;
                                                p_final_local[1] = std::sin(phi2) * r2;
                                                p_final_local[2] = 0;

                                                // Arc 2 : Héliocentrique utilise les vecteurs absolus
                                                auto [traj2_start, traj2_final] = lambert_universal(
                                                    p_intermediate_start_helio, p_intermediate_final_helio, tof2, mu_sun);

                                                double v_earth_y = std::sqrt(mu_sun / pos_terre_x);

                                                double delta_v_soi_earth = std::sqrt(
                                                    std::pow(traj2_start[0] - (traj1_final[0] + 0.0), 2) +
                                                    std::pow(traj2_start[1] - (traj1_final[1] + v_earth_y), 2) +
                                                    std::pow(traj2_start[2] - traj1_final[2], 2)
                                                );

                                                double v_fp_speed = std::sqrt(mu_sun / distFinalSun);
                                                double v_fp_x = -std::sin(phiFinalPlanet) * v_fp_speed;
                                                double v_fp_y = std::cos(phiFinalPlanet) * v_fp_speed;

                                                // Arc 3 : Lambert autour de la planète finale utilise les vecteurs locaux
                                                auto [traj3_start, traj3_final] = lambert_universal(
                                                    p_intermediate_final_local, p_final_local, tof3, mu_final);

                                                double delta_v_soi_final = std::sqrt(
                                                    std::pow(traj2_final[0] - (traj3_start[0] + v_fp_x), 2) +
                                                    std::pow(traj2_final[1] - (traj3_start[1] + v_fp_y), 2) +
                                                    std::pow(traj2_final[2] - traj3_start[2], 2)
                                                );

                                                double v_circ_final_speed = std::sqrt(mu_final / r2);
                                                glm::dvec2 u_r_final = glm::normalize(glm::dvec2(p_final_local[0], p_final_local[1]));
                                                std::array<double, 3> v_circ_final = { -u_r_final.y * v_circ_final_speed, u_r_final.x * v_circ_final_speed, 0.0 };

                                                double delta_v_arrival = std::sqrt(
                                                    std::pow(traj3_final[0] - v_circ_final[0], 2) +
                                                    std::pow(traj3_final[1] - v_circ_final[1], 2) +
                                                    std::pow(traj3_final[2] - v_circ_final[2], 2)
                                                );

                                                double delta_v_total = delta_v_start + delta_v_soi_earth + delta_v_soi_final + delta_v_arrival;

                                                if (delta_v_total < delta_v_best) {
                                                    delta_v_best = delta_v_total;
                                                    phiFinalPlanet_best = phiFinalPlanet;
                                                    tof1_best = tof1;
                                                    tof2_best = tof2;
                                                    tof3_best = tof3;
                                                    phiStart_best = phiStart;
                                                    phiFinal_best = phiFinal;
                                                    r1_best = r1;
                                                    phi1_best = phi1;
                                                    r2_best = r2;
                                                    phi2_best = phi2;

                                                    // Sauvegarde des vitesses d'impulsion pour Python
                                                    v1_start_x_best = traj1_start[0];
                                                    v1_start_y_best = traj1_start[1];
                                                    v2_start_x_best = traj2_start[0];
                                                    v2_start_y_best = traj2_start[1];
                                                    v3_start_x_best = traj3_start[0];
                                                    v3_start_y_best = traj3_start[1];
                                                }
                                            }
                                    }
                                }
                        }
        }
    }

    res[0] = delta_v_best;
    res[1] = phiFinalPlanet_best;
    res[2] = tof1_best;
    res[3] = tof2_best;
    res[4] = tof3_best;
    res[5] = phiStart_best;
    res[6] = phiFinal_best;
    res[7] = r1_best;
    res[8] = phi1_best;
    res[9] = r2_best;
    res[10] = phi2_best;

    // Ajouts des impulsions héliocentriques
    res[11] = v1_start_x_best;
    res[12] = v1_start_y_best;
    res[13] = v2_start_x_best;
    res[14] = v2_start_y_best;
    res[15] = v3_start_x_best;
    res[16] = v3_start_y_best;

    return res;
}

std::array<double, 11> orbit::lambert_by_parts_OpenMP(
    double mu_start, double mu_final, double mu_sun,
    double distFinalSun,
    double tof1_min, double tof1_max,
    double tof2_min, double tof2_max,
    double tof3_min, double tof3_max,
    double r1_min, double r1_max,
    double theta1_min, double theta1_max,
    double r2_min, double r2_max,
    double theta2_min, double theta2_max,
    double RStart,
    double thetaStart_min, double thetaStart_max,
    double RFinal,
    double thetaFinal_min, double thetaFinal_max,
    double thetaFinalPlanet_min, double thetaFinalPlanet_max,
    double Nt1, double Nt2, double Nt3,
    double Nr1,
    double Nth1,
    double Nr2,
    double Nth2,
    double NthStart,
    double NthFinal,
    double NthFinalPlanet) {

    std::array<double, 11> res{};

    // Variables globales qui stockeront le meilleur résultat absolu
    double delta_v_best_global = std::numeric_limits<double>::infinity();
    double thetaFinalPlanet_best_global = thetaFinalPlanet_min;
    double tof1_best_global = tof1_min;
    double tof2_best_global = tof2_min;
    double tof3_best_global = tof3_min;
    double thetaStart_best_global = thetaStart_min;
    double thetaFinal_best_global = thetaFinal_min;
    double r1_best_global = r1_min;
    double theta1_best_global = theta1_min;
    double r2_best_global = r2_min;
    double theta2_best_global = theta2_min;

    // en m
    constexpr double pos_terre_x = AU_to_meters(1.);
    // en m
    constexpr double pos_terre_y = 0;

    // Pas de discrétisation
    double d1 = (thetaFinalPlanet_max - thetaFinalPlanet_min) / NthFinalPlanet;
    double d2 = (tof1_max - tof1_min) / Nt1;
    double d3 = (tof2_max - tof2_min) / Nt2;
    double d4 = (tof3_max - tof3_min) / Nt3;
    double d5 = (thetaStart_max - thetaStart_min) / NthStart;
    double d6 = (thetaFinal_max - thetaFinal_min) / NthFinal;
    double d7 = (r1_max - r1_min) / Nr1;
    double d8 = (theta1_max - theta1_min) / Nth1;
    double d9 = (r2_max - r2_min) / Nr2;
    double d10 = (theta2_max - theta2_min) / Nth2;

    // Calcul du nombre d'itérations entières pour la boucle OpenMP
    int iters_outer = static_cast<int>(NthFinalPlanet);

    #pragma omp parallel
    {
        // Variables locales au thread
        double delta_v_best_local = std::numeric_limits<double>::infinity();
        double thetaFinalPlanet_best_local = thetaFinalPlanet_min;
        double tof1_best_local = tof1_min;
        double tof2_best_local = tof2_min;
        double tof3_best_local = tof3_min;
        double thetaStart_best_local = thetaStart_min;
        double thetaFinal_best_local = thetaFinal_min;
        double r1_best_local = r1_min;
        double theta1_best_local = theta1_min;
        double r2_best_local = r2_min;
        double theta2_best_local = theta2_min;

        // Parcours parallèle sur la boucle externe avec un index entier
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < iters_outer; ++i) {

            double thetaFinalPlanet = thetaFinalPlanet_min + i * d1;

            std::array<double, 3> pos_finalPlanet;
            {
                double cos_theta = std::cos(thetaFinalPlanet);
                double sin_theta = std::sin(thetaFinalPlanet);

                pos_finalPlanet[0] = cos_theta * distFinalSun;
                pos_finalPlanet[1] = sin_theta * distFinalSun;
                pos_finalPlanet[2] = 0;
            }

            for (double tof1 = tof1_min; tof1 <= tof1_max; tof1 += d2)
                for (double tof2 = tof2_min; tof2 <= tof2_max; tof2 += d3)
                    for (double tof3 = tof3_min; tof3 <= tof3_max; tof3 += d4)
                        for (double thetaStart = thetaStart_min; thetaStart < thetaStart_max; thetaStart += d5)
                        {
                            std::array<double, 3> p_intermediate_start;
                            {
                                double cos_theta = std::cos(thetaStart);
                                double sin_theta = std::sin(thetaStart);

                                p_intermediate_start[0] = cos_theta * RStart + pos_terre_x;
                                p_intermediate_start[1] = sin_theta * RStart + pos_terre_y;
                                p_intermediate_start[2] = 0;
                            }

                            for (double r1 = r1_min; r1 <= r1_max; r1 += d7)
                                for (double theta1 = theta1_min; theta1 <= theta1_max; theta1 += d8)
                                {
                                    std::array<double, 3> p_initial;
                                    {
                                        double cos_theta = std::cos(theta1);
                                        double sin_theta = std::sin(theta1);

                                        p_initial[0] = cos_theta * r1;
                                        p_initial[1] = sin_theta * r1;
                                        p_initial[2] = 0;
                                    }

                                    std::array<double, 3> p_inter_start_local = { std::cos(thetaStart) * RStart, std::sin(thetaStart) * RStart, 0 };

                                    auto [traj1_start, traj1_final] = lambert_universal(p_initial, p_inter_start_local, tof1, mu_start);



                                    std::array<double, 3> v_circ_init;
                                    {
                                        double speed = std::sqrt(mu_start / r1);
                                        glm::dvec2 u_r = glm::normalize(glm::dvec2(
                                            p_initial[0] - pos_terre_x,
                                            p_initial[1] - pos_terre_y
                                        ));

                                        glm::dvec2 velocity = speed * glm::dvec2(-u_r.y, u_r.x);

                                        v_circ_init[0] = velocity.x;
                                        v_circ_init[1] = velocity.y;
                                        v_circ_init[2] = 0;
                                    }

                                    double delta_v_start = std::sqrt(
                                        (traj1_start[0] - v_circ_init[0]) * (traj1_start[0] - v_circ_init[0]) +
                                        (traj1_start[1] - v_circ_init[1]) * (traj1_start[1] - v_circ_init[1]) +
                                        (traj1_start[2] - v_circ_init[2]) * (traj1_start[2] - v_circ_init[2])
                                    );

                                    for (double thetaFinal = thetaFinal_min; thetaFinal <= thetaFinal_max; thetaFinal += d6)
                                    {
                                        std::array<double, 3> p_intermediate_final;
                                        {
                                            double cos_theta = std::cos(thetaFinal);
                                            double sin_theta = std::sin(thetaFinal);

                                            p_intermediate_final[0] = cos_theta * RFinal + pos_finalPlanet[0];
                                            p_intermediate_final[1] = sin_theta * RFinal + pos_finalPlanet[1];
                                            p_intermediate_final[2] = 0;
                                        }

                                        for (double r2 = r2_min; r2 <= r2_max; r2 += d9)
                                            for (double theta2 = theta2_min; theta2 <= theta2_max; theta2 += d10)
                                            {
                                                std::array<double, 3> p_final;
                                                {
                                                    double cos_theta = std::cos(theta2);
                                                    double sin_theta = std::sin(theta2);

                                                    p_final[0] = cos_theta * r2;
                                                    p_final[1] = sin_theta * r2;
                                                    p_final[2] = 0;
                                                }

                                                // Arc 2 : héliocentrique
                                                auto [traj2_start, traj2_final] = lambert_universal(
                                                    p_intermediate_start, p_intermediate_final, tof2, mu_sun);

                                                // Vitesse héliocentrique de la Terre
                                                double v_earth_y = std::sqrt(mu_sun / pos_terre_x);

                                                // Δv au SOI Terre
                                                double delta_v_soi_earth = std::sqrt(
                                                    (traj2_start[0] - (traj1_final[0] + 0.0)) * (traj2_start[0] - (traj1_final[0] + 0.0)) +
                                                    (traj2_start[1] - (traj1_final[1] + v_earth_y)) * (traj2_start[1] - (traj1_final[1] + v_earth_y)) +
                                                    (traj2_start[2] - (traj1_final[2])) * (traj2_start[2] - (traj1_final[2]))
                                                );

                                                // Vitesse héliocentrique de la planète finale
                                                double v_fp_speed = std::sqrt(mu_sun / distFinalSun);
                                                double v_fp_x = -std::sin(thetaFinalPlanet) * v_fp_speed;
                                                double v_fp_y = std::cos(thetaFinalPlanet) * v_fp_speed;

                                                // Arc 3 : approche de la planète finale
                                                std::array<double, 3> p_inter_final_local = { std::cos(thetaFinal) * RFinal, std::sin(thetaFinal), 0 };

                                                auto [traj3_start, traj3_final] = lambert_universal(
                                                    p_inter_final_local, p_final, tof3, mu_final);

                                                // Δv au SOI planète finale
                                                double delta_v_soi_final = std::sqrt(
                                                    (traj2_final[0] - (traj3_start[0] + v_fp_x)) * (traj2_final[0] - (traj3_start[0] + v_fp_x)) +
                                                    (traj2_final[1] - (traj3_start[1] + v_fp_y)) * (traj2_final[1] - (traj3_start[1] + v_fp_y)) +
                                                    (traj2_final[2] - (traj3_start[2])) * (traj2_final[2] - (traj3_start[2]))
                                                );

                                                // Vitesse de l'orbite circulaire finale
                                                double v_circ_final_speed = std::sqrt(mu_final / r2);
                                                glm::dvec2 u_r_final = glm::normalize(glm::dvec2(
                                                    p_final[0],
                                                    p_final[1]
                                                ));
                                                std::array<double, 3> v_circ_final{};
                                                v_circ_final[0] = -u_r_final.y * v_circ_final_speed;
                                                v_circ_final[1] = u_r_final.x * v_circ_final_speed;
                                                v_circ_final[2] = 0.0;

                                                // Δv arrivée : insertion en orbite circulaire
                                                double delta_v_arrival = std::sqrt(
                                                    (traj3_final[0] - v_circ_final[0]) * (traj3_final[0] - v_circ_final[0]) +
                                                    (traj3_final[1] - v_circ_final[1]) * (traj3_final[1] - v_circ_final[1]) +
                                                    (traj3_final[2] - v_circ_final[2]) * (traj3_final[2] - v_circ_final[2])
                                                );

                                                // ── Mise à jour du meilleur LOCAL ──────────────────────────────────────
                                                double delta_v_total = delta_v_start + delta_v_soi_earth + delta_v_soi_final + delta_v_arrival;

                                                if (delta_v_total < delta_v_best_local) {
                                                    delta_v_best_local = delta_v_total;
                                                    thetaFinalPlanet_best_local = thetaFinalPlanet;
                                                    tof1_best_local = tof1;
                                                    tof2_best_local = tof2;
                                                    tof3_best_local = tof3;
                                                    thetaStart_best_local = thetaStart;
                                                    thetaFinal_best_local = thetaFinal;
                                                    r1_best_local = r1;
                                                    theta1_best_local = theta1;
                                                    r2_best_local = r2;
                                                    theta2_best_local = theta2;
                                                }
                                            }
                                    }
                                }
                        }
        } // Fin de la boucle #pragma omp for

        // Mise à jour globale sécurisée (une seule fois par thread à la fin de son travail)
        #pragma omp critical
        {
            if (delta_v_best_local < delta_v_best_global) {
                delta_v_best_global = delta_v_best_local;
                thetaFinalPlanet_best_global = thetaFinalPlanet_best_local;
                tof1_best_global = tof1_best_local;
                tof2_best_global = tof2_best_local;
                tof3_best_global = tof3_best_local;
                thetaStart_best_global = thetaStart_best_local;
                thetaFinal_best_global = thetaFinal_best_local;
                r1_best_global = r1_best_local;
                theta1_best_global = theta1_best_local;
                r2_best_global = r2_best_local;
                theta2_best_global = theta2_best_local;
            }
        }
    } // Fin du bloc #pragma omp parallel

    // Remplissage du tableau de résultat avec les valeurs globales
    res[0] = delta_v_best_global;
    res[1] = thetaFinalPlanet_best_global;
    res[2] = tof1_best_global;
    res[3] = tof2_best_global;
    res[4] = tof3_best_global;
    res[5] = thetaStart_best_global;
    res[6] = thetaFinal_best_global;
    res[7] = r1_best_global;
    res[8] = theta1_best_global;
    res[9] = r2_best_global;
    res[10] = theta2_best_global;

    return res;
}





glm::dvec3 forceAttractionGrav(const SimuCore::Structures::Entity& from, const SimuCore::Structures::Entity& on) {
    const glm::dvec3 r = from.position - on.position;
    const double distSqr = glm::dot(r, r) + SimuCore::constants::softening * SimuCore::constants::softening;
    const double dist = std::sqrt(distSqr);
    const double F = SimuCore::constants::G * from.mass * on.mass / distSqr;

    return F / dist * r; // produit scalaire par une direction unitaire
}



std::pair<double, double> calcul_perige_et_apoge (
    double distance_to_central_body,
    double system_mass,
    double mu_central_body,
    double constante_des_aires,
    double system_energy,
    bool* is_trajectory_elliptic
) {
    if (system_energy >= 0) {
        *is_trajectory_elliptic = false;
        return { -1, -1 };
    }

    double C2 = constante_des_aires * constante_des_aires;
    double mu2 = mu_central_body * mu_central_body;
    double minimum_energie_effective = -0.5 * system_mass * (mu2 / C2);

    if (system_energy < minimum_energie_effective) {
        std::cerr << "Erreur, ce n'est pas cense etre possible physiquement |\tLIGNE :" << __LINE__ << " | fichier :" << __FILE__ << std::endl;
        std::abort();
    }


	double energy_times_2 = 2 * system_energy;
    double delta = system_mass * (system_mass * mu2 + energy_times_2 * C2);
	double sqrt_delta = std::sqrt(delta);
    
    double r1 = -1, r2 = -1;

	double inverse_energy_times_2 = 1 / energy_times_2;
	double moins_mu_fois_masse = - mu_central_body * system_mass;

    double r_min = (moins_mu_fois_masse + sqrt_delta) * inverse_energy_times_2;
    double r_max = (moins_mu_fois_masse - sqrt_delta) * inverse_energy_times_2;

	if (r_min > r_max) {
        std::cerr << "Erreur, r_min > r_max |\tLIGNE :" << __LINE__ << " | fichier :" << __FILE__ << std::endl;
        std::abort();
	}

	if (r_min < 0) {
        std::cerr << "Erreur, r_min < 0 |\tLIGNE :" << __LINE__ << " | fichier :" << __FILE__ << std::endl;
        std::abort();
	}

	if (r_max < 0) {
        std::cerr << "Erreur, r_max < 0 |\tLIGNE :" << __LINE__ << " | fichier :" << __FILE__ << std::endl;
        std::abort();
	}

	*is_trajectory_elliptic = true;
	return { r_min, r_max };
}
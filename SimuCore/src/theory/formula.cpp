#include "pch.h"
#include "theory/formula.h"
#include "constants.h"
#ifdef _OPENMP
#include <omp.h>
#endif


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

            double norm_v1 = std::sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
            double norm_v2 = std::sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
            deltaV[idx] = norm_v1 + norm_v2;
        }
        catch (...) {
            // en cas d'exception, garder NaN
            deltaV[idx] = std::numeric_limits<double>::quiet_NaN();
        }
    }

    return deltaV;
}




glm::dvec3 forceAttractionGrav(const SimuCore::Structures::Entity& from, const SimuCore::Structures::Entity& on) {
    const glm::dvec3 r = from.position - on.position;
    const double distSqr = glm::dot(r, r) + SimuCore::constants::softening * SimuCore::constants::softening;
    const double dist = std::sqrt(distSqr);
    const double F = SimuCore::constants::G * from.mass * on.mass / distSqr;

    return F / dist * r; // produit scalaire par une direction unitaire
}

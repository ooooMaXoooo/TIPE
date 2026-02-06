#include "pch.h"
#include <SimuCore/theory/formula.h>
#include <SimuCore/constants.h>

#pragma comment(lib, "SimuCore.lib") 

namespace py = pybind11;


#ifdef _DEBUG
#define MODULE_NAME TIPE_SimuOrbit_d
#else
#define MODULE_NAME TIPE_SimuOrbit
#endif

PYBIND11_MODULE(MODULE_NAME, m) {
    // Ajout du sous-module 'stumpff'
    pybind11::module_ stumpff_mod = m.def_submodule("stumpff", "Fonctions de Stumpff");
    stumpff_mod.def("C", &stumpff::C, R"pbdoc(
        Fonction Stumpff C(z)
        Calcule C(z) selon la valeur de z (positive, négative ou proche de zéro)
        )pbdoc"
    );

    stumpff_mod.def("S", &stumpff::S, R"pbdoc(
        Fonction Stumpff S(z)
        Calcule S(z) avec traitement des cas limites
        )pbdoc"
    );

    py::module_ orbit_mod = m.def_submodule("orbit", "Fonctions orbitales");
    orbit_mod.def("lambert_universal", &orbit::lambert_universal,
        py::arg("r1"), py::arg("r2"), py::arg("tof"), py::arg("mu") = SimuCore::constants::mu,
        "Résout le problème de Lambert universel");

    orbit_mod.def("lambert_batch", &orbit::lambert_batch,
        py::arg("r1"), py::arg("r2_list"), py::arg("tof_list"), py::arg("mu") = SimuCore::constants::mu,
        "Calcule DeltaV (||v1||+||v2||) pour une liste de (r2, tof) en parallèle (OpenMP)");

    orbit_mod.def("lambert_batch_numpy",
        [](const std::array<double, 3>& r1,
            const std::vector<std::array<double, 3>>& r2_list,
            const std::vector<double>& tof_list,
            double mu) -> py::array_t<double>
        {
            size_t N = r2_list.size();
            if (N != tof_list.size()) {
                throw std::invalid_argument("r2_list et tof_list doivent avoir la même taille");
            }

            // Crée directement un numpy.ndarray de taille N
            py::array_t<double> result(N);
            auto r = result.mutable_unchecked<1>();  // accès direct (1D)

            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(N); i++) {
                try {
                    auto res = orbit::lambert_universal(r1, r2_list[i], tof_list[i], mu);
                    const auto& v1 = res.first;
                    const auto& v2 = res.second;

                    if (std::isfinite(v1[0]) && std::isfinite(v2[0])) {

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
                        std::array<double, 3> vecteur_vitesse_orbite_fin = calcul_vecteur_cercle(r2_list[i]);


                        r(i) = norme(substract(v1, vecteur_vitesse_orbite_depart)) + norme(substract(vecteur_vitesse_orbite_fin, v2));
                    }
                    else {
                        r(i) = NAN;
                    }
                }
                catch (...) {
                    r(i) = NAN;
                }
            }
            return result;
        },
        py::arg("r1"), py::arg("r2_list"), py::arg("tof_list"), py::arg("mu") = SimuCore::constants::mu,
        "Version optimisée: retourne directement un numpy.ndarray (sans copie)"
    );

}

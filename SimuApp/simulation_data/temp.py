try:
    import TIPE_SimuOrbit as cpp  # version Release (rapide)
    cpp_version = "Release"
except ImportError:
    try:
        import TIPE_SimuOrbit_d as cpp  # fallback vers Debug (pour dev)
        cpp_version = "Debug"
    except ImportError:
        raise ImportError("\x1B[38;5;202m Neither TIPE_SimuOrbit.pyd nor TIPE_SimuOrbit_d.pyd found \033[0m")

import Donnees_astres
import ImportData
import ATS
import numpy as np
import matplotlib.pyplot as plt
import cmath

#plt.style.use('dark_background')


def simulate_individual (dossier, gen, verbose=False) :

    base_filepath = dossier + "/gen_" + str(gen)

    (
      tof, max_time,
      thetaMin, thetaMax, thetaMean,
      NbTours,
      RMin, RMax, RMean,
      startPlanet, finalPlanet,
      p0, v0,
      finalPlanetStartIndice,
      imp0, impulsions, dates,
      genome
   ) = ImportData.lire_donnees_new(base_filepath + "_rockets_data.txt")

    X_rocket, Y_rocket = cpp.genetic.GetIndividualTrajectory__configD32_2_2(
        dt_seconds = 3600,
        max_time_days = max_time,
        simulation_duration_days= max_time*2,
        mass_kg = 1e3,
        startPlanet = startPlanet,
        finalPlanet = finalPlanet,
        genome = genome
    )

    X_rocket, Y_rocket = ATS.convert_list_to_km(X_rocket, Y_rocket)
    
    fig, ax = ATS.affiche_fichier(dossier, gen, "Jsp", verbose)

    ax.plot(X_rocket, Y_rocket, '+')


dossier = "simu_31_05_2026_17_19_38"
gen = 257

simulate_individual(dossier, gen, verbose=True)

plt.show()
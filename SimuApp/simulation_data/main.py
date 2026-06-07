from ATS import affiche_fichier, affiche_couple_generations
import matplotlib.pyplot as plt
import os
import simu_ind

plt.close("all")

dossier = "simu_07_06_2026_15_43_01"
generation = 31
verbose = False

fig, ax = affiche_fichier(dossier, generation, window_title="Terre_jupiter", verbose=verbose)
simu_ind.simulate_individual(dossier, generation, 0, fig, ax, verbose=verbose, is_best=True, trajColor="#9A228C")

plt.show()
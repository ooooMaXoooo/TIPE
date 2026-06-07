from ATS import affiche_fichier, affiche_couple_generations
import matplotlib.pyplot as plt
import os
import simu_ind

plt.close("all")

dossier = "simu_07_06_2026_16_24_29"
generation = 301
verbose = True

fig, ax = affiche_fichier(dossier, generation, window_title="Terre_jupiter", verbose=verbose)
simu_ind.simulate_individual(dossier, generation, 0, fig, ax, verbose=verbose, is_best=True, trajColor="#9A228C", alpha=1.3)

plt.show()
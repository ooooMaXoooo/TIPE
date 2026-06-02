from ATS import affiche_fichier, affiche_couple_generations
import matplotlib.pyplot as plt
import os

plt.close("all")

dossier = "simu_02_06_2026_17_50_01"
generation = 123

fig, ax = affiche_fichier(dossier, generation, window_title="Terre_jupiter")
plt.show()
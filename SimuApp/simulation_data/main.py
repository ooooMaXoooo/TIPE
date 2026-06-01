from ATS import affiche_fichier, affiche_couple_generations
import matplotlib.pyplot as plt
from graphics import barre_chargement
import os

plt.close("all")

dossier = "simu_31_05_2026_23_18_12"
generation = 15091

fig, ax = affiche_fichier(dossier, generation, window_title="Terre_jupiter")
plt.show()
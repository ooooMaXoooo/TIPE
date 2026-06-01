from ATS import affiche_couple_generations
import matplotlib.pyplot as plt
from graphics import barre_chargement
import os
import simu_ind

plt.close("all")

dir1 = "simu_01_06_2026_22_37_31"
dir2 = "simu_01_06_2026_22_37_31"

gen1 = 300
gen2 = 300

idx1 = 1
idx2 = 2

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6 * 2, 6 * 1))

fig, ax1 = simu_ind.simulate_individual(dir1, gen1, idx1, fig, ax1, verbose=False)
fig, ax2 = simu_ind.simulate_individual(dir2, gen2, idx2, fig, ax2, verbose=False)


plt.show()
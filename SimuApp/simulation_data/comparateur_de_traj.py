from ATS import affiche_couple_generations
from stats import distance_generation
import matplotlib.pyplot as plt
from graphics import barre_chargement
import os

plt.close("all")

dir1 = "simu_18_05_2026_15_57_44__cool"
dir2 = "simu_19_05_2026_11_11_35"

gen1 = 42
gen2 = 1


chercher = True
distance = distance_generation(gen1, gen2, dir1, dir2)
max_gen2 = 200

while chercher and (distance > 12 or distance < 0) and (gen2 <= max_gen2):
    distance = distance_generation(gen1, gen2, dir1, dir2)
    gen2 += 1

if gen2 > max_gen2 :
    raise "Pas de traj trouvee"


fig, _, _ = affiche_couple_generations(dir1, dir2, (gen1, gen2), window_title=f"Score : {distance:.2f} | {(gen1, gen2)}", verbose=False)

plt.show()
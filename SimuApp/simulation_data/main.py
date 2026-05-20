from ATS import affiche_fichier, affiche_couple_generations
from stats import distance_generation
import matplotlib.pyplot as plt
from graphics import barre_chargement
import os

plt.close("all")

affiche_traj_simple = False

if affiche_traj_simple :
    dossier = "simu_17_05_2026_15_12_55__super"
    generation = 9

    fig, ax = affiche_fichier(dossier, generation, window_title="Terre_jupiter")
    plt.show()

else :
    dir1 = "simu_19_05_2026_10_54_59"
    dir2 = "simu_18_05_2026_14_37_28__super"
    
    max_generation_1 = 200
    max_generation_2 = 180



    d_max = -1
    d_min = 10_000_000_000

    i_j_min = -1, -1
    i_j_max = -1, -1

    for i in range(1, max_generation_1+1) :
        for j in range(1, max_generation_2+1) :
            assert(distance_generation(i, j, dir1, dir2) == distance_generation(j, i, dir2, dir1))

            dist = distance_generation(i, j, dir1, dir2)

            if (dist > d_max) :
                d_max = dist
                i_j_max = i, j
            if(dist < d_min) :
                d_min = dist
                i_j_min = i, j
    
    print(f"dist min : {d_min}, atteint pour (i, j) : {i_j_min}")
    print(f"dist max : {d_max}, atteint pour (i, j) : {i_j_max}")


    fig1, _, _ = affiche_couple_generations(dir1, dir2, i_j_min, window_title=f"Min - {d_min:.2f}", verbose=False)
    fig2, _, _ = affiche_couple_generations(dir1, dir2, i_j_max, window_title=f"Max - {d_max:.2f}", verbose=False)

    output_dir = "distances_trajectoires_distv1/" + dir1 + "-" + dir2 + '/'

    i_min, j_min = i_j_min
    i_max, j_max = i_j_max
    os.makedirs(output_dir, exist_ok=True) # création du dossier
    fig1.savefig(output_dir + f"min_{i_min}_{j_min}-{d_min:.3f}.png", format='png')
    fig2.savefig(output_dir + f"max_{i_max}_{j_max}-{d_max:.3f}.png", format='png')

    plt.show()
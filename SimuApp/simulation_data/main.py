from ATS import affiche_fichier
from stats import distance_generation
import matplotlib.pyplot as plt

plt.close("all")

affiche_traj_simple = False

if affiche_traj_simple :
    dossier = "simu_19_05_2026_12_10_47"
    generation = 200

    fig, ax = affiche_fichier(dossier, generation)
    plt.show()

else :
    dir1 = "simu_19_05_2026_10_54_59"
    dir2 = "simu_19_05_2026_12_10_47"

    
    gen1 = 1
    gen2 = 2

    assert(distance_generation(gen1, gen2, dir1, dir2) == distance_generation(gen2, gen1, dir2, dir1))

    fig1, ax1 = affiche_fichier(dir1, gen1, window_title="1", verbose=False)
    fig2, ax2 = affiche_fichier(dir2, gen2, window_title="2", verbose=False)

    print(f"Distance entre les trajectoires : {distance_generation(gen1, gen2, dir1, dir2)}")

    plt.show()
    

    """
    d_max = -1
    d_min = 10_000_000_000

    i_j_min = -1, -1
    i_j_max = -1, -1

    for i in range(1, 201) :
        for j in range(i+1, 201) :
            gen1 = i
            gen2 = j
            assert(distance_generation(gen1, gen2, dir1, dir2) == distance_generation(gen2, gen1, dir2, dir1))

            dist = distance_generation(gen1, gen2, dir1, dir2)

            #print(f"Distance entre les trajectoires {i} et {j} : {dist}")

            if (dist > d_max) :
                d_max = dist
                i_j_max = i, j
            if(dist < d_min) :
                d_min = dist
                i_j_min = i, j
    
    print(f"dist min : {d_min}, atteint pour (i, j) : {i_j_min}")
    print(f"dist max : {d_max}, atteint pour (i, j) : {i_j_max}")
    """

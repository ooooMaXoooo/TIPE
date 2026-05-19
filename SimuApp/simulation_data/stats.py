import numpy as np
from ImportData import *

def distance_fusees_2Impulsions(f1, f2) :
    pos_1, imp1_1, t1_1, imp2_1, t2_1 = f1
    pos_2, imp1_2, t1_2, imp2_2, t2_2 = f2

    return np.linalg.norm(pos_1 - pos_2) + np.linalg.norm(imp1_1 - imp1_2) + np.abs(t1_1 - t1_2) + np.linalg.norm(imp2_1 - imp2_2) + np.abs(t2_1 - t2_2)


def distance_generation (gen1, gen2, dir) :
    filepath_1 = dir + "\\gen_" + gen1 + "_rocket.txt"
    filepath_2 = dir + "\\gen_" + gen2 + "_rocket.txt"

    with open(filepath_1) as f1:
        pos_1 = lire_vecteur(f1.readline().strip())

    with open(filepath_2) as f2:
        pos_2 = lire_vecteur(f2.readline().strip())


    tof, max_time, dt, delta_v1, planete_depart, planete_arrivee, dimension, pos_init_depart1, pos_init_arrivee1, pos_final_depart1, pos_final_arrivee1, vitesse_finale1, impulsions1, est_etat_lie1, r_min1, r_max1 = lire_donnees(dir + "/gen_" + str(gen1) + "_physics.txt")
    tof, max_time, dt, delta_v2, planete_depart, planete_arrivee, dimension, pos_init_depart2, pos_init_arrivee2, pos_final_depart2, pos_final_arrivee2, vitesse_finale2, impulsions2, est_etat_lie2, r_min2, r_max2 = lire_donnees(dir + "/gen_" + str(gen2) + "_physics.txt")
    assert(dimension == 2)

    f1 = (pos_1, impulsions1[0][1], impulsions1[0][0], impulsions1[1][1], impulsions1[1][0])
    f2 = (pos_2, impulsions2[0][1], impulsions2[0][0], impulsions2[1][1], impulsions2[1][0])
    
    return distance_fusees_2Impulsions(f1, f2)
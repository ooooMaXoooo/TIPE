import numpy as np
import cmath
from ImportData import *
from 

def dist_position(pos1, pos2, norme_max, alpha=10) : 
    return alpha * np.linalg.norm(pos1 - pos2) / norme_max

def norme_temps (instant, temps_max) :
    return instant / temps_max

def dist_temps (t1, tmax1, t2, tmax2, alpha = 5) :
    return alpha * np.abs(norme_temps(t1, tmax1) - norme_temps(t2, tmax2))

def distance_fusees_2Impulsions(f1, f2, tmax1, tmax2, norme_max_pos) :
    pos_1, imp1_1, t1_1, imp2_1, t2_1 = f1
    pos_2, imp1_2, t1_2, imp2_2, t2_2 = f2

    norme_t1 = dist_temps(t1_1, tmax1, t1_2, tmax2, alpha=5)
    norme_t2 = dist_temps(t2_1, tmax1, t2_2, tmax2, alpha=5)
    norme_temps = norme_t1 + norme_t2

    norme_position = dist_position(pos_1, pos_2, norme_max_pos, alpha=10)

    return norme_temps + norme_position + np.linalg.norm(imp1_1 - imp1_2)  + np.linalg.norm(imp2_1 - imp2_2)

def tourne_vecteur (vec, cos_angle, sin_angle) :
    v = [0,0]
    v[0] = vec[0] * cos_angle - vec[1] * sin_angle
    v[1] = vec[0] * sin_angle + vec[1] * cos_angle
    return v

def rotate_rocket(rocket, pos_terre) :
    pos, imp1, t1, imp2, t2 = rocket

    _, angle_terre = cmath.polar(complex(pos_terre[0], pos_terre[1]))
    angle_rad = - angle_terre
	
	# Rotation autour de l'axe z
    cos_angle = np.cos(angle_rad)
    sin_angle = np.sin(angle_rad)

    pos = tourne_vecteur(pos, cos_angle, sin_angle)
    imp1 = tourne_vecteur(imp1, cos_angle, sin_angle)
    imp2 = tourne_vecteur(imp2, cos_angle, sin_angle)

    return pos, imp1, t1, imp2, t2






def distance_generation (gen1, gen2, dir1, dir2) :
    filepath_1 = dir1 + "/gen_" + str(gen1) + "_rocket.txt"
    filepath_2 = dir2 + "/gen_" + str(gen2) + "_rocket.txt"

    with open(filepath_1) as f1:
        pos_1 = lire_vecteur(f1.readline().strip())
        pos_1 = np.array(pos_1)

    with open(filepath_2) as f2:
        pos_2 = lire_vecteur(f2.readline().strip())
        pos_2 = np.array(pos_2)


    tof, max_time, dt, delta_v1, planete_depart, planete_arrivee, dimension, pos_init_depart1, pos_init_arrivee1, pos_final_depart1, pos_final_arrivee1, vitesse_finale1, impulsions1, est_etat_lie1, r_min1, r_max1 = lire_donnees(dir1 + "/gen_" + str(gen1) + "_physics.txt")
    tof, max_time, dt, delta_v2, planete_depart, planete_arrivee, dimension, pos_init_depart2, pos_init_arrivee2, pos_final_depart2, pos_final_arrivee2, vitesse_finale2, impulsions2, est_etat_lie2, r_min2, r_max2 = lire_donnees(dir2 + "/gen_" + str(gen2) + "_physics.txt")
    assert(dimension == 2)

    imp1_1 = np.array(impulsions1[0][1])
    t1_1 = impulsions1[0][0]
    imp2_1 = np.array(impulsions1[1][1])
    t2_1 = impulsions1[1][0]

    imp1_2 = np.array(impulsions2[0][1])
    t1_2 = impulsions2[0][0]
    imp2_2 = np.array(impulsions2[1][1])
    t2_2 = impulsions2[1][0]

    f1 = (pos_1, imp1_1, t1_1, imp2_1, t2_1)
    f2 = (pos_2, impulsions2[0][1], impulsions2[0][0], impulsions2[1][1], impulsions2[1][0])

    rotate_rocket(f1, pos_init_depart1)
    rotate_rocket(f2, pos_init_depart2)
    
    return distance_fusees_2Impulsions(f1, f2)
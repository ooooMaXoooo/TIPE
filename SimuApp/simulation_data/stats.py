import numpy as np
import cmath
from ImportData import *
from Donnees_astres import *


def dist_position(pos1, pos2, norme_max, delta=10) : 
    return 0.5 * delta * np.linalg.norm(pos1 - pos2) / norme_max

def norme_temps (instant, temps_max) :
    return instant / temps_max

def dist_temps (t1, tmax1, t2, tmax2, tau = 5) :
    return tau * np.abs(norme_temps(t1, tmax1) - norme_temps(t2, tmax2)) / 2

def dist_impulsion(imp1, imp2, norme_max, beta=15) :
    return beta * np.linalg.norm(imp1 - imp2) / (2 * norme_max)

def norme_max_position(planete_depart) :    # en km
    return rayon_max(planete_depart, k)

def get_alpha_norme_imp2(nbImpulsion) :
    return np.clip(0.3 / nbImpulsion, 0.02, 0.15)

def distance_fusees_2Impulsions(f1, f2, tmax1, tmax2, norme_max_pos, norme_max_imp1, norme_max_imp2, verbose=False) :
    pos_1, imp1_1, imp2_1, t1 = f1
    pos_2, imp1_2, imp2_2, t2 = f2

    norme_temps = dist_temps(t1, tmax1, t2, tmax2, tau=5) # compris entre 0 et tau

    norme_pos = dist_position(pos_1, pos_2, norme_max_pos, delta=10) # compris entre 0 et delta
    
    norme_imp1 = dist_impulsion(imp1_1, imp1_2, norme_max_imp1, beta=15) # compris entre 0 et beta
    norme_imp2 = dist_impulsion(imp2_1, imp2_2, norme_max_imp2, beta=15) # compris entre 0 et beta

    if verbose :
        print(f"norme_temps : {norme_temps}")
        print(f"norme_pos : {norme_pos}")
        print(f"norme_imp1 : {norme_imp1}")
        print(f"norme_imp2 : {norme_imp2}")

    return norme_temps + norme_pos + norme_imp1 + norme_imp2

def tourne_vecteur (vec, cos_angle, sin_angle) :
    v = [0,0]
    v[0] = vec[0] * cos_angle - vec[1] * sin_angle
    v[1] = vec[0] * sin_angle + vec[1] * cos_angle
    return np.array(v)

def rotate_rocket(rocket, pos_terre) :
    pos, imp1, imp2, t = rocket

    _, angle_terre = cmath.polar(complex(pos_terre[0], pos_terre[1]))
    angle_rad = - angle_terre
	
	# Rotation autour de l'axe z
    cos_angle = np.cos(angle_rad)
    sin_angle = np.sin(angle_rad)

    pos = tourne_vecteur(pos, cos_angle, sin_angle)
    imp1 = tourne_vecteur(imp1, cos_angle, sin_angle)
    imp2 = tourne_vecteur(imp2, cos_angle, sin_angle)

    return pos, imp1, imp2, t



def distance_generation (gen1, gen2, dir1, dir2, verbose=False) :
    filepath_1 = dir1 + "/gen_" + str(gen1) + "_rocket.txt"
    filepath_2 = dir2 + "/gen_" + str(gen2) + "_rocket.txt"

    with open(filepath_1) as f1:
        pos_1 = lire_vecteur(f1.readline().strip())
        pos_1 = np.array(pos_1)

    with open(filepath_2) as f2:
        pos_2 = lire_vecteur(f2.readline().strip())
        pos_2 = np.array(pos_2)


    tof1, max_time1, dt, delta_v1, planete_depart, planete_arrivee, dimension, pos_init_depart1, pos_init_arrivee1, pos_final_depart1, pos_final_arrivee1, vitesse_finale1, impulsions1, est_etat_lie1, r_min1, r_max1 = lire_donnees(dir1 + "/gen_" + str(gen1) + "_physics.txt")
    tof2, max_time2, dt, delta_v2, planete_depart, planete_arrivee, dimension, pos_init_depart2, pos_init_arrivee2, pos_final_depart2, pos_final_arrivee2, vitesse_finale2, impulsions2, est_etat_lie2, r_min2, r_max2 = lire_donnees(dir2 + "/gen_" + str(gen2) + "_physics.txt")
    assert(dimension == 2)

    imp1_1 = np.array(impulsions1[0][1])
    imp2_1 = np.array(impulsions1[1][1])
    t1 = impulsions1[1][0]

    imp1_2 = np.array(impulsions2[0][1])
    imp2_2 = np.array(impulsions2[1][1])
    t2 = impulsions2[1][0]


    f1 = (pos_1, imp1_1, imp2_1, t1)
    f2 = (pos_2, imp1_2, imp2_2, t2)

    f1 = rotate_rocket(f1, pos_init_depart1)
    f2 = rotate_rocket(f2, pos_init_depart2)
    
    # On passe les positions dans le referentiel de la planete de depart
    pos_1, imp1_1, imp2_1, t1 = f1
    pos_2, imp1_2, imp2_2, t2 = f2
    
    pos_terre = np.array([1, 0]) # position de la terre en AU après la rotation
    
    pos_1 = pos_1 - pos_terre
    pos_2 = pos_2 - pos_terre
    
    pos_1 *= AU # en km
    pos_2 *= AU # en km
    
    
    norme_max_imp1_1 = 2 * vitesse_de_liberation(np.linalg.norm(pos_1 * 1e3), Masses[3] * G) * 1e-3
    norme_max_imp1_2 = 2 * vitesse_de_liberation(np.linalg.norm(pos_2 * 1e3), Masses[3] * G) * 1e-3
    norme_max_imp1 = max(norme_max_imp1_1, norme_max_imp1_2)
    
    norme_max_imp2 = get_alpha_norme_imp2(max(len(impulsions1), len(impulsions2))) * vitesse_orbite_circulaire(AU*1e3, mu_soleil) * 1e-3
    
    norme_max_pos = norme_max_position(planete_depart)

    if verbose :
        print(f"norme_max_imp1 : {norme_max_imp1} km/s")
        print(f"norme_max_imp2 : {norme_max_imp2} km/s")
        print(f"norme_max_pos : {norme_max_pos} km")
        print(f"||pos_1|| : {np.linalg.norm(pos_1)} km")
        print(f"||pos_2|| : {np.linalg.norm(pos_2)} km")
        print(f"planete_depart : {planete_depart} ")

    
    # on annule potentiellement la deuxième impulsion si sa date est après tof
    if (t1 > tof1) :
        imp2_1 = np.array([0,0])
        # à voir si on change les borne des distances

    if (t2 > tof2) :
        imp2_2 = np.array([0,0])
        # à voir si on change les borne des distances
    


    f1 = (pos_1, imp1_1, imp2_1, t1)
    f2 = (pos_2, imp1_2, imp2_2, t2)
    
    return distance_fusees_2Impulsions(f1, f2, max_time1, max_time2, norme_max_pos, norme_max_imp1, norme_max_imp2, verbose)

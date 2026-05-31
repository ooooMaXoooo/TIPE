import numpy as np
import matplotlib.pyplot as plt

import ImportData
import graphics
from Donnees_astres import *

def affiche_console (tof, max_time, dt, delta_v, planete_depart, planete_arrivee, dimension, pos_init_depart, pos_init_arrivee, pos_final_depart, pos_final_arrivee, vitesse_final, impulsions, est_etat_lie, r_min, r_max, nbPointsTraj) :
    def get_vec_string (vec) :
        return f"({vec[0]}, {vec[1]})"

    def get_planet_name (id) :
        names = ["Soleil", "Mercure", "Venus", "Terre", "Mars", "Jupiter", "Saturne", "Uranus", "Neptune"]
        return names[id]
    
    print(f"tof:\t{tof} days")
    print(f"delta v:\t{delta_v} km/s")
    print(f"temps de simulation max:\t{max_time} days")
    print(f"dt:\t{dt} sec")
    print(f"Start planet :\t{get_planet_name(planete_depart)}")
    print(f"Final planet :\t{get_planet_name(planete_arrivee)}")
    print(f"dimension:\t{dimension}")
    print(f"Start Planet initial position:\t{get_vec_string(pos_init_depart)} UA")
    print(f"Start Planet final position:\t{get_vec_string(pos_final_depart)} UA")
    print(f"Final Planet initial position:\t{get_vec_string(pos_init_arrivee)} UA")
    print(f"Final Planet final position:\t{get_vec_string(pos_final_arrivee)} UA")
    print(f"Rocket's final velocity:\t{get_vec_string(vitesse_final)} km/s")
    print(f"Etat lie ? :\t{"True" if est_etat_lie else "False"}")
    print(f"Rayon minimum de l'encadrement de l'ellipse final de la fusee:\t{r_min} km")
    print(f"Rayon maximum de l'encadrement de l'ellipse final de la fusee:\t{r_max} km")

    i=1
    for t, dv in impulsions:
        index = int((t * 86400) / dt)

        if index < nbPointsTraj:
            print(f"delta_v {i} de norme {np.linalg.norm(dv)} km/s")
            i+=1


def convert_vec_to_km (vec_ua) :
    return tuple([x * AU for x in vec_ua])

def convert_list_to_km(xs, ys):
    return [x * AU for x in xs], [y * AU for y in ys]

def conversion_data (pos_init_depart, pos_final_depart, pos_init_arrivee, pos_final_arrivee, x_depart, y_depart, x_final, y_final, x_fusee, y_fusee) :
    pos_init_depart = convert_vec_to_km(pos_init_depart)
    pos_final_depart = convert_vec_to_km(pos_final_depart)
    pos_init_arrivee = convert_vec_to_km(pos_init_arrivee)
    pos_final_arrivee = convert_vec_to_km(pos_final_arrivee)
    x_depart, y_depart = convert_list_to_km(x_depart, y_depart)
    x_final, y_final = convert_list_to_km(x_final, y_final)
    x_fusee, y_fusee = convert_list_to_km(x_fusee, y_fusee)

    return (pos_init_depart, pos_final_depart, pos_init_arrivee, pos_final_arrivee, x_depart, y_depart, x_final, y_final, x_fusee, y_fusee)

def affiche_orbite(fig, ax, pos_init_depart, pos_final_depart, pos_init_arrivee, pos_final_arrivee, x_depart, y_depart, x_final, y_final, x_fusee, y_fusee, planete_depart, planete_arrivee, est_etat_lie, r_min, r_max, impulsions, dt, generation) :
    ax.plot(x_depart, y_depart, '-', color='blue',label='Trajectoire planète départ')
    ax.plot(x_final, y_final, '-', color='red',label='Trajectoire planète final')
    ax.plot(x_fusee, y_fusee, '+', color='#07b00a',label='Trajectoire fusée')
    ax.plot(0, 0, 'yo', label='Soleil')  # centre

    graphics.DrawRing(planete_depart, k, ax, "Anneau_planete_depart_position_initiale", "skyblue", pos_init_depart)
    graphics.DrawRing(planete_depart, k, ax, "Anneau_planete_depart_position_final", "blue", pos_final_depart)
    graphics.DrawRing(planete_arrivee, k, ax, "Anneau_planete_final_position_initiale", "lightcoral", pos_init_arrivee)
    graphics.DrawRing(planete_arrivee, k, ax, "Anneau_planete_final_position_final", "red", pos_final_arrivee)

    if (est_etat_lie) :
        graphics.DrawRealRing(r_min, r_max, pos_final_arrivee, ax, "#0A680C", "encadrement ellipse final")


    graphics.affiche_impulsions(impulsions, dt, x_fusee, y_fusee, ax)

    ax.plot([], [], 'ko', label='Impulsions')
    ax.plot([], [], color='green', label='Delta-v')

    ax.plot(x_fusee[-1], y_fusee[-1], 'o')


    # Ajuster les limites pour afficher entièrement le grand cercle
    ax.set_aspect('equal')
    ax.set_title("Trajectoire des objets du système - gen " + str(generation))
    ax.set_xlabel("x (km)")
    ax.set_ylabel("y (km)")
    ax.grid(False)
    #ax.legend()
    return fig, ax

def _affiche_generation_sur_ax(fig, ax, dossier, generation, verbose=False):
    x_depart, y_depart, x_final, y_final, x_fusee, y_fusee = [], [], [], [], [], []

    ImportData.init_liste(dossier + "/start_planet.txt", x_depart, y_depart)
    ImportData.init_liste(dossier + "/final_planet.txt", x_final, y_final)
    ImportData.init_liste(dossier + "/gen_" + str(generation) + "_rocket.txt", x_fusee, y_fusee)
    tof, max_time, dt, delta_v, planete_depart, planete_arrivee, dimension, pos_init_depart, pos_init_arrivee, pos_final_depart, pos_final_arrivee, vitesse_final, impulsions, est_etat_lie, r_min, r_max = ImportData.lire_donnees(dossier + "/gen_" + str(generation) + "_physics.txt")

    if verbose:
        affiche_console(tof, max_time, dt, delta_v, planete_depart, planete_arrivee, dimension, pos_init_depart, pos_init_arrivee, pos_final_depart, pos_final_arrivee, vitesse_final, impulsions, est_etat_lie, r_min, r_max, len(x_fusee))

    pos_init_depart, pos_final_depart, pos_init_arrivee, pos_final_arrivee, x_depart, y_depart, x_final, y_final, x_fusee, y_fusee = conversion_data(pos_init_depart, pos_final_depart, pos_init_arrivee, pos_final_arrivee, x_depart, y_depart, x_final, y_final, x_fusee, y_fusee)

    affiche_orbite(fig, ax, pos_init_depart, pos_final_depart, pos_init_arrivee, pos_final_arrivee, x_depart, y_depart, x_final, y_final, x_fusee, y_fusee, planete_depart, planete_arrivee, est_etat_lie, r_min, r_max, impulsions, dt, generation)


def affiche_grille(dossiers_generations, layout, window_title="Trajectoires", verbose=False):
    """
    dossiers_generations : liste de tuples (dossier, generation)
    layout               : tuple (nrows, ncols)
    """
    nrows, ncols = layout
    assert len(dossiers_generations) <= nrows * ncols, "Plus de graphiques que de cases dans le layout"

    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 6 * nrows))
    fig.canvas.manager.set_window_title(window_title)

    # Normalise axes en liste 1D quelle que soit la forme du layout
    axes_flat = np.array(axes).flatten()

    for ax, (dossier, generation) in zip(axes_flat, dossiers_generations):
        _affiche_generation_sur_ax(fig, ax, dossier, generation, verbose=verbose)

    # Masque les axes inutilisés si le layout n'est pas rempli
    for ax in axes_flat[len(dossiers_generations):]:
        ax.set_visible(False)

    fig.tight_layout()
    return fig, axes_flat

def affiche_fichier(dossier, generation, window_title="Trajectoire", verbose=True):
    fig, axes = affiche_grille([(dossier, generation)], layout=(1, 1), window_title=window_title, verbose=verbose)
    return fig, axes[0]

def affiche_couple_generations(dir1, dir2, couple_gen, window_title="Title", verbose=False):
    entries = [(dir1, couple_gen[0]), (dir2, couple_gen[1])]
    fig, axes = affiche_grille(entries, layout=(1, 2), window_title=window_title, verbose=verbose)
    return fig, axes[0], axes[1]
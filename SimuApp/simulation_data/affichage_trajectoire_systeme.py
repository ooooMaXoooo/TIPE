# Jules METAIREAU
# TIPE
# 12/02/2026

import numpy as np
import matplotlib.pyplot as plt

import ImportData
import graphics
from Donnees_astres import *


dossier = "simu_04_04_2026_16_13_42"
generation = 17501




x_depart = []
y_depart = []

x_finale = []
y_finale = []

x_fusee = []
y_fusee = []


#init_liste("simu_jour_mois_annee_heure_minute_seconde/file", x_fusee, y_fusee)
                
ImportData.init_liste(dossier + "/start_planet.txt", x_depart, y_depart)
ImportData.init_liste(dossier + "/final_planet.txt", x_finale, y_finale)
ImportData.init_liste(dossier + "/gen_" + str(generation) + "_rocket.txt", x_fusee, y_fusee)




# Recuperation des données
tof, dt, planete_depart, planete_arrivee, dimension, pos_init_depart, pos_init_arrivee, pos_final_depart, pos_final_arrivee, vitesse_finale, impulsions = ImportData.lire_donnees(dossier + "/gen_" + str(generation) + "_physics.txt")

# Affichage anneaux

def get_vec_string (vec) :
    return f"({vec[0]}, {vec[1]})"

def get_planet_name (id) :
   names = ["Soleil", "Mercure", "Venus", "Terre", "Mars", "Jupiter", "Saturne", "Uranus", "Neptune"]
   return names[id]

print(f"tof:\t{tof} days")
print(f"dt:\t{dt} sec")
print(f"Start planet :\t{get_planet_name(planete_depart)}")
print(f"Final planet :\t{get_planet_name(planete_arrivee)}")
print(f"dimension:\t{dimension}")
print(f"Start Planet initial position:\t{get_vec_string(pos_init_depart)} UA")
print(f"Start Planet final position:\t{get_vec_string(pos_final_depart)} UA")
print(f"Start Planet initial position:\t{get_vec_string(pos_init_arrivee)} UA")
print(f"Final Planet final position:\t{get_vec_string(pos_final_arrivee)} UA")
print(f"Rocket's final velocity:\t{get_vec_string(vitesse_finale)} km/s")


AU = 1.496e8  # km

def convert_vec_to_km (vec_ua) :
    return tuple([x * AU for x in vec_ua])

def convert_list_to_km(xs, ys):
    return [x * AU for x in xs], [y * AU for y in ys]


pos_init_depart = convert_vec_to_km(pos_init_depart)
pos_final_depart = convert_vec_to_km(pos_final_depart)
pos_init_arrivee = convert_vec_to_km(pos_init_arrivee)
pos_final_arrivee = convert_vec_to_km(pos_final_arrivee)
x_depart, y_depart = convert_list_to_km(x_depart, y_depart)
x_finale, y_finale = convert_list_to_km(x_finale, y_finale)
x_fusee, y_fusee = convert_list_to_km(x_fusee, y_fusee)

# Affichage de l'orbite
fig, ax = plt.subplots(figsize=(6, 6))

ax.plot(x_depart, y_depart, '-', color='blue',label='Trajectoire planète départ')
ax.plot(x_finale, y_finale, '-', color='red',label='Trajectoire planète finale')
ax.plot(x_fusee, y_fusee, '-', color='#07b00a',label='Trajectoire fusée')
ax.plot(0, 0, 'yo', label='Soleil')  # centre

graphics.DrawRing(planete_depart, k, ax, "Anneau_planete_depart_position_initiale", "skyblue", pos_init_depart)
graphics.DrawRing(planete_depart, k, ax, "Anneau_planete_depart_position_finale", "blue", pos_final_depart)
graphics.DrawRing(planete_arrivee, k, ax, "Anneau_planete_finale_position_initiale", "lightcoral", pos_init_arrivee)
graphics.DrawRing(planete_arrivee, k, ax, "Anneau_planete_finale_position_finale", "red", pos_final_arrivee)


# Affichage impulsions

for t, dv in impulsions:
    index = int((t * 86400) / dt)

    if index < len(x_fusee):
        x = x_fusee[index]
        y = y_fusee[index]

        # point d'impulsion
        ax.plot(x, y, 'ko')  # point noir

        # vecteur impulsion
        scale = 1  # ajuste selon visibilité

        ax.quiver(
            x, y,
            dv[0]*scale, dv[1]*scale,
            angles='xy',
            scale_units='xy',
            color='green'
        )

ax.plot([], [], 'ko', label='Impulsions')
ax.plot([], [], color='green', label='Delta-v')





# Ajuster les limites pour afficher entièrement le grand cercle
ax.set_aspect('equal')
ax.set_title("Trajectoire des objets du système - gen " + str(generation))
ax.set_xlabel("x (km)")
ax.set_ylabel("y (km)")
ax.grid(True)
#ax.legend()
plt.show()
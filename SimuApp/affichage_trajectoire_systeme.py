# Jules METAIREAU
# TIPE
# 12/02/2026

#import numpy as np
import matplotlib.pyplot as plt

x_depart = []
y_depart = []

x_finale = []
y_finale = []

x_fusee = []
y_fusee = []

def init_liste(nom, x, y):
    with open(nom, "r", encoding="utf-8") as f:
        for ligne in f:
            ligne = ligne.strip()
            if ligne:  # ignore lignes vides
                x_val, y_val = ligne.split(";")
                x.append(float(x_val))
                y.append(float(y_val))
                
init_liste("simulation_data/simu_13_02_2026_23_21_29_gen_1_final.txt", x_depart, y_depart)
init_liste("simulation_data/simu_13_02_2026_23_21_29_gen_1_start.txt", x_finale, y_finale)
init_liste("simulation_data/simu_13_02_2026_23_23_04_gen_401_rocket.txt", x_fusee, y_fusee)

# simu_jour_mois_annee_heure_minute_seconde_gen_id_rocket.txt

# Affichage de l'orbite
fig, ax = plt.subplots(figsize=(8, 8))

ax.plot(x_depart, y_depart, '-', color='blue',label='Trajectoire planète départ')

ax.plot(x_finale, y_finale, '-', color='red',label='Trajectoire planète finale')

ax.plot(x_fusee, y_fusee, '-', color='#07b00a',label='Trajectoire fusée')

ax.plot(0, 0, 'yo', label='Soleil')  # centre

"""
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
"""


# Ajuster les limites pour afficher entièrement le grand cercle
ax.set_aspect('equal')
ax.set_title("Trajectoire des objets du système")
ax.set_xlabel("x (UA)")
ax.set_ylabel("y (UA)")
ax.grid(True)
ax.legend()
plt.show()
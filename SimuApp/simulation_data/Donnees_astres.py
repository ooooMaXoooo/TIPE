import numpy as np

#                           Soleil   Mercure    Venus      Terre     Mars     Jupiter   Saturne   Uranus   Neptune
Masses               = [1.988E+30,   3.3E+23,  4.87E+24, 5.97E+24, 6.42E+23,  1.9E+27, 5.68E+26, 8.68E+25, 1.02E+26]    # en kg
Rayons_astres        = [   696340,      2440,      6050,     6378,     3400,    71500,    60300,    25600,    24800]    # en km
Exobases             = [        0,         0,       175,      450,      200,     2000,     3500,     6000,     3000]    # en km
Distances_au_soleil  = [        0,   5.81E+07, 1.08E+08, 1.50E+08, 2.28E+08, 7.80E+08, 1.43E+09, 2.88E+09, 4.50E+09]    # en km


k = 100 # facteur de domination : l'attraction d'une planète sur l'attraction du soleil




def altitude_max(planete, k): # en km
   """
   Donne l'altitude maximale autorisée autour d'une planète donnée afin que l'attraction de la planète soit plus forte que k fois l'attraction du soleil.
   Le résulat est en km.
   """
   return ((Distances_au_soleil[planete] * 1000) / (1 + np.sqrt(k * Masses[0] / Masses[planete])) - (Rayons_astres[planete] * 1000)) * 0.001

def epaisseur_anneau(planete, k):   # en km
   """
   Donne l'épaisseur de l'anneau dont l'altidude minimal est l'exobase et l'altitude maximale est celle de la fonction qui la calcule.
   Le résultat est en km.
   """
   return altitude_max(planete, k) - Exobases[planete]

def rayon_max(planete, k):  # en km
   """
   Calcul le rayon correspondant à l'altitude maximale (voir la fonction liée).
   Le résultat est en km.
   """
   return altitude_max(planete, k) + Rayons_astres[planete]
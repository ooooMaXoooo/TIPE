import matplotlib.pyplot as plt
from Donnees_astres import *


def DrawRing(planete, k, ax, nom, couleur, position): 
   x, y = position
   ax.plot(x, y, color = couleur, marker = 'o')  # centre
   
   anneau_min = plt.Circle(position, Rayons_astres[planete] + Exobases[planete], color=couleur, linestyle='-', fill=False, label=nom+'_min')
   anneau_max = plt.Circle(position, Rayons_astres[planete] + altitude_max(planete, k), color=couleur, linestyle='-', fill=False, label=nom+'_max')
   ax.add_artist(anneau_min)
   ax.add_artist(anneau_max)
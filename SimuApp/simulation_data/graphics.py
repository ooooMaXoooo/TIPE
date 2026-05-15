import matplotlib.pyplot as plt
from Donnees_astres import *

def DrawRealRing (r_min, r_max, pos, ax, color, name) :
   anneau_min = plt.Circle(pos, r_min, color=color, linestyle='-', fill=False, label=name+'_min')
   anneau_max = plt.Circle(pos, r_max, color=color, linestyle='-', fill=False, label=name+'_max')
   ax.add_artist(anneau_min)
   ax.add_artist(anneau_max)


def DrawRing(planete, k, ax, nom, couleur, position): 
   x, y = position
   ax.plot(x, y, color = couleur, marker = 'o')  # centre
   
   anneau_min = plt.Circle(position, Rayons_astres[planete] + Exobases[planete], color=couleur, linestyle='-', fill=False, label=nom+'_min')
   anneau_max = plt.Circle(position, Rayons_astres[planete] + altitude_max(planete, k), color=couleur, linestyle='-', fill=False, label=nom+'_max')
   ax.add_artist(anneau_min)
   ax.add_artist(anneau_max)
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

def affiche_impulsions (impulsions, dt, x_fusee, y_fusee, ax) :
   for t, dv in impulsions:
      index = int((t * 86400) / dt)

      if index < len(x_fusee):

         x = x_fusee[index]
         y = y_fusee[index]

         # point d'impulsion
         ax.plot(x, y, 'ko')  # point noir

         # vecteur impulsion
         scale = 3*1e7  # ajuste selon visibilité

         ax.quiver(
               x, y,
               dv[0]*scale, dv[1]*scale,
               angles='xy',
               scale_units='xy',
               scale=1,
               color='green'
         )

def barre_chargement (value, max_value, size, fill_char='#', empty_char='_') :
   if (value >= max_value) :
      print() # sautement de ligne
   else :
      proportion = value / max_value
      nb_filled = int(proportion * size)

      print('\r', end='')
      for _ in range(nb_filled) :
         print(fill_char, end='')
      for _ in range(size-nb_filled) :
         print(empty_char, end='')

      print(f" - {proportion * 100}%", end='')
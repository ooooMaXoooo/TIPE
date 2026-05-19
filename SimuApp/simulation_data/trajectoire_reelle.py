import numpy as np

def erreur(x, y, z, nom):
    r_xy = np.sqrt(np.array(x)**2 + np.array(y)**2)
    z_error = np.abs(z)
    
    print(nom)
    print("Erreur max : ", np.max(z_error / r_xy))
    print("Erreur moyenne : ", np.mean(z_error / r_xy))
    print("Erreur :", np.mean(z_error) / np.mean(r_xy))
    
"""
Terre : #2383c2
Mars : #ff0038
Jupiter : #ecc46c
Saturne : #cc8c45
Uranus : #30b7d2
Neptune : #095197

Trajectoire : #0cda51
"""

def croisement(x, y, R, ax, color): # determiner l'indice de la position qui croise le cercle de rayon R et l'affiche sur le graphe
    r = np.sqrt(np.array(x)**2 + np.array(y)**2)
    diff = np.abs(r - R)
    i = np.argmin(diff)
    ax.plot(x[i], y[i], 'o',color=color)
   
#croisement(x_1, y_1, Distances_au_soleil[3], ax, "black")

def augmenter_fenetre(x, y, ax):
    rmax = int(np.max(np.abs(x + y)))
    marge = rmax * 0.5  # 5% de marge
    ax.set_xlim(-rmax - marge, rmax + marge)
    ax.set_ylim(-rmax - marge, rmax + marge)

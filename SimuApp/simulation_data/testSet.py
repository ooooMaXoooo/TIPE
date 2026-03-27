from Donnees_astres import *

def donnee_initiale(planete_depart, planete_finale):
   position_initiale_depart = (0, Distances_au_soleil[planete_depart])
   position_initiale_finale = (Distances_au_soleil[planete_finale], 0)
   position_finale_depart = (- Distances_au_soleil[planete_depart], 0)
   position_finale_finale = (0, Distances_au_soleil[planete_finale])
   return position_initiale_depart, position_initiale_finale, position_finale_depart, position_finale_finale
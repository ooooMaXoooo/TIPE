def init_liste(filepath, x, y):
   """
   Ouvre un fichier et lit les coordonnées des vecteurs écrit dedans.
   Précondition : le fichier a besoin d'être formatté avec un vecteur par ligne, et ses coordonnées séparées d'un seul point-virgule
   """
   with open(filepath, "r", encoding="utf-8") as f:
      for ligne in f:
         ligne = ligne.strip()
         if ligne:  # ignore lignes vides
            x_val, y_val = ligne.split(";")
            x.append(float(x_val))
            y.append(float(y_val))


def donnees_physiques(filepath, x, y):
   with open(filepath, "r", encoding="utf-8") as f:
      for ligne in f:
         ligne = ligne.strip()
         if ligne:  # ignore lignes vides
            x_val, y_val = ligne.split(";")
            x.append(float(x_val))
            y.append(float(y_val))


def lire_vecteur(ligne):
   return tuple(float(x) for x in ligne.split(";"))


def lire_donnees(nom_fichier):
   with open(nom_fichier, "r", encoding="utf-8") as f:
      lignes = [ligne.strip() for ligne in f if ligne.strip()]

   i = 0

   tof = float(lignes[i])
   i += 1
   dt = float(lignes[i])
   i += 1
   planete_depart = int(lignes[i])
   i += 1
   planete_arrivee = int(lignes[i])
   i += 1
   dimension = int(lignes[i])
   i += 1
   nb_impulsions = int(lignes[i])
   i += 1

   pos_init_depart = lire_vecteur(lignes[i])
   i += 1
   pos_init_arrivee = lire_vecteur(lignes[i])
   i += 1
   pos_final_depart = lire_vecteur(lignes[i])
   i += 1
   pos_final_arrivee = lire_vecteur(lignes[i])
   i += 1

   vitesse_finale = lire_vecteur(lignes[i])
   i += 1

   impulsions = []
   for _ in range(nb_impulsions):
      temps = float(lignes[i])
      i += 1
      vect = lire_vecteur(lignes[i])
      i += 1

      impulsions.append((temps, vect))

   return (tof, dt, planete_depart, planete_arrivee, dimension,
         pos_init_depart, pos_init_arrivee,
         pos_final_depart, pos_final_arrivee,
         vitesse_finale, impulsions)
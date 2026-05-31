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

   max_time = float(lignes[i])
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

   etat_lie = int(lignes[i]) == 1
   i += 1

   r_min = float(lignes[i])
   i += 1

   r_max = float(lignes[i])
   i += 1

   tof = float(lignes[i])
   i += 1

   delta_v = float(lignes[i])
   i += 1

   return (tof, max_time, dt, delta_v, planete_depart, planete_arrivee, dimension,
         pos_init_depart, pos_init_arrivee,
         pos_final_depart, pos_final_arrivee,
         vitesse_finale, impulsions, etat_lie, r_min, r_max)

def load_horizons_file(filepath, x_list, y_list):
   """
   Charge un fichier Horizons format :
   GEOMETRIC cartesian states

   et remplit :
   x_list
   y_list

   Les coordonnées sont supposées en KM.
   """

   with open(filepath, "r", encoding="utf-8") as f:
      lines = f.readlines()

   in_data = False

   for line in lines:

      line = line.strip()

      # Début données
      if line.startswith("$$SOE"):
         in_data = True
         continue

      # Fin données
      if line.startswith("$$EOE"):
         break

      if not in_data:
         continue

      # Ignore lignes vides / headers
      if not line or line.startswith("*"):
         continue

      # Format CSV Horizons :
      #
      # JD, DATE, X, Y, Z,
      #
      parts = [p.strip() for p in line.split(",")]

      if len(parts) < 5:
         continue

      try:
         x = float(parts[2])
         y = float(parts[3])

         x_list.append(x)
         y_list.append(y)

      except:
         continue


def lire_donnees_new(nom_fichier):
   with open(nom_fichier, "r", encoding="utf-8") as f:
      lignes = [ligne.strip() for ligne in f if ligne.strip()]

   i = 0

   tof = float(lignes[i])
   i += 1

   max_time = float(lignes[i])
   i += 1


   thetaMin = float(lignes[i])
   i += 1
   thetaMax = float(lignes[i])
   i += 1
   thetaMean = float(lignes[i])

   i += 1
   NbTours = float(lignes[i])
   i += 1

   RMin = float(lignes[i])
   i += 1
   RMax = float(lignes[i])
   i += 1
   RMean = float(lignes[i])
   i+=1

   startPlanet = int(lignes[i])
   i+=1
   finalPlanet = int(lignes[i])
   i+=1

   p0 = list(lire_vecteur(lignes[i]))
   i += 1
   v0 = list(lire_vecteur(lignes[i]))
   i += 1

   finalPlanetStartIndice = int(lignes[i])
   i+=1

   imp0 = list(lire_vecteur(lignes[i]))
   i+=1

   nbImpulsions = int(lignes[i])
   i+=1

   impulsions, dates = [], []
   for j in range(nbImpulsions-1) :
      impulsions.append(list(lire_vecteur(lignes[i])))
      i+=1

      dates.append(float(lignes[i]))
      i+=1

   nbChromosomes = int(lignes[i])
   i += 1
   nbGenes = int(lignes[i])
   i += 1


   genome = []
   for k in range(nbChromosomes) :
      chromo = [int(x) for x in lignes[i].split(";")]
      genome.append(chromo)
      print(chromo)
      i+=1



   return (
      tof, max_time,
      thetaMin, thetaMax, thetaMean,
      NbTours,
      RMin, RMax, RMean,
      startPlanet, finalPlanet,
      p0, v0,
      finalPlanetStartIndice,
      imp0, impulsions, dates,
      genome
   )
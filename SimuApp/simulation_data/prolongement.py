try:
    import TIPE_SimuOrbit as cpp  # version Release (rapide)
    cpp_version = "Release"
except ImportError:
    try:
        import TIPE_SimuOrbit_d as cpp  # fallback vers Debug (pour dev)
        cpp_version = "Debug"
    except ImportError:
        raise ImportError("\x1B[38;5;202m Neither TIPE_SimuOrbit.pyd nor TIPE_SimuOrbit_d.pyd found \033[0m")

import Donnees_astres
import ImportData
import ATS
import numpy as np
import matplotlib.pyplot as plt
import cmath


def find_closest_indice(pos, X, Y):
    pos = np.array(pos[:2])
    _, y = pos

    if y >= 0:
        i_min, i_max = 0, len(X) // 2
    else:
        i_min, i_max = len(X) // 2 + 1, len(X) - 1

    best_idx = i_min
    best_dist = np.linalg.norm(np.array([X[i_min], Y[i_min]]) - pos)

    for i in range(i_min + 1, i_max + 1):
        d = np.linalg.norm(np.array([X[i], Y[i]]) - pos)
        if d < best_dist:
            best_dist = d
            best_idx = i

    return best_idx


# en km/s
def initial_speed (r_UA, start_planet=3) :
    # Donnees_astres.AU = 1.496e8  # km 
    r_km = (Donnees_astres.AU * np.array(r_UA)) - np.array([Donnees_astres.Distances_au_soleil[start_planet], 0]) # tableau Distances_au_soleil en km

    dist = np.linalg.norm(r_km)
    print(f"dist : {dist}")

    mu_terre = Donnees_astres.G * Donnees_astres.Masses[start_planet] # tableau Masses en kg
    norm = np.sqrt(mu_terre / (dist * 1e3)) * 1e-3

    u_r = (1/dist) * np.array(r_km)
    u_theta = np.array([-u_r[1], u_r[0], 0])


    v_terre = 1e-3 * Donnees_astres.vecteur_vitesse_orbite_circulaire(1e3 * Donnees_astres.AU * np.array([1, 0]), Donnees_astres.mu_soleil)
    v_terre = np.array([v_terre[0], v_terre[1], 0])

    return norm * u_theta + v_terre


def rotate_trajectory (theta, X, Y) :
    cos_theta, sin_theta = np.cos(theta), np.sin(theta)

    R_theta = np.array([[cos_theta, -sin_theta], [sin_theta, cos_theta]])
    Traj = np.array([X, Y])

    return np.dot(R_theta, Traj)


def rotate_vector (theta, vec) :
    cos_theta, sin_theta = np.cos(theta), np.sin(theta)

    x = vec[0]
    y = vec[1]

    new_vec = [0, 0]
    new_vec[0] = x * cos_theta - y * sin_theta
    new_vec[1] = x * sin_theta + y * cos_theta

    return new_vec







def prolongement (tempsSupplementaire, dossier, gen, verbose=False) :
    x_depart, y_depart, x_final, y_final = [], [], [], []
    ImportData.init_liste(dossier + "/start_planet.txt", x_depart, y_depart)
    ImportData.init_liste(dossier + "/final_planet.txt", x_final, y_final)

    tof, max_time, dt, _, planete_depart, planete_arrivee, _, pos_init_depart, pos_init_arrivee, pos_final_depart, pos_final_arrivee, _, impulsions, est_etat_lie, r_min, r_max = ImportData.lire_donnees(dossier + "/gen_" + str(gen) + "_physics.txt")

    x_fusee, y_fusee = [], []
    ImportData.init_liste(dossier + "/gen_" + str(gen) + "_rocket.txt", x_fusee, y_fusee)



    # Détermination de l'angle de rotation :
    pos_terre = pos_init_depart
    _, angle_terre = cmath.polar(complex(pos_terre[0], pos_terre[1]))
    theta = -angle_terre

    # rotation de la trajectoire
    traj = rotate_trajectory(theta, x_fusee, y_fusee).tolist()
    x_fusee, y_fusee = traj[0], traj[1]


    # rotation des planètes
    pos_init_depart = rotate_vector(theta, pos_init_depart)
    pos_init_arrivee = rotate_vector(theta, pos_init_arrivee)

    pos_final_depart = rotate_vector(theta, pos_final_depart)
    pos_final_arrivee = rotate_vector(theta, pos_final_arrivee)

    # rotation des impulsions
    for i in range(len(impulsions)) :
        t, dv = impulsions[i]
        dv = rotate_vector(theta, dv)

        impulsions[i] = (t, dv)


    ####### prolongement :


    indice_position = find_closest_indice(pos_init_arrivee, x_final, y_final)


    impulsions_, dates_ = [], []

    for t, dv in impulsions:
        index = int((t * 86400) / dt)

        if index < len(x_fusee):
            impulsions_.append((dv[0], dv[1], 0))
            dates_.append(t)


    initial_pos = (x_fusee[0], y_fusee[0], 0)
    vitesse_initiale = initial_speed((initial_pos[0], initial_pos[1]))
    approx_vitesse_initiale = (Donnees_astres.AU/dt) * (np.array([x_fusee[1], y_fusee[1]]) - np.array([x_fusee[0], y_fusee[0]]))
    approx_vitesse_initiale = (approx_vitesse_initiale[0], approx_vitesse_initiale[1], 0)

    """
    v0 = 29.780 * np.array([0, 1])  # vitesse héliocentrique de la Terre en (1,0)

    AU = Donnees_astres.AU * 1e3  # en mètres
    mu_start = Donnees_astres.G * Donnees_astres.Masses[planete_depart]

    ref_terre = (np.array(initial_pos) - np.array([1, 0, 0])) * AU  # Terre→fusée en mètres

    ref_terre_2d = ref_terre[:2]
    jsp = np.array([-ref_terre_2d[1], ref_terre_2d[0]])
    orbit_velocity = jsp / np.linalg.norm(jsp)
    orbit_velocity = orbit_velocity * 1e-3 * np.sqrt(mu_start / np.linalg.norm(ref_terre))  # km/s

    v0 = orbit_velocity + v0
    v0 = [v0[0], v0[1], 0]
    """

    #v0 = [np.float64(1.9416986420472568), np.float64(31.160700016016634), 0]
    v0 = [np.float64(2.0116986420472568), np.float64(31.100000016016634), 0]

    if verbose :
        print(f"dt : {dt} s")
        print(f"max_time : {max_time} j")
        print(f"tof : {tof} j")
        print(f"initial_pos : {initial_pos} UA")
        print(f"vitesse initiale : {vitesse_initiale} km/s")
        print(f"approx_vitesse_initiale : {approx_vitesse_initiale} km/s")
        print(f"v0 : {v0} km/s")


    X_test, Y_test = cpp.integrator.Simulate(
        dt_seconds=dt,
        max_time_days=tof + tempsSupplementaire,
        mass_kg=1e3,
        initial_pos_UA=initial_pos,
        initial_velocity_km_s=approx_vitesse_initiale,
        final_planet_pos_indice=indice_position,
        impulsions_km_s=impulsions_,
        dates_days=dates_,
        startPlanet=planete_depart,
        finalPlanet=planete_arrivee
    )
    X_test, Y_test = ATS.convert_list_to_km(X_test, Y_test)




    pos_init_depart, pos_final_depart, pos_init_arrivee, pos_final_arrivee, \
    x_depart, y_depart, x_final, y_final, x_fusee, y_fusee = \
        ATS.conversion_data(
            pos_init_depart, pos_final_depart, pos_init_arrivee, pos_final_arrivee,
            x_depart, y_depart, x_final, y_final, x_fusee, y_fusee
        )


    # affichage trajectoire tournée
    fig, ax = plt.subplots()
    fig, ax = ATS.affiche_orbite(
        fig, ax,
        pos_init_depart,
        pos_final_depart,
        pos_init_arrivee,
        pos_final_arrivee,
        x_depart,
        y_depart,
        x_final,
        y_final,
        x_fusee,
        y_fusee,
        planete_depart,
        planete_arrivee,
        est_etat_lie,
        r_min, r_max, 
        impulsions,
        dt, gen
    )

    ax.plot(X_test, Y_test, '+')

    pos1 = np.array([X_test[1], X_test[1]])
    pos1 = np.array([x_fusee[1], y_fusee[1]])

    print()


def prolong2 (dossier, gen, verbose=False) :
    (
      tof,
      thetaMin, thetaMax, thetaMean,
      NbTours,
      RMin, RMax, RMean,
      startPlanet, finalPlanet,
      p0, v0,
      finalPlanetStartIndice,
      imp0, impulsions, dates
   ) = ImportData.lire_donnees_new(dossier + "/gen_" + str(gen) + "_rockets_data.txt")
    
    def vec_to3D (vec) :
        return [vec[0], vec[1], 0]
    
    for i in range(len(impulsions)) :
        impulsions[i] = vec_to3D(impulsions[i])

    impulsions_ = [vec_to3D(imp0)]
    dates_ = [1e-15]

    for i in range(len(impulsions)):
        impulsions_.append(impulsions[i])
        dates_.append(dates[i])

    print(impulsions_)
    print(dates_)
    print(p0)
    print(v0)

    X_test, Y_test = cpp.integrator.Simulate(
        dt_seconds=3600,
        max_time_days=tof * 3,
        mass_kg=1e3,
        initial_pos_UA=vec_to3D(p0),
        initial_velocity_km_s=vec_to3D(v0),
        final_planet_pos_indice=finalPlanetStartIndice,
        impulsions_km_s=impulsions_,
        dates_days=dates_,
        startPlanet=startPlanet,
        finalPlanet=finalPlanet
    )
    X_test, Y_test = ATS.convert_list_to_km(X_test, Y_test)
    
    fig, ax = ATS.affiche_fichier(dossier, gen, "Jsp", verbose)

    ax.plot(X_test, Y_test, '+')





dossier = "simu_21_05_2026_23_55_59"
gen = 2

#prolongement(0, dossier, gen, True)
prolong2(dossier, gen)

plt.show()
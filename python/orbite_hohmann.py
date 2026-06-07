# Jules METAIREAU
# TIPE

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton as scipy_newton

#plt.close('all')

#__________________________________________________________#
#--------------------Calcul_Energetique--------------------#
#__________________________________________________________#

def v_instantanee_ellipse(r, a, mu):
    return np.sqrt(mu * (2 / r - 1 / a))

def v_instantanee_cercle(r, mu):
    return np.sqrt(mu / r)

def delta_v_1(r1, a1, mu):
    return np.sqrt(mu * (2 / r1 - 1 / a1)) - np.sqrt(mu / r1)

def delta_v_2(r2, a1, mu):
    return np.sqrt(mu / r2) - np.sqrt(mu * (2 / r2 - 1 / a1))

def delta_u_1(r1, a2, mu):
    return np.sqrt(mu * (2 / r1 - 1 / a2)) - np.sqrt(mu / r1)

def delta_u_2(r3, a2, a3, mu):
    return np.sqrt(mu * (2 / r3 - 1 / a3)) - np.sqrt(mu * (2 / r3 - 1 / a2))

def delta_u_3(r2, a3, mu):
    return np.sqrt(mu / r2) - np.sqrt(mu * (2 / r2 - 1 / a3))

def delta_v(r1, r2, mu):
    a1 = (r1 + r2) / 2
    return np.abs(delta_v_1(r1, a1, mu)) + np.abs(delta_v_2(r2, a1, mu))

def delta_u(r1, r2, x, mu):
    a2 = (r1 + x) / 2 
    a3 = (r2 + x) / 2
    return np.abs(delta_u_1(r1, a2, mu)) + np.abs(delta_u_2(x, a2, a3, mu)) + np.abs(delta_u_3(r2, a3, mu))

# On cherche a minimiser la variation de masse
# Il faut connaitre la masse intiale mi de la fusee,
# la variation de vitesse lors de l'impulsion
# la vitesse d'ejection ve des gaz la fusee (constante)
# Retourne la variation de masse delta_m lors de l'impulsion et la masse finale de la fusee
def tsiolkovski(mi, delta_v, ve):
    mf = mi * np.exp(- delta_v / ve)
    return (mf - mi), mf

#__________________________________________________________#
#----------------Energie-one-tangent-burn------------------#
#__________________________________________________________#
def trouver_phi(p, e, r2):
    return np.arccos((p - r2) / (e * r2))

def delta_v_otb(v_ellipse, v_cercle, alpha):
    return np.sqrt(v_ellipse**2 + v_cercle**2 - 2 * v_ellipse * v_cercle * np.cos(alpha))

def calcul_alpha(e, phi):
    return np.arctan2(e * np.sin(phi), 1 + e * np.cos(phi))

#__________________________________________________________#
#__________________________________________________________#
#------------------------Simulation------------------------#
#__________________________________________________________#
#__________________________________________________________#

#__________________________________________________________#
#------------------------Parametres------------------------#
#__________________________________________________________#

def periode(a, masse_central):
    return np.sqrt(4 * PI**2 * a**3 / (G * masse_central))

# Implementation de newton
def newton_kepler(T, t, e, psi_init):
    M = 2 * PI * t / T  # Anomalie moyenne
    def f(psi):  # Équation de Kepler
        return psi - e * np.sin(psi) - M
    def df(psi):  # Dérivée de f par rapport à psi
        return 1 - e * np.cos(psi)  
    return scipy_newton(f, psi_init, fprime=df, tol=1e-6, maxiter=10000)

def conv_psi_en_phi(e, psi):
    return 2 * np.arctan2(np.sqrt(1 + e) * np.sin(psi / 2), np.sqrt(1 - e) * np.cos(psi / 2))

def calcul_rayon(e, phi, p):
    return p / (1 + e * np.cos(phi - 0))

def calcul_parametres(r1, r2, masse_central):
    a = (r1 + r2) / 2
    e = abs((r1 - r2) / (r1 + r2))
    p = a * (1 - e**2)
    T = periode(a, masse_central)
    return a, e, p, T

#__________________________________________________________#
#-------------------------Affichage------------------------#
#__________________________________________________________#

def affichage_hohmann(cartesien_x, cartesien_y, r1, r2):
    # Affichage de l'orbite
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.plot(cartesien_x, cartesien_y, '-', color='#0cda51', linewidth = 6,label='Trajectoire')
    ax.plot(0, 0, 'o',color='#ffde21', label='Soleil')  # centre

    # Cercles de rayon r1 et r2
    cercle_r1 = plt.Circle((0, 0), r1, color='#2383c2', linewidth = 4, linestyle='--', fill=False, label=f'Orbite r1 = {r1 * 1e-3:.2e} km')
    cercle_r2 = plt.Circle((0, 0), r2, color='#ff0038', linewidth = 4, linestyle='--', fill=False, label=f'Orbite r2 = {r2 * 1e-3:.2e} km')
    ax.add_artist(cercle_r1)
    ax.add_artist(cercle_r2)

    # Ajuster les limites pour afficher entièrement le grand cercle
    rmax = max(int(np.max(np.abs(cartesien_x + cartesien_y))), r2, r1)
    marge = rmax * 0.05  # 5% de marge
    ax.set_xlim(-rmax - marge, rmax + marge)
    ax.set_ylim(-rmax - marge, rmax + marge)
    ax.set_aspect('equal')
    #ax.set_title("Trajectoire entre les orbites r1 et r2")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.grid(False)
    #ax.legend()
    ax.axis('off')
    plt.show()

def hohmann_ellipse(r1, r2, masse_central, premiere_iteration, nombre_iteration, pas, psi, a , e, p, T):
    
    cartesien_x = []
    cartesien_y = []
    temps = []
    
    for i in range(premiere_iteration, nombre_iteration):
        t = i * pas
        psi = newton_kepler(T, t, e, psi)
        phi = conv_psi_en_phi(e, psi)
        rayon = calcul_rayon(e, phi, p)
        x = rayon * np.cos(phi)
        y = rayon * np.sin(phi)
        cartesien_x.append(x)
        cartesien_y.append(y)
        temps.append(t)
        
    return cartesien_x, cartesien_y, psi

#__________________________________________________________#
#--------------------------Simple--------------------------#
#__________________________________________________________#

def simu_ellipse(r1, r2, masse_central, nombre_iteration, pas, masse_fusee):
    print("Hohmann simple")
    
    a, e, p, T = calcul_parametres(r1, r2, masse_central)
    print("a = {:.2e} km".format(a * 1e-3))
    print("e = {:.2f}".format(e))
    print("demi-periode = {:.0f} jours".format(T / (2 * 24 * 3600)))    # jours
    print("demi-periode = {:.0f} annees".format(T / (2 * 24 * 3600 * 365)))

    cartesien_x, cartesien_y, _ = hohmann_ellipse(r1, r2, masse_central, 0, int(T // (2 * pas)), pas, 0, a, e, p, T)
    affichage_hohmann(cartesien_x, cartesien_y, r1, r2)
    
    # Variation totale de vitesse
    mu = G * masse_central
    # delta_v = variation_vitesse(r1, masse_central, a) + variation_vitesse(r2, masse_central, a)
    d_v = delta_v(r1, r2, mu)
    print("Δv : {:.2f} km/s".format(d_v*1e-3))
    
    # Equation de tsiolkovski
    ve = 2000 # m / s, la valeur est plutot arbitraire (trouve sur internet)
    delta_m, mf = tsiolkovski(masse_fusee, d_v, ve)
    #print("Δm : {:.2f} kg avec son signe".format(delta_m))
    print("mi = {:.2f}".format(masse_fusee))
    print("mf = {:.2f}".format(mf))
    
    # Encadre pour comparaison
    a1 = (r1 + r2) / 2
    impulsion = delta_v_1(r1, a1, mu)
    print("")   # saut de ligne
    print("Donnee comparaison")
    print("impulsion Δv : (0, {:.4f} km/s) ".format(impulsion*1e-3))
    print("temps total de voyage : {:.4f} jours".format(T / (2 * 24 * 3600)))
    
    
#__________________________________________________________#
#-----------------------Bi_elliptique----------------------#
#__________________________________________________________#

def simu_bielliptique(r1, r2, r3, masse_central, nombre_iteration, pas, masse_fusee):
    
    # Ellipse 1
    a1, e1, p1, T1 = calcul_parametres(r1, r3, masse_central)
    print("\nEllipse 1 :")
    print("a1 = {:.2e} km".format(a1 * 1e-3))
    print("e1 = {:.2f}".format(e1))
    print("T1 = {:.0f} jours".format(T1 / (24 * 3600)))    # jours
    
    cartesien_x1, cartesien_y1, psi1 = hohmann_ellipse(r1, r3, masse_central, 0, int(T1 // (2 * pas)), pas, 0, a1, e1, p1, T1)
    
    # Ellipse 2
    a2, e2, p2, T2 = calcul_parametres(r3, r2, masse_central)
    print("\nEllipse 2 :")
    print("a2 = {:.2e} km".format(a2 * 1e-3))
    print("e2 = {:.2f}".format(e2))
    print("T2 = {:.0f} jours".format(T2 / (24 * 3600)))    # jours
    
    cartesien_x2, cartesien_y2, _ = hohmann_ellipse(r3, r2, masse_central, int(T2 // (2 * pas)), int(T2 // (pas)), pas, psi1, a2, e2, p2, T2)
    
    affichage_hohmann(cartesien_x1 + cartesien_x2, cartesien_y1 + cartesien_y2, r1, r2)
    
    # Donnees orbite totale
    print("\nHohmann bi-elliptique")

    T = (T1 / 2 + T2 / 2)
    print("temps de voyage = {:.0f} jours".format((T1 / 2 + T2 / 2) / (24 * 3600)))    # jours
    print("temps de voyage = {:.0f} annees".format((T1 / 2 + T2 / 2) / (24 * 3600 * 365)))
    
    # Variation totale de vitesse
    mu = G * masse_central
    #a3 = (r3 + r2) / 2
    d_u = delta_u(r1, r2, r3, mu)
    print("Δu : {:.2f} km/s".format(d_u*1e-3))
    
    # Encadre pour comparaison
    a2 = (r1 + r3) / 2 
    a3 = (r2 + r3) / 2
    impulsion1 = delta_u_1(r1, a2, mu)
    impulsion2 = delta_u_2(r3, a2, a3, mu)
    impulsion3 = delta_u_3(r2, a3, mu)
    print("")   # saut de ligne
    print("Donnee comparaison")
    print("impulsion1 Δu1 : (0, {:.4f} km/s) ".format(impulsion1*1e-3))
    print("impulsion2 Δu2 : (0, {:.4f} km/s) ".format(impulsion2*1e-3))
    print("impulsion2 Δu3 : (0, {:.4f} km/s) ".format(impulsion3*1e-3))
    print("a l'instant t : {:.4f} jours".format(T1 / (2 * 24 * 3600)))   
    print("temps total de voyage : {:.4f} jours".format(T / (24 * 3600)))

#__________________________________________________________#
#----------------Simple_avec_angle_initial-----------------#
#__________________________________________________________#


#__________________________________________________________#
#---------------------One_tangent_burn---------------------#
#__________________________________________________________#


def determination_phi(p, e , r):
    return np.arccos((p - r) / (e * r))

def equa_temps_kepler(T, psi, e):   # retourne l'instant
    return T * (psi - e * np.sin(psi)) / (2 * PI)

def calcul_parametres_otb(r1, r2, masse_central, a):
    c = a - r1
    e = c / a
    p = a * (1 - e**2)
    T = periode(a, masse_central)
    return a, e, p, T

def conv_phi_en_psi(e, phi):
    return 2 * np.arctan2(np.sqrt(1 - e) * np.sin(phi / 2), np.sqrt(1 + e) * np.cos(phi / 2))

def simu_one_tangent_burn(r1, r2, masse_central, nombre_iteration, pas, masse_fusee, a):
    print("\nHohmann one tangent burn")
    
    _, e, p, T = calcul_parametres_otb(r1, r2, masse_central, a)
    print("a = {:.2e} km".format(a * 1e-3))
    print("e = {:.2f}".format(e))
    print("demi-periode = {:.0f} jours".format(T / (2 * 24 * 3600)))    # jours
    print("demi-periode = {:.0f} annees".format(T / (2 * 24 * 3600 * 365)))
    
    t_init = equa_temps_kepler(T, conv_phi_en_psi(e, 0), e)
    print("t_init = ", t_init)
    t_final = equa_temps_kepler(T, conv_phi_en_psi(e, determination_phi(p, e, r2)), e)
    print("t_final = ", t_final)
    print("temps de voyage = {:.0f} jours".format((t_final - t_init) / (24 * 3600)))
    print("temps de voyage = {:.0f} annees".format((t_final - t_init) / (24 * 3600 * 365)))

    cartesien_x, cartesien_y, _ = hohmann_ellipse(r1, r2, masse_central, int(t_init // pas), int(t_final // pas), pas, 0, a, e, p, T)
    affichage_hohmann(cartesien_x, cartesien_y, r1, r2)
    
    # Variation totale de vitesse
    mu = G * masse_central
    d_v1 = v_instantanee_ellipse(r1, a, mu) - v_instantanee_cercle(r1, mu)  # norme de la vitesse
    d_v2 = delta_v_otb(v_instantanee_ellipse(r2, a, mu), v_instantanee_cercle(r2, mu), calcul_alpha(e, trouver_phi(p, e, r2)))  # norme de la vitesse
    d_vtot = np.abs(d_v1) + np.abs(d_v2)
    print("Δv : {:.2f} km/s".format(d_vtot*1e-3))

#__________________________________________________________#
#--------------------------Donnees-------------------------#
#__________________________________________________________#

# Constantes
G = 6.67430e-11  # constante gravitationnelle
PI = np.pi

UA = 150e9

# Donnees
# Venus : 0.72
# Terre : 1.0
# Mars : 1.52
# Jupiter : 5.21
# Saturne : 9.54
# Uranus : 19.18
# Neptune : 30.11

# Parametres
# ATTENTION : ces donnees sont utiles pour les fonctions simu
# MAIS r1 et r2 servent aussi pour les comparaisons des delta_v
r1 = 1.00 * UA   # m
r2 = 1.52 * UA
r3 = 5 * UA
masse_central = 1.989e30  # kg
masse_fusee = 400e3 # kg
nombre_iteration = 300 # jours
pas = 1 * 24 * 3600 # 1 jour en seconde
teta0 = 0
a = 10 * UA     # ATTENTION : a > r2 environ


simu_ellipse(r1, r2, masse_central, nombre_iteration, pas, masse_fusee)
#simu_bielliptique(r1, r2, r3, masse_central, nombre_iteration, pas, masse_fusee)
#simu_one_tangent_burn(r1, r2, masse_central, nombre_iteration, pas, masse_fusee, a)
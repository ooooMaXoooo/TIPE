import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. CONSTANTES PHYSIQUES ET PARAMÈTRES
# ==========================================
DAY_IN_SEC = 86400.0
AU = 149597870700.0  # m

# Paramètres gravitationnels (à vérifier avec ton C++)
MU_EARTH = 3.986004418e14
MU_SUN = 1.32712440042e20

# ==========================================
# 2. DONNÉES C++ INJECTÉES
# ==========================================
# Temps de vol (convertis en secondes)
tof1_sec = 1.74 * DAY_IN_SEC
tof2_sec = 298 * DAY_IN_SEC
tof3_sec = 13.8 * DAY_IN_SEC

# Angles dans le plan (Azimut / Est-Ouest)
phiFinalPlanet = 1.83 
phiStart = 0.0395        
phiFinal = 5.29        
phi1 = 4.62            
phi2 = 6.2            

# Rayons
r1 = 1.33e08       # Départ Terre
r2 = 1.17e10       # Arrivée Planète
RStart = 2.59e8    # SOI Terre
RFinal = 2.34e10   # SOI Planète d'arrivée

distFinalSun = 7.78e11   #rayon orbital planète cible

# Vitesses initiales des arcs (m/s)
traj1_start = np.array([2.42e+03, -67.3])
traj2_start = np.array([4.51e+03, 4.51e+04])
traj3_start = np.array([-818, 1.56e+04])


MU_TARGET = 1.27e+17 

# ==========================================
# 3. MOTEUR PHYSIQUE (LEAPFROG)
# ==========================================
def get_acceleration(r_vec, mu):
    r_mag = np.linalg.norm(r_vec)
    return -mu / (r_mag**3) * r_vec

def leapfrog_integrate(r0, v0, mu, tof, num_steps=15000):
    dt = tof / num_steps
    r_hist = np.zeros((num_steps, 2))
    v_hist = np.zeros((num_steps, 2))
    
    r_hist[0] = r0
    v_hist[0] = v0
    
    a_current = get_acceleration(r_hist[0], mu)
    v_half = v_hist[0] + a_current * (dt / 2.0)
    
    for i in range(num_steps - 1):
        r_hist[i+1] = r_hist[i] + v_half * dt
        a_next = get_acceleration(r_hist[i+1], mu)
        v_hist[i+1] = v_half + a_next * (dt / 2.0)
        v_half = v_half + a_next * dt
        
    return r_hist

# ==========================================
# 4. PRÉPARATION DES CONDITIONS INITIALES
# ==========================================

# --- ARC 1 : Autour de la Terre ---
# Référentiel centré sur la Terre
r0_arc1 = np.array([r1 * np.cos(phi1), r1 * np.sin(phi1)])
v0_arc1 = traj1_start

# --- ARC 2 : Héliocentrique ---
# Référentiel centré sur le Soleil
pos_terre_x, pos_terre_y = AU, 0.0
r0_arc2 = np.array([
    pos_terre_x + RStart * np.cos(phiStart),
    pos_terre_y + RStart * np.sin(phiStart)
])
v0_arc2 = traj2_start

# --- ARC 3 : Autour de la planète cible ---
# Référentiel centré sur la planète cible
r0_arc3 = np.array([RFinal * np.cos(phiFinal), RFinal * np.sin(phiFinal)])
v0_arc3 = traj3_start

# ==========================================
# 5. EXÉCUTION DES SIMULATIONS
# ==========================================
print("Intégration de l'Arc 1 (Terre)...")
traj_arc1 = leapfrog_integrate(r0_arc1, v0_arc1, MU_EARTH, tof1_sec)

print("Intégration de l'Arc 2 (Soleil)...")
traj_arc2 = leapfrog_integrate(r0_arc2, v0_arc2, MU_SUN, tof2_sec)

print("Intégration de l'Arc 3 (Planète cible)...")
traj_arc3 = leapfrog_integrate(r0_arc3, v0_arc3, MU_TARGET, tof3_sec)

# ==========================================
# 6. AFFICHAGE (3 GRAPHIQUES)
# ==========================================
plt.style.use('dark_background')
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle("Profil de Mission en 3 Phases (Coniques Raccordées)", fontsize=16)

# --- Graphique 1 : Départ Terre ---
ax1.set_title(f"Arc 1 : Évasion Terrestre ({tof1_sec/DAY_IN_SEC:.1f} jours)")
ax1.plot(0, 0, 'bo', markersize=10, label="Terre")
ax1.plot(traj_arc1[:, 0], traj_arc1[:, 1], 'g-', label="Trajectoire Arc 1")
ax1.plot(traj_arc1[0, 0], traj_arc1[0, 1], 'go', label="Départ (r1)")
ax1.plot(traj_arc1[-1, 0], traj_arc1[-1, 1], 'rx', label="Arrivée SOI")
# Cercle représentant la SOI
theta_circle = np.linspace(0, 2*np.pi, 100)
ax1.plot(RStart * np.cos(theta_circle), RStart * np.sin(theta_circle), 'w--', alpha=0.3, label="Limite SOI")
ax1.axis('equal')
ax1.grid(alpha=0.2)
ax1.legend()

# --- Graphique 2 : Transfert Héliocentrique ---
ax2.set_title(f"Arc 2 : Héliocentrique ({tof2_sec/DAY_IN_SEC:.1f} jours)")
ax2.plot(0, 0, 'yo', markersize=15, label="Soleil")
ax2.plot(traj_arc2[:, 0], traj_arc2[:, 1], 'g-', label="Trajectoire Arc 2")
ax2.plot(traj_arc2[0, 0], traj_arc2[0, 1], 'go', label="Départ (Bord SOI Terre)")
ax2.plot(traj_arc2[-1, 0], traj_arc2[-1, 1], 'rx', label="Arrivée (Bord SOI Cible)")
# Orbites des planètes pour référence
ax2.plot(AU * np.cos(theta_circle), AU * np.sin(theta_circle), 'b--', alpha=0.5, label="Orbite Terre")
ax2.plot(distFinalSun * np.cos(theta_circle), distFinalSun * np.sin(theta_circle), 'r--', alpha=0.5, label="Orbite Cible")
ax2.axis('equal')
ax2.grid(alpha=0.2)
ax2.legend()

# --- Graphique 3 : Arrivée Planète Cible ---
ax3.set_title(f"Arc 3 : Insertion ({tof3_sec/DAY_IN_SEC:.1f} jours)")
ax3.plot(0, 0, 'ro', markersize=10, label="Planète Cible")
ax3.plot(traj_arc3[:, 0], traj_arc3[:, 1], 'g-', label="Trajectoire Arc 3")
ax3.plot(traj_arc3[0, 0], traj_arc3[0, 1], 'go', label="Départ SOI Cible")
ax3.plot(traj_arc3[-1, 0], traj_arc3[-1, 1], 'rx', label="Arrivée (r2)")
# Cercle représentant la SOI d'arrivée
ax3.plot(RFinal * np.cos(theta_circle), RFinal * np.sin(theta_circle), 'w--', alpha=0.3, label="Limite SOI")
ax3.axis('equal')
ax3.grid(alpha=0.2)
ax3.legend()
plt.tight_layout()


# ==========================================
# 7. VUE GLOBALE HÉLIOCENTRIQUE (ZOOMABLE)
# ==========================================

# -- Changement de référentiel (Local -> Héliocentrique) --
# L'Arc 2 est déjà en héliocentrique.
# Pour l'Arc 1, on le décale à la position de la Terre
traj_arc1_helio = traj_arc1 + np.array([AU, 0.0])

# Pour l'Arc 3, on le décale à la position de la planète cible
target_pos_x = distFinalSun * np.cos(phiFinalPlanet)
target_pos_y = distFinalSun * np.sin(phiFinalPlanet)
traj_arc3_helio = traj_arc3 + np.array([target_pos_x, target_pos_y])

# -- Création du graphique global --
fig_global = plt.figure(figsize=(10, 10))
plt.title("Vue Globale Héliocentrique\n(Utilisez l'outil loupe pour zoomer sur les planètes et les raccords)", fontsize=14)

# Le Soleil au centre
plt.plot(0, 0, 'yo', markersize=15, label="Soleil")

# Les trajectoires superposées
plt.plot(traj_arc1_helio[:, 0], traj_arc1_helio[:, 1], 'c-', linewidth=3, label="Arc 1 (Évasion Terrestre)")
plt.plot(traj_arc2[:, 0], traj_arc2[:, 1], 'g-', linewidth=2, label="Arc 2 (Transfert Héliocentrique)")
plt.plot(traj_arc3_helio[:, 0], traj_arc3_helio[:, 1], 'm-', linewidth=3, label="Arc 3 (Insertion Planétaire)")

# Les positions des planètes
plt.plot(AU, 0, 'bo', markersize=8, label="Terre")
plt.plot(target_pos_x, target_pos_y, 'ro', markersize=8, label="Planète Cible")

# Les orbites planétaires en filigrane
plt.plot(AU * np.cos(theta_circle), AU * np.sin(theta_circle), 'b--', alpha=0.3)
plt.plot(distFinalSun * np.cos(theta_circle), distFinalSun * np.sin(theta_circle), 'r--', alpha=0.3)

# Les Sphères d'Influence (SOI)
plt.plot(AU + RStart * np.cos(theta_circle), RStart * np.sin(theta_circle), 'w:', alpha=0.6, label="SOI Terre")
plt.plot(target_pos_x + RFinal * np.cos(theta_circle), target_pos_y + RFinal * np.sin(theta_circle), 'w:', alpha=0.6, label="SOI Cible")

plt.axis('equal')
plt.grid(alpha=0.2)
plt.legend(loc="upper right")

plt.show()



plt.show()
import numpy as np
import time
import threading
import matplotlib.pyplot as plt

try:
    import TIPE_SimuOrbit  # version Release (rapide)
    cpp_version = "Release"
except ImportError:
    try:
        import TIPE_SimuOrbit_d as TIPE_SimuOrbit  # fallback vers Debug (pour dev)
        cpp_version = "Debug"
    except ImportError:
        raise ImportError("\x1B[38;5;202m Neither TIPE_SimuOrbit.pyd nor TIPE_SimuOrbit_d.pyd found \033[0m")

# Test rapide du module C++
try:
    test_r1 = [1.5e11, 0.0, 0.0]
    test_r2 = [0.0, 2.28e11, 0.0]
    test_result = TIPE_SimuOrbit.orbit.lambert_universal(test_r1, test_r2, 200*86400)
except Exception as e:
    raise



AU = 149597870700.0 # en m


# === Constantes ===
mu = 1.32712440018e20  # du soleil

R_Mercure = 0.39 * AU
R_Venus = 0.72 * AU
R_Terre = AU
R_Mars = 1.5237 * AU
R_Jupiter = 5.21 * AU
R_Saturne = 9.54 * AU
R_Uranus  = 19.18 * AU
R_Neptune = 30.11 * AU 

R_final = R_Mars
R_depart = R_Terre

theta   = 179.4     # deg
tof     = 258.3       # jours



r1 = np.array([R_depart, 0.0, 0.0])



# === fonctions ===
def clamp (x, a, b) :
    min = a if a < b else b
    max = a+b-min

    x = x if x > a else a
    x = x if x < b else b

    return x


def calcul_vecteur_vitesse_circulaire (pos) :
    r = np.linalg.norm(pos)
    norme = np.sqrt(mu / r)

    u_r = pos / r
    u_theta = np.array([-u_r[1], u_r[0], 0])

    return u_theta * norme



def compute_trajectory(theta, tof_days):
    r2 = R_final * np.array([np.cos(theta), np.sin(theta), 0.0])
    
    try:
        v1_chosen, v2_chosen = TIPE_SimuOrbit.orbit.lambert_universal(
            r1.tolist(), r2.tolist(), tof_days*86400, mu
        )
        v1_chosen = np.array(v1_chosen)
        v2_chosen = np.array(v2_chosen)
        
        if np.any(np.isnan(v1_chosen)) or np.any(np.isnan(v2_chosen)):
            result = (np.empty((0,3)), r2, np.nan)
            return result
        
        n_traj_points = 200
        dt = tof_days * 86400 / n_traj_points
        r = r1.copy()
        v = v1_chosen.copy()
        traj = [r.copy()]
        
        v_half = v + 0.5 * dt * (-mu * r / np.linalg.norm(r)**3)
        for _ in range(n_traj_points):
            r = r + dt * v_half
            a_new = -mu * r / np.linalg.norm(r)**3
            v_half = v_half + dt * a_new
            traj.append(r.copy())
        
        vitesse_depart = calcul_vecteur_vitesse_circulaire(traj[0])
        vitesse_finale = calcul_vecteur_vitesse_circulaire(traj[-1])

        deltav = np.linalg.norm(v1_chosen - vitesse_depart) + np.linalg.norm(v2_chosen - vitesse_finale)
        result = (np.array(traj), r2, deltav, v1_chosen, v2_chosen, vitesse_depart, vitesse_finale)
        return result
        
    except Exception as e:
        result = (np.empty((0,3)), r2, np.nan, np.empty((0,3)), np.empty((0,3)), np.empty((0,3)), np.empty((0,3)))
        return result

def compute_circle (r) :
    N = 1000
    X = []
    Y = []
    for i in range(N) :
        theta = 2 * np.pi * i / N
        X.append(r * np.cos(theta))
        Y.append(r * np.sin(theta))

    return (X, Y)


def plot_vecteur (pos, vec, color, label = '') :
    (x, y) = pos

    # vecteur impulsion
    scale = 4e6  # ajuste selon visibilité

    if label != '' :
        plt.quiver(
            x, y,
            vec[0]*scale, vec[1]*scale,
            angles='xy',
            scale_units='xy',
            scale=1,
            color=color,
            label = label
        )
    else : 
        plt.quiver(
            x, y,
            vec[0]*scale, vec[1]*scale,
            angles='xy',
            scale_units='xy',
            scale=1,
            color=color
        )


traj_3D, r2, deltav, vi, vf, vcd, vca = compute_trajectory(np.deg2rad(theta), tof)

X = []
Y = []

for pos in traj_3D :
    X.append(pos[0])
    Y.append(pos[1])

plt.close("all")
#plt.style.use('dark_background')

plt.plot(X, Y, linewidth=3)

plt.plot(0, 0, 'oy')
plt.plot(X[0], Y[0], 'og')
plt.plot(X[-1], Y[-1], 'or')

X_depart, Y_depart = compute_circle(np.linalg.norm([X[0], Y[0]]))
X_final, Y_final = compute_circle(np.linalg.norm([X[-1], Y[-1]]))

plt.plot(X_depart, Y_depart, ':g', linewidth=3)
plt.plot(X_final, Y_final, ':r', linewidth=3)



pi = (X[0], Y[0])
pf = (X[-1], Y[-1])

plot_vecteur(pi, vi, "#018BFC", label='vitesse voulue')
plot_vecteur(pi, vcd, "#12E647", label='vitesse initiale')

plot_vecteur(pf, vf, "#F8550A", label='vitesse finale')
plot_vecteur(pf, vca, "#CD10B4", label='vitesse voulue')



plt.axis('equal')
plt.legend()
#plt.grid()
plt.show()
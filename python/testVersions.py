import print_color as pc
import numpy as np
import time
pc.print_color("numpy imported (1/6)", pc.Color.Green)

import matplotlib.pyplot as plt
pc.print_color("matplotlib.pyplot imported (2/6)", pc.Color.Green)

from mpl_toolkits.mplot3d import Axes3D
pc.print_color("mpl_toolkits.mplot3d imported (3/6)", pc.Color.Green)

from matplotlib.widgets import Slider
pc.print_color("matplotlib.widgets::Slider imported (4/6)", pc.Color.Green)

from scipy.optimize import brentq
pc.print_color("scipy.optimize::brentq imported (5/6)", pc.Color.Green)

try:
    import TIPE_SimuOrbit  # version Release (rapide)
    pc.print_color("Using RELEASE version of TIPE_SimuOrbit (6/6)", pc.Color.Green)
    cpp_version = "Release"
except ImportError:
    try:
        import TIPE_SimuOrbit_d as TIPE_SimuOrbit  # fallback vers Debug (pour dev)
        pc.print_color("WARNING: Using DEBUG version of TIPE_SimuOrbit (6/6)", pc.Color.Yellow)
        cpp_version = "Debug"
    except ImportError:
        pc.print_color("ERROR: No version of TIPE_SimuOrbit found!", pc.Color.Red)
        raise ImportError("\x1B[38;5;202m Neither TIPE_SimuOrbit.pyd nor TIPE_SimuOrbit_d.pyd found \033[0m")

# Test rapide du module C++
try:
    test_r1 = [1.5e11, 0.0, 0.0]
    test_r2 = [0.0, 2.28e11, 0.0]
    test_result = TIPE_SimuOrbit.orbit.lambert_universal(test_r1, test_r2, 200*86400)
    pc.print_color("Module C++ test: OK", pc.Color.Green)
except Exception as e:
    pc.print_color(f"Module C++ test failed: {e}", pc.Color.Red)
    raise

# === Constantes ===
mu = 1.32712440018e20
R_Mars = 2.28e11
r1 = np.array([1.5e11, 0.0, 0.0])

# === Paramètres utilisateur ===
n_points = 500
tof_min_days, tof_max_days = 20, 500
theta_min_deg, theta_max_deg = 0, 360
tof_chosen_days, theta_chosen_deg = 300, 45

print("=== Paramètres utilisateur ===")
print(f"Nombre de points : {n_points}")
print(f"Intervalle temps (jours) : {tof_min_days} à {tof_max_days}")
print(f"Intervalle angle θ (deg) : {theta_min_deg} à {theta_max_deg}")
print(f"Temps choisi (jours) : {tof_chosen_days}")
print(f"Angle choisi θ (deg) : {theta_chosen_deg}\n")

# === Préparation des grilles ===
tof_list = np.linspace(tof_min_days*86400, tof_max_days*86400, n_points)
theta_list = np.linspace(theta_min_deg, theta_max_deg, n_points)
Theta, ToF = np.meshgrid(theta_list, tof_list/86400)

# Construction des entrées pour lambert_batch
r2_list = []
tof_flat = []

for i in range(n_points):
    for j in range(n_points):
        theta_rad = np.deg2rad(Theta[i, j])
        r2 = R_Mars * np.array([np.cos(theta_rad), np.sin(theta_rad), 0.0])
        r2_list.append(r2.tolist())
        tof_flat.append(tof_list[i])

# === Appel C++ batch ===
print("=== Calcul de la surface Δv avec le module C++ (batch OpenMP) ===")
start_time = time.time()

DeltaV_flat = TIPE_SimuOrbit.orbit.lambert_batch(
    r1.tolist(), r2_list, tof_flat, mu
)

end_time = time.time()
computation_time = end_time - start_time

DeltaV = np.array(DeltaV_flat).reshape(n_points, n_points)

valid_computations = np.sum(np.isfinite(DeltaV))
total_computations = DeltaV.size

print(f"\n=== Résultats ===")
pc.print_color(f"Calcul terminé en {computation_time:.2f} secondes!", pc.Color.Green)
pc.print_color(f"Vitesse: {total_computations/computation_time:.1f} calculs/seconde", pc.Color.Cyan)
pc.print_color(f"Solutions valides: {valid_computations}/{total_computations} "
               f"({100*valid_computations/total_computations:.1f}%)", pc.Color.Blue)

# === Figure Δv ===
fig = plt.figure(figsize=(16,7))
ax1 = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(122, projection='3d')

DeltaV_masked = np.ma.masked_invalid(DeltaV)
surf = ax1.plot_surface(Theta, ToF, DeltaV_masked, cmap='viridis', alpha=0.7)
point_chosen = ax1.scatter(theta_chosen_deg, tof_chosen_days, 0, color='red', s=50, label="Δv choisi")
ax1.set_xlabel("Angle θ (deg)")
ax1.set_ylabel("Temps (jours)")
ax1.set_zlabel("Δv (m/s)")
ax1.set_title(f"Δv(θ, temps) - C++ {cpp_version}\n"
              f"{computation_time:.2f}s, {total_computations/computation_time:.0f} calc/s")
ax1.legend()
fig.colorbar(surf, ax=ax1, shrink=0.5, aspect=10, label="Δv (m/s)")

# === Fonction trajectoire (reste identique à ton script) ===
def compute_trajectory(theta_deg, tof_days):
    theta_rad = np.deg2rad(theta_deg)
    r2 = R_Mars * np.array([np.cos(theta_rad), np.sin(theta_rad), 0.0])
    try:
        v1_chosen, v2_chosen = TIPE_SimuOrbit.orbit.lambert_universal(
            r1.tolist(), r2.tolist(), tof_days*86400, mu
        )
    except:
        return np.empty((0,3)), r2
    v1_chosen = np.array(v1_chosen)
    v2_chosen = np.array(v2_chosen)
    if np.any(np.isnan(v1_chosen)) or np.any(np.isnan(v2_chosen)):
        return np.empty((0,3)), r2

    n_traj_points = 300
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
    return np.array(traj), r2

traj_init, r2_init = compute_trajectory(theta_chosen_deg, tof_chosen_days)

# === Affichage trajectoire ===
theta_orbite = np.linspace(0, 2*np.pi, 500)
orbit_mars = R_Mars * np.array([np.cos(theta_orbite), np.sin(theta_orbite), np.zeros_like(theta_orbite)])
ax2.plot(orbit_mars[0], orbit_mars[1], orbit_mars[2], linestyle='--', color='orange', label="Orbite Mars")
if traj_init.size > 0:
    traj_line = ax2.plot(traj_init[:,0], traj_init[:,1], traj_init[:,2], color='red', label="Trajectoire fusée")[0]
else:
    traj_line = ax2.plot([], [], [], color='red', label="Trajectoire fusée")[0]
point_A = ax2.scatter(*r1, color='green', s=50, label="Terre")
point_B = ax2.scatter(*r2_init, color='blue', s=50, label="Point B")
point_Sun = ax2.scatter(0, 0, 0, color='yellow', s=100, label="Soleil")
ax2.set_xlabel("X (m)")
ax2.set_ylabel("Y (m)")
ax2.set_zlabel("Z (m)")
ax2.set_title("Trajectoire 3D fusée")
ax2.legend()

# === Sliders ===
ax_slider_theta = plt.axes([0.25, 0.02, 0.5, 0.03])
slider_theta = Slider(ax_slider_theta, "Angle θ (°)", theta_min_deg, theta_max_deg, valinit=theta_chosen_deg)
ax_slider_tof = plt.axes([0.25, 0.06, 0.5, 0.03])
slider_tof = Slider(ax_slider_tof, "Temps (jours)", tof_min_days, tof_max_days, valinit=tof_chosen_days)

def update(val):
    theta = slider_theta.val
    tof = slider_tof.val
    theta_rad = np.deg2rad(theta)
    r2 = R_Mars * np.array([np.cos(theta_rad), np.sin(theta_rad), 0.0])
    try:
        v1, v2 = TIPE_SimuOrbit.orbit.lambert_universal(r1.tolist(), r2.tolist(), tof*86400, mu)
        v1, v2 = np.array(v1), np.array(v2)
        if not (np.any(np.isnan(v1)) or np.any(np.isnan(v2))):
            cost = np.linalg.norm(v1) + np.linalg.norm(v2)
        else:
            cost = np.nan
    except:
        cost = np.nan
    point_chosen._offsets3d = ([theta], [tof], [cost])
    traj, r2_new = compute_trajectory(theta, tof)
    if traj.size > 0:
        traj_line.set_data(traj[:,0], traj[:,1])
        traj_line.set_3d_properties(traj[:,2])
        point_B._offsets3d = ([r2_new[0]], [r2_new[1]], [r2_new[2]])
    fig.canvas.draw_idle()

slider_theta.on_changed(update)
slider_tof.on_changed(update)

print("\n=== Interface prête ===")
pc.print_color("Utilisez les sliders pour explorer les trajectoires!", pc.Color.Magenta)
plt.show()

import cProfile
cProfile.run('update(0)')

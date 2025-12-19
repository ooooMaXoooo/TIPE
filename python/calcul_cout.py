import print_color as pc
import numpy as np
pc.print_color("numpy imported", pc.Color.Green)

import matplotlib.pyplot as plt
pc.print_color("matplotlib.pyplot imported", pc.Color.Green)

from mpl_toolkits.mplot3d import Axes3D
pc.print_color("mpl_toolkits.mplot3d imported", pc.Color.Green)


from matplotlib.widgets import Slider
pc.print_color("matplotlib.widgets::Slider imported", pc.Color.Green)


from scipy.optimize import brentq
pc.print_color("scipy.optimize::brentq imported", pc.Color.Green)


try:
    import TIPE_SimuOrbit  # version Release (rapide)
    pc.print_color("Using RELEASE version of TIPE_SimuOrbit", pc.Color.Green)
except ImportError:
    try:
        import TIPE_SimuOrbit_d as TIPE_SimuOrbit  # fallback vers Debug (pour dev)
        pc.print_color("WARNING: Using DEBUG version of TIPE_SimuOrbit", pc.Color.Yellow)
    except ImportError:
        pc.print_color("ERROR: No version of TIPE_SimuOrbit found!", pc.Color.Red)
        raise ImportError("\x1B[38;5;202m Neither TIPE_SimuOrbit.pyd nor TIPE_SimuOrbit_d.pyd found \033[0m")

# Maintenant vous pouvez utiliser TIPE_SimuOrbit normalement
#print("Available modules:", dir(TIPE_SimuOrbit))



# --- Constantes ---
mu = 1.32712440018e20
R_Mars = 2.28e11
r1 = np.array([1.5e11, 0.0, 0.0])

# --- Stumpff robustes ---
def stumpff_C(z):
    if z > 1e-8: return (1 - np.cos(np.sqrt(z))) / z
    elif z < -1e-8: return (np.cosh(np.sqrt(-z)) - 1) / -z
    else: return 0.5 - z/24 + z**2/720

def stumpff_S(z):
    if z > 1e-8: s=np.sqrt(z); return (s-np.sin(s))/(s**3)
    elif z < -1e-8: s=np.sqrt(-z); return (np.sinh(s)-s)/(s**3)
    else: return 1/6 - z/120 + z**2/5040

# --- Lambert universel robuste ---
def lambert_universal(r1,r2,tof,mu=mu):
    r1_norm=np.linalg.norm(r1); r2_norm=np.linalg.norm(r2)
    cos_dtheta=np.clip(np.dot(r1,r2)/(r1_norm*r2_norm),-1,1)
    dtheta=np.arccos(cos_dtheta)
    if abs(dtheta)<1e-8: return np.full(3,np.nan), np.full(3,np.nan)
    A=np.sin(dtheta)*np.sqrt(r1_norm*r2_norm/(1-cos_dtheta))
    if A==0: return np.full(3,np.nan), np.full(3,np.nan)

    def tof_of_z(z):
        C=stumpff_C(z); S=stumpff_S(z)
        if C<=0: return np.nan
        y=r1_norm+r2_norm + A*(z*S-1)/np.sqrt(C)
        if y<=0: return np.nan
        chi=np.sqrt(y/C)
        return (chi**3*S + A*np.sqrt(y))/np.sqrt(mu)

    z_grid=np.linspace(-4*np.pi**2,4*np.pi**2,1200)
    F_vals=np.array([tof_of_z(z)-tof for z in z_grid])
    bracket=None
    for i in range(len(F_vals)-1):
        if np.isfinite(F_vals[i]) and np.isfinite(F_vals[i+1]) and F_vals[i]*F_vals[i+1]<0:
            bracket=(z_grid[i],z_grid[i+1]); break
    if bracket is None: return np.full(3,np.nan), np.full(3,np.nan)
    z=brentq(lambda z_:tof_of_z(z_)-tof,bracket[0],bracket[1],xtol=1e-10)
    C=stumpff_C(z); S=stumpff_S(z)
    y=r1_norm+r2_norm + A*(z*S-1)/np.sqrt(C)
    chi=np.sqrt(y/C)
    f=1-y/r1_norm; g=A*np.sqrt(y/mu); gdot=1-y/r2_norm
    if abs(g)<1e-6: return np.full(3,np.nan), np.full(3,np.nan)
    return (r2-f*r1)/g, (gdot*r2-r1)/g

# --- Paramètres utilisateur ---
n_points = 50
tof_min_days=20; tof_max_days=500
theta_min_deg=0; theta_max_deg=360
tof_chosen_days=300; theta_chosen_deg=45

# --- Affichage console ---
print("=== Paramètres utilisateur ===")
print(f"Nombre de points : {n_points}")
print(f"Intervalle temps (jours) : {tof_min_days} à {tof_max_days}")
print(f"Intervalle angle θ (deg) : {theta_min_deg} à {theta_max_deg}")
print(f"Temps choisi (jours) : {tof_chosen_days}")
print(f"Angle choisi θ (deg) : {theta_chosen_deg}\n")
print("=== Conditions initiales ===")

print(f"Position départ r1 (Terre) : {r1}")
print(f"Rayon orbite Mars : {R_Mars}\n")

# --- Pré-calcul Δv(θ,tof) ---
tof_list=np.linspace(tof_min_days*86400,tof_max_days*86400,n_points)
theta_list=np.linspace(theta_min_deg,theta_max_deg,n_points)
Theta, ToF = np.meshgrid(theta_list, tof_list/86400)
DeltaV=np.zeros_like(Theta)

for i in range(n_points):
    for j in range(n_points):
        theta_rad=np.deg2rad(Theta[i,j])
        r2=R_Mars*np.array([np.cos(theta_rad),np.sin(theta_rad),0.0])
        v1,v2=lambert_universal(r1,r2,tof_list[i])
        DeltaV[i,j]=np.linalg.norm(v1)+np.linalg.norm(v2) if not np.any(np.isnan(v1)) else np.nan

# --- Figure et axes ---
fig = plt.figure(figsize=(14,6))
ax1 = fig.add_subplot(121, projection='3d') # Δv(θ,tof)
ax2 = fig.add_subplot(122, projection='3d') # Trajectoire

# --- Surface Δv ---
surf=ax1.plot_surface(Theta,ToF,DeltaV,cmap='viridis',alpha=0.7)
point_chosen=ax1.scatter(theta_chosen_deg,tof_chosen_days,0,color='red',s=50,label="Δv choisi")
ax1.set_xlabel("Angle θ (deg)"); ax1.set_ylabel("Temps (jours)"); ax1.set_zlabel("Δv (m/s)")
ax1.set_title("Δv(θ, temps)")
ax1.legend()
fig.colorbar(surf, ax=ax1, shrink=0.5, aspect=10, label="Δv (m/s)")

# --- Trajectoire initiale ---
def compute_trajectory(theta_deg, tof_days):
    theta_rad=np.deg2rad(theta_deg)
    r2 = R_Mars*np.array([np.cos(theta_rad),np.sin(theta_rad),0.0])
    v1_chosen, v2_chosen = lambert_universal(r1,r2,tof_days*86400)
    if np.any(np.isnan(v1_chosen)): return np.empty((0,3)), r2
    n_traj_points=300; dt=tof_days*86400/n_traj_points
    r=r1.copy(); v=v1_chosen.copy(); traj=[r.copy()]
    v_half=v+0.5*dt*(-mu*r/np.linalg.norm(r)**3)
    for _ in range(n_traj_points):
        r=r+dt*v_half
        a_new=-mu*r/np.linalg.norm(r)**3
        v_half=v_half+dt*a_new
        traj.append(r.copy())
    return np.array(traj), r2

traj_init, r2_init=compute_trajectory(theta_chosen_deg,tof_chosen_days)
theta_orbite=np.linspace(0,2*np.pi,500)
orbit_mars=R_Mars*np.array([np.cos(theta_orbite),np.sin(theta_orbite),np.zeros_like(theta_orbite)])
ax2.plot(orbit_mars[0],orbit_mars[1],orbit_mars[2],linestyle='--',color='orange',label="Orbite Mars")
traj_line=ax2.plot(traj_init[:,0],traj_init[:,1],traj_init[:,2],color='red',label="Trajectoire fusée")[0]
point_A=ax2.scatter(*r1,color='green',s=50,label="Terre")
point_B=ax2.scatter(*r2_init,color='blue',s=50,label="Point B")
point_Sun=ax2.scatter(0,0,0,color='yellow',s=100,label="Soleil")
ax2.set_xlabel("X (m)"); ax2.set_ylabel("Y (m)"); ax2.set_zlabel("Z (m)")
ax2.set_title("Trajectoire 3D fusée")
ax2.legend()

# --- Sliders ---
ax_slider_theta = plt.axes([0.25,0.02,0.5,0.03])
slider_theta = Slider(ax_slider_theta,"Angle θ (°)",theta_min_deg,theta_max_deg,valinit=theta_chosen_deg)
ax_slider_tof = plt.axes([0.25,0.06,0.5,0.03])
slider_tof = Slider(ax_slider_tof,"Temps (jours)",tof_min_days,tof_max_days,valinit=tof_chosen_days)

# --- Update ---
def update(val):
    theta = slider_theta.val
    tof = slider_tof.val

    # Δv choisi
    theta_rad = np.deg2rad(theta)
    r2 = R_Mars*np.array([np.cos(theta_rad),np.sin(theta_rad),0.0])
    v1,v2 = lambert_universal(r1,r2,tof*86400)
    cost = np.linalg.norm(v1)+np.linalg.norm(v2) if not np.any(np.isnan(v1)) else np.nan
    point_chosen._offsets3d = ([theta],[tof],[cost])

    # Trajectoire
    traj, r2_new = compute_trajectory(theta,tof)
    if traj.size>0:
        traj_line.set_data(traj[:,0],traj[:,1])
        traj_line.set_3d_properties(traj[:,2])
        point_B._offsets3d = ([r2_new[0]],[r2_new[1]],[r2_new[2]])
    fig.canvas.draw_idle()

slider_theta.on_changed(update)
slider_tof.on_changed(update)

plt.show()

# VERSION DEBUG - pour identifier les problèmes

import print_color as pc
import numpy as np
import time
import dash
from dash import dcc, html, callback, Input, Output
import plotly.graph_objects as go
import threading
import webbrowser

try:
    import TIPE_SimuOrbit
    cpp_version = "Release"
    pc.print_color("Using RELEASE version", pc.Color.Green)
except ImportError:
    try:
        import TIPE_SimuOrbit_d as TIPE_SimuOrbit
        cpp_version = "Debug"
        pc.print_color("Using DEBUG version", pc.Color.Yellow)
    except ImportError:
        raise ImportError("No TIPE_SimuOrbit found!")

# === Constantes ===
mu = 1.32712440018e20
R_Mars = 2.28e11
r1 = np.array([1.5e11, 0.0, 0.0])

# === Paramètres réduits pour debug ===
n_points = 50  # TRÈS réduit pour debug
tof_min_days, tof_max_days = 50, 350
theta_min_deg, theta_max_deg = 0, 720  # Réduit aussi
tof_chosen_days, theta_chosen_deg = 250, 90

print(f"DEBUG: Grille réduite {n_points}x{n_points} = {n_points**2} points")

def compute_deltav_surface():
    """Version debug avec plus de diagnostics"""
    print("=== DEBUG: Calcul surface Δv ===")
    
    tof_list = np.linspace(tof_min_days*86400, tof_max_days*86400, n_points)
    theta_list = np.linspace(theta_min_deg, theta_max_deg, n_points)
    Theta, ToF = np.meshgrid(theta_list, tof_list/86400)
    
    print(f"DEBUG: Theta shape = {Theta.shape}, range = [{Theta.min():.1f}, {Theta.max():.1f}]")
    print(f"DEBUG: ToF shape = {ToF.shape}, range = [{ToF.min():.1f}, {ToF.max():.1f}]")
    
    r2_list = []
    tof_flat = []
    
    for i in range(n_points):
        for j in range(n_points):
            theta_rad = np.deg2rad(Theta[i, j])
            r2 = R_Mars * np.array([np.cos(theta_rad), np.sin(theta_rad), 0.0])
            r2_list.append(r2.tolist())
            tof_flat.append(tof_list[i])
    
    print(f"DEBUG: r2_list length = {len(r2_list)}")
    print(f"DEBUG: Premier r2 = {r2_list[0]}")
    print(f"DEBUG: Dernier r2 = {r2_list[-1]}")
    
    start_time = time.time()
    try:
        DeltaV_flat = TIPE_SimuOrbit.orbit.lambert_batch(
            r1.tolist(), r2_list, tof_flat, mu
        )
        print(f"DEBUG: lambert_batch OK, résultats = {len(DeltaV_flat)}")
    except Exception as e:
        print(f"DEBUG: Erreur lambert_batch: {e}")
        return None, None, None, 0, 0, 0
    
    end_time = time.time()
    computation_time = end_time - start_time
    
    DeltaV = np.array(DeltaV_flat).reshape(n_points, n_points)
    print(f"DEBUG: DeltaV shape = {DeltaV.shape}")
    print(f"DEBUG: DeltaV range = [{np.nanmin(DeltaV):.0f}, {np.nanmax(DeltaV):.0f}]")
    
    valid_computations = np.sum(np.isfinite(DeltaV))
    total_computations = DeltaV.size
    
    print(f"DEBUG: Solutions valides = {valid_computations}/{total_computations}")
    
    return Theta, ToF, DeltaV, computation_time, valid_computations, total_computations

# Calcul initial avec diagnostics
result = compute_deltav_surface()
if result[0] is None:
    print("ERREUR: Impossible de calculer la surface!")
    exit(1)

Theta, ToF, DeltaV, computation_time, valid_computations, total_computations = result

def compute_trajectory_debug(theta_deg, tof_days):
    """Version debug de compute_trajectory"""
    print(f"DEBUG: Calcul trajectoire θ={theta_deg}°, t={tof_days}j")
    
    theta_rad = np.deg2rad(theta_deg)
    r2 = R_Mars * np.array([np.cos(theta_rad), np.sin(theta_rad), 0.0])
    
    print(f"DEBUG: r1 = {r1}")
    print(f"DEBUG: r2 = {r2}")
    
    try:
        v1_chosen, v2_chosen = TIPE_SimuOrbit.orbit.lambert_universal(
            r1.tolist(), r2.tolist(), tof_days*86400, mu
        )
        v1_chosen = np.array(v1_chosen)
        v2_chosen = np.array(v2_chosen)
        
        print(f"DEBUG: v1 = {v1_chosen}")
        print(f"DEBUG: v2 = {v2_chosen}")
        
        if np.any(np.isnan(v1_chosen)) or np.any(np.isnan(v2_chosen)):
            print("DEBUG: Vitesses NaN!")
            return np.empty((0,3)), r2, np.nan
        
        # Trajectoire simplifiée
        n_traj_points = 50
        dt = tof_days * 86400 / n_traj_points
        r = r1.copy()
        v = v1_chosen.copy()
        traj = [r.copy()]
        
        for i in range(n_traj_points):
            a = -mu * r / np.linalg.norm(r)**3
            v = v + dt * a
            r = r + dt * v
            traj.append(r.copy())
            
        deltav = np.linalg.norm(v1_chosen) + np.linalg.norm(v2_chosen)
        
        print(f"DEBUG: Trajectoire OK, {len(traj)} points, Δv={deltav:.0f}")
        return np.array(traj), r2, deltav
        
    except Exception as e:
        print(f"DEBUG: Erreur trajectoire: {e}")
        return np.empty((0,3)), r2, np.nan

# === APPLICATION DASH SIMPLE ===
app = dash.Dash(__name__)

app.layout = html.Div([
    html.H1("🔍 DEBUG - Lambert Trajectory Explorer"),
    
    html.Div([
        html.P(f"Grille: {n_points}x{n_points} | Calcul: {computation_time:.2f}s | "
               f"Valides: {valid_computations}/{total_computations}"),
    ]),
    
    html.Div([
        html.Label("Temps (jours):"),
        dcc.Slider(
            id='time-slider',
            min=tof_min_days,
            max=tof_max_days,
            value=tof_chosen_days,
            marks={i: f'{i}' for i in range(tof_min_days, tof_max_days+1, 50)},
            tooltip={"placement": "bottom", "always_visible": True}
        )
    ]),
    
    html.Div([
        html.Label("Angle (deg):"),
        dcc.Slider(
            id='angle-slider',
            min=theta_min_deg,
            max=theta_max_deg,
            value=theta_chosen_deg,
            marks={i: f'{i}' for i in range(0, 181, 30)},
            tooltip={"placement": "bottom", "always_visible": True}
        )
    ]),
    
    html.Div([
        dcc.Graph(id='surface-plot', style={'height': '400px', 'width': '50%', 'display': 'inline-block'}),
        dcc.Graph(id='trajectory-plot', style={'height': '400px', 'width': '50%', 'display': 'inline-block'})
    ]),
    
    html.Div(id='debug-info')
])

@app.callback(
    [Output('surface-plot', 'figure'),
     Output('trajectory-plot', 'figure'),
     Output('debug-info', 'children')],
    [Input('time-slider', 'value'),
     Input('angle-slider', 'value')]
)
def update_debug(tof_days, theta_deg):
    print(f"\n=== CALLBACK: t={tof_days}, θ={theta_deg} ===")
    
    # Calcul trajectoire
    traj, r2, deltav = compute_trajectory_debug(theta_deg, tof_days)
    
    # === SURFACE 3D ===
    fig_surface = go.Figure()
    
    print(f"DEBUG: Ajout surface, DeltaV shape = {DeltaV.shape}")
    print(f"DEBUG: DeltaV sample values = {DeltaV[0:3, 0:3]}")
    
    # Vérifier que les données sont valides
    valid_mask = np.isfinite(DeltaV)
    print(f"DEBUG: {np.sum(valid_mask)} valeurs finies sur {DeltaV.size}")
    
    fig_surface.add_trace(go.Surface(
        x=Theta,
        y=ToF,
        z=DeltaV,
        colorscale='Viridis',
        name='Surface Δv'
    ))
    
    # Point sélectionné
    if np.isfinite(deltav):
        fig_surface.add_trace(go.Scatter3d(
            x=[theta_deg],
            y=[tof_days],
            z=[deltav],
            mode='markers',
            marker=dict(size=8, color='red'),
            name='Point sélectionné'
        ))
    
    fig_surface.update_layout(
        title='Surface Δv (DEBUG)',
        scene=dict(
            xaxis_title='Angle θ (deg)',
            yaxis_title='Temps (jours)',
            zaxis_title='Δv (m/s)'
        )
    )
    
    # === TRAJECTOIRE 3D ===
    fig_traj = go.Figure()
    
    # Orbite Mars
    theta_orbit = np.linspace(0, 2*np.pi, 50)
    mars_x = R_Mars * np.cos(theta_orbit)
    mars_y = R_Mars * np.sin(theta_orbit)
    mars_z = np.zeros_like(theta_orbit)
    
    fig_traj.add_trace(go.Scatter3d(
        x=mars_x,
        y=mars_y,
        z=mars_z,
        mode='lines',
        line=dict(color='orange', dash='dash'),
        name='Orbite Mars'
    ))
    
    # Trajectoire
    if traj.size > 0:
        print(f"DEBUG: Trajectoire shape = {traj.shape}")
        fig_traj.add_trace(go.Scatter3d(
            x=traj[:,0],
            y=traj[:,1],
            z=traj[:,2],
            mode='lines',
            line=dict(color='red'),
            name='Trajectoire'
        ))
    else:
        print("DEBUG: Pas de trajectoire à afficher")
    
    # Points
    fig_traj.add_trace(go.Scatter3d(
        x=[r1[0]], y=[r1[1]], z=[r1[2]],
        mode='markers',
        marker=dict(size=8, color='green'),
        name='Terre'
    ))
    
    fig_traj.add_trace(go.Scatter3d(
        x=[r2[0]], y=[r2[1]], z=[r2[2]],
        mode='markers',
        marker=dict(size=8, color='blue'),
        name='Mars'
    ))
    
    fig_traj.add_trace(go.Scatter3d(
        x=[0], y=[0], z=[0],
        mode='markers',
        marker=dict(size=12, color='yellow'),
        name='Soleil'
    ))
    
    fig_traj.update_layout(
        title='Trajectoire 3D (DEBUG)',
        scene=dict(
            xaxis_title='X (m)',
            yaxis_title='Y (m)',
            zaxis_title='Z (m)',
            aspectmode='cube'
        )
    )
    
    # Info debug
    debug_text = f"θ={theta_deg}°, t={tof_days}j, Δv={deltav:.0f} m/s, Traj points={len(traj) if traj.size > 0 else 0}"
    
    return fig_surface, fig_traj, debug_text

if __name__ == '__main__':
    print("\n=== LANCEMENT DEBUG ===")
    threading.Timer(1, lambda: webbrowser.open('http://127.0.0.1:8050/')).start()
    
    try:
        app.run(debug=True, host='127.0.0.1', port=8050)  # Mode debug activé
    except AttributeError:
        app.run_server(debug=True, host='127.0.0.1', port=8050)
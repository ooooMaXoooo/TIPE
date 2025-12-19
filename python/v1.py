import print_color as pc

Nb_modules = 8

import numpy as np
pc.print_color(f"numpy imported (1/{Nb_modules})", pc.Color.Green)
import time
pc.print_color(f"time imported (2/{Nb_modules})", pc.Color.Green)

import dash
from dash import dcc, html, callback, Input, Output, ctx
pc.print_color(f"dash imported (3/{Nb_modules})", pc.Color.Green)

import plotly.graph_objects as go
from plotly.subplots import make_subplots
pc.print_color(f"plotly imported (4/{Nb_modules})", pc.Color.Green)

import threading
pc.print_color(f"threading imported (5/{Nb_modules})", pc.Color.Green)
import webbrowser
pc.print_color(f"webbrowser imported (6/{Nb_modules})", pc.Color.Green)
import json
pc.print_color(f"json imported (7/{Nb_modules})", pc.Color.Green)

try:
    import TIPE_SimuOrbit  # version Release (rapide)
    pc.print_color(f"Using RELEASE version of TIPE_SimuOrbit ({Nb_modules}/{Nb_modules})", pc.Color.Green)
    cpp_version = "Release"
except ImportError:
    try:
        import TIPE_SimuOrbit_d as TIPE_SimuOrbit  # fallback vers Debug (pour dev)
        pc.print_color(f"WARNING: Using DEBUG version of TIPE_SimuOrbit ({Nb_modules}/{Nb_modules})", pc.Color.Yellow)
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

pc.print_color("All imports successful!", pc.Color.Green)

# === Constantes ===
mu = 1.32712440018e20
R_Mars = 2.28e11
r1 = np.array([1.5e11, 0.0, 0.0])

# === Paramètres utilisateur ===
n_points = 300  # Réduit pour de meilleures performances web
tof_min_days, tof_max_days = 30, 120
theta_min, theta_max = 11*np.pi/6, 13*np.pi/6
tof_chosen_days, theta_chosen = 100, 11.7*np.pi / 6

print("=== Paramètres utilisateur ===")
print(f"Nombre de points : {n_points}")
print(f"Intervalle temps (jours) : {tof_min_days} à {tof_max_days}")
print(f"Intervalle angle θ (deg) : {np.rad2deg(theta_min)} à {np.rad2deg(theta_max)}")
print(f"Temps choisi (jours) : {tof_chosen_days}")
print(f"Angle choisi θ (deg) : {np.rad2deg(theta_chosen)}\n")

# === Cache pour optimiser les performances ===
trajectory_cache = {}

def compute_deltav_surface():
    """Calcule la surface Δv avec le module C++"""
    print("=== Calcul de la surface Δv avec le module C++ (batch OpenMP) ===")
    
    tof_list = np.linspace(tof_min_days * 86400, tof_max_days * 86400, n_points)
    theta_list = np.linspace(theta_min, theta_max, n_points)
    Theta, ToF = np.meshgrid(theta_list, tof_list/86400)
    
    r2_list = []
    tof_flat = []
    
    for i in range(n_points):
        for j in range(n_points):
            theta = Theta[i, j]
            r2 = R_Mars * np.array([np.cos(theta), np.sin(theta), 0.0])
            r2_list.append(r2.tolist())
            tof_flat.append(tof_list[i])
    
    start_time = time.time()
    DeltaV_flat = TIPE_SimuOrbit.orbit.lambert_batch(
        r1.tolist(), r2_list, tof_flat, mu
    )
    end_time = time.time()
    computation_time = end_time - start_time
    
    DeltaV = np.array(DeltaV_flat).reshape(n_points, n_points)
    valid_computations = np.sum(np.isfinite(DeltaV))
    total_computations = DeltaV.size
    
    pc.print_color(f"Calcul terminé en {computation_time:.2f} secondes!", pc.Color.Green)
    pc.print_color(f"Vitesse: {total_computations/computation_time:.1f} trajectoire/seconde", pc.Color.Cyan)
    pc.print_color(f"Solutions valides: {valid_computations}/{total_computations} "
                   f"({100*valid_computations/total_computations:.1f}%)", pc.Color.Blue)
    
    return Theta, ToF, DeltaV, computation_time, valid_computations, total_computations

# Calcul initial
Theta, ToF, DeltaV, computation_time, valid_computations, total_computations = compute_deltav_surface()

def compute_trajectory(theta, tof_days):
    """Calcule une trajectoire spécifique avec cache"""
    cache_key = (round(theta, 5), round(tof_days, 1))
    
    if cache_key in trajectory_cache:
        return trajectory_cache[cache_key]
    
    r2 = R_Mars * np.array([np.cos(theta), np.sin(theta), 0.0])
    
    try:
        v1_chosen, v2_chosen = TIPE_SimuOrbit.orbit.lambert_universal(
            r1.tolist(), r2.tolist(), tof_days*86400, mu
        )
        v1_chosen = np.array(v1_chosen)
        v2_chosen = np.array(v2_chosen)
        
        if np.any(np.isnan(v1_chosen)) or np.any(np.isnan(v2_chosen)):
            result = (np.empty((0,3)), r2, np.nan)
            trajectory_cache[cache_key] = result
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
            
        deltav = np.linalg.norm(v1_chosen) + np.linalg.norm(v2_chosen)
        result = (np.array(traj), r2, deltav)
        trajectory_cache[cache_key] = result
        return result
        
    except Exception as e:
        result = (np.empty((0,3)), r2, np.nan)
        trajectory_cache[cache_key] = result
        return result

# === APPLICATION DASH ===
app = dash.Dash(__name__)

app.layout = html.Div([
    html.H1("🚀 Lambert Trajectory Explorer", 
            style={'textAlign': 'center', 'color': '#2c3e50', 'marginBottom': '20px'}),
    
    html.Div([
        html.Div([
            html.H4(f"⚡ Performance C++ {cpp_version}", style={'margin': 0, 'color': '#27ae60'}),
            html.P(f"Calcul: {computation_time:.2f}s | "
                   f"Vitesse: {total_computations/computation_time:.0f} calc/s | "
                   f"Solutions valides: {100*valid_computations/total_computations:.1f}%",
                   style={'margin': 0, 'color': '#7f8c8d'})
        ], style={'backgroundColor': '#ecf0f1', 'padding': '15px', 'borderRadius': '8px'})
    ], style={'marginBottom': '20px'}),
    
    html.Div([
        html.Div([
            html.Label("🕐 Temps de vol (jours):", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
            dcc.Slider(
                id='time-slider',
                min=tof_min_days,
                max=tof_max_days,
                value=tof_chosen_days,
                marks={i: f'{i}j' for i in range(tof_min_days, tof_max_days+1, 50)},
                tooltip={"placement": "bottom", "always_visible": True},
                updatemode='drag'
            )
        ], style={'width': '48%', 'display': 'inline-block', 'paddingRight': '2%'}),
        
        html.Div([
            html.Label("🎯 Angle θ (degrés):", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
            dcc.Slider(
                id='angle-slider',
                min=theta_min,
                max=theta_max,
                value=theta_chosen,
                marks={np.deg2rad(i): f"{i}°" for i in range(0, 361, 45)},
                tooltip={"placement": "bottom", "always_visible": True},
                updatemode='drag'
            )
        ], style={'width': '48%', 'display': 'inline-block', 'paddingLeft': '2%'})
    ], style={'backgroundColor': '#f8f9fa', 'padding': '20px', 'borderRadius': '8px', 'marginBottom': '20px'}),
    
    html.Div([
        html.Button("🎯 Minimum Δv", id='btn-min-deltav', n_clicks=0,
                   style={'margin': '5px', 'backgroundColor': '#3498db', 'color': 'white', 'border': 'none', 'padding': '8px 16px', 'borderRadius': '4px'}),
        html.Button("⚡ Transfert rapide", id='btn-fast-transfer', n_clicks=0,
                   style={'margin': '5px', 'backgroundColor': '#2ecc71', 'color': 'white', 'border': 'none', 'padding': '8px 16px', 'borderRadius': '4px'}),
        html.Button("🔄 Reset", id='btn-reset', n_clicks=0,
                   style={'margin': '5px', 'backgroundColor': '#e74c3c', 'color': 'white', 'border': 'none', 'padding': '8px 16px', 'borderRadius': '4px'})
    ], style={'textAlign': 'center', 'marginBottom': '20px'}),
    
    html.Div([
        html.Div(id='current-values', 
                style={'textAlign': 'center', 'fontSize': '18px', 'fontWeight': 'bold'})
    ], style={'marginBottom': '20px'}),
    
    html.Div([
        html.Div([
            dcc.Loading(
                dcc.Graph(id='deltav-surface', style={'height': '500px'}),
                type="default"
            )
        ], style={'width': '50%', 'display': 'inline-block'}),
        
        html.Div([
            dcc.Loading(
                dcc.Graph(id='trajectory-3d', style={'height': '500px'}),
                type="default"
            )
        ], style={'width': '50%', 'display': 'inline-block'})
    ]),
    
    html.Div([
        html.Label("Options d'affichage:", style={'fontWeight': 'bold'}),
        dcc.Checklist(
            id='display-options',
            options=[
                {'label': ' Orbite Mars', 'value': 'mars_orbit'},
                {'label': ' Point sélectionné', 'value': 'selected_point'},
                {'label': ' Échelle log Δv', 'value': 'log_scale'},
                {'label': ' Grille de référence', 'value': 'grid_lines'},
                {'label': ' Animation fluide', 'value': 'smooth_animation'}
            ],
            value=['mars_orbit', 'selected_point', 'smooth_animation'],
            inline=True,
            style={'marginLeft': '10px'}
        )
    ], style={'backgroundColor': '#f1f2f6', 'padding': '15px', 'borderRadius': '8px', 'marginTop': '20px'}),
    
    html.Div([
        html.Div(id='detailed-info',
                style={'fontSize': '14px', 'color': '#34495e'})
    ], style={'backgroundColor': '#ecf0f1', 'padding': '15px', 'borderRadius': '8px', 'marginTop': '20px'})

], style={'fontFamily': 'Arial, sans-serif', 'padding': '20px', 'backgroundColor': '#ffffff'})

# --- Callbacks ---
@app.callback(
    [Output('time-slider', 'value'),
     Output('angle-slider', 'value')],
    [Input('btn-min-deltav', 'n_clicks'),
     Input('btn-fast-transfer', 'n_clicks'),
     Input('btn-reset', 'n_clicks')]
)
def handle_buttons(min_clicks, fast_clicks, reset_clicks):
    if not ctx.triggered:
        return float(tof_chosen_days), float(theta_chosen)
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    if button_id == 'btn-min-deltav':
        min_idx = np.unravel_index(np.nanargmin(DeltaV), DeltaV.shape)
        return float(ToF[min_idx]), float(Theta[min_idx])
    elif button_id == 'btn-fast-transfer':
        return float(tof_min_days + 50), float(np.pi - 1e-6)
    elif button_id == 'btn-reset':
        return float(tof_chosen_days), float(theta_chosen)
    
    return float(tof_chosen_days), float(theta_chosen)

@app.callback(
    [Output('deltav-surface', 'figure'),
     Output('trajectory-3d', 'figure'),
     Output('current-values', 'children'),
     Output('detailed-info', 'children')],
    [Input('time-slider', 'value'),
     Input('angle-slider', 'value'),
     Input('display-options', 'value')]
)
def update_plots(tof_days, theta_deg, display_options):
    # --- Calcul trajectoire ---
    traj, r2, deltav_current = compute_trajectory(theta_deg, tof_days)
    
    from numpy import nan, log10, isfinite
    
    # --- Surface Δv ---
    DeltaV_masked = np.where(isfinite(DeltaV), DeltaV, nan)
    fig_surface = go.Figure()
    
    if 'log_scale' in display_options:
        DeltaV_log = np.where(DeltaV_masked > 0, log10(DeltaV_masked), nan)
        z_data = DeltaV_log
        zaxis_type = "linear"
        color_title = "log₁₀(Δv)"
    else:
        zaxis_type = "linear"
        z_data = DeltaV_masked
        color_title = "Δv (m/s)"
    
    customdata_str = np.where(np.isnan(DeltaV_masked), "n/a", DeltaV_masked.astype(int).astype(str))
    
    fig_surface.add_trace(go.Surface(
        x=Theta,
        y=ToF,
        z=z_data,
        customdata=customdata_str,
        colorscale='Viridis',
        name='Surface Δv',
        colorbar=dict(title=dict(text=color_title, side="right"), len=0.8),
        hovertemplate='Angle: %{x:.1f}°<br>Temps: %{y:.1f}j<br>Δv: %{customdata} m/s<extra></extra>'
    ))
    
    if 'selected_point' in display_options:
        fig_surface.add_trace(go.Scatter3d(
            x=[theta_deg],
            y=[tof_days],
            z=[deltav_current if 'log_scale' not in display_options else log10(deltav_current)],
            mode='markers',
            marker=dict(size=6, color='red'),
            name='Point sélectionné'
        ))
    
    fig_surface.update_layout(
        scene=dict(
            xaxis_title='θ (rad)',
            yaxis_title='ToF (jours)',
            zaxis_title=color_title
        ),
        margin=dict(l=0, r=0, t=30, b=0)
    )
    
    # --- Trajectoire 3D ---
    fig_traj = go.Figure()
    if traj.size > 0:
        fig_traj.add_trace(go.Scatter3d(
            x=traj[:,0], y=traj[:,1], z=traj[:,2],
            mode='lines',
            line=dict(color='blue', width=2),
            name='Trajectoire fusée'
        ))
    if 'mars_orbit' in display_options:
        theta_mars = np.linspace(0, 2*np.pi, 200)
        mars_orbit = np.array([R_Mars*np.cos(theta_mars), R_Mars*np.sin(theta_mars), np.zeros_like(theta_mars)]).T
        fig_traj.add_trace(go.Scatter3d(
            x=mars_orbit[:,0], y=mars_orbit[:,1], z=mars_orbit[:,2],
            mode='lines',
            line=dict(color='orange', width=2, dash='dash'),
            name='Orbite Mars'
        ))
    
    fig_traj.update_layout(
        scene=dict(
            xaxis_title='X (m)',
            yaxis_title='Y (m)',
            zaxis_title='Z (m)'
        ),
        margin=dict(l=0, r=0, t=30, b=0)
    )
    
    current_values_str = f"Δv = {deltav_current:.2f} m/s | θ = {np.rad2deg(theta_deg):.2f}° | ToF = {tof_days:.1f} jours"
    
    detailed_info_str = f"Point cible Mars: {r2} m\nTrajectoire points: {traj.shape[0]}"
    
    return fig_surface, fig_traj, current_values_str, detailed_info_str

# --- Ouverture navigateur automatique ---
def open_browser():
    time.sleep(2)
    webbrowser.open('http://127.0.0.1:8050/')

if __name__ == '__main__':
    threading.Timer(2, open_browser).start()
    try:
        app.run(debug=False, host='127.0.0.1', port=8050)
    except AttributeError:
        app.run_server(debug=False, host='127.0.0.1', port=8050)

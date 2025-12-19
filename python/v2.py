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
mu = 1.32712440018e20  # du soleil
#R_Mars = 2.9e12 #2.28e11
#R_Mars = 7.78e11
R_Mars = 2.28e11
r1 = np.array([1.5e11, 0.0, 0.0])


# === fonctions ===
def clamp (x, a, b) :
    min = a if a < b else b
    max = a+b-min

    x = x if x > a else a
    x = x if x < b else b

    return x


# === Paramètres utilisateur ===
n_points = 100
tof_min_days, tof_max_days = 50, 500
theta_min, theta_max = 0.1, 2*np.pi
tof_chosen_days, theta_chosen = 259, np.deg2rad(179.99)

tof_chosen_days = clamp(tof_chosen_days, tof_min_days, tof_max_days)
theta_chosen    = clamp(theta_chosen, theta_min, theta_max)


print("=== Paramètres utilisateur ===")
print(f"Nombre de points : {n_points}")
print(f"Intervalle temps (jours) : {tof_min_days} à {tof_max_days}")
print(f"Intervalle angle θ (deg) : {np.rad2deg(theta_min)} à {np.rad2deg(theta_max)}")
print(f"Temps choisi (jours) : {tof_chosen_days}")
print(f"Angle choisi θ (deg) : {np.rad2deg(theta_chosen)}\n")

# === Cache pour optimiser les performances ===
trajectory_cache = {}

def compute_deltav_surface():
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

Theta, ToF, DeltaV, computation_time, valid_computations, total_computations = compute_deltav_surface()

def compute_trajectory(theta, tof_days):
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
                min=np.rad2deg(theta_min),
                max=np.rad2deg(theta_max),
                value=np.rad2deg(theta_chosen),
                marks={i: f"{i}°" for i in range(0, 361, 45)},
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

# === CALLBACKS ===
@app.callback(
    [Output('time-slider', 'value'),
     Output('angle-slider', 'value')],
    [Input('btn-min-deltav', 'n_clicks'),
     Input('btn-fast-transfer', 'n_clicks'),
     Input('btn-reset', 'n_clicks')]
)
def handle_buttons(min_clicks, fast_clicks, reset_clicks):
    if not ctx.triggered:
        return tof_chosen_days, np.rad2deg(theta_chosen)
    
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    if button_id == 'btn-min-deltav':
        min_idx = np.unravel_index(np.nanargmin(DeltaV), DeltaV.shape)
        return ToF[min_idx], np.rad2deg(Theta[min_idx])
    elif button_id == 'btn-fast-transfer':
        return tof_min_days + 50, 180
    elif button_id == 'btn-reset':
        return tof_chosen_days, np.rad2deg(theta_chosen)
    
    return tof_chosen_days, np.rad2deg(theta_chosen)

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
    theta_rad = np.deg2rad(theta_deg)
    traj, r2, deltav_current = compute_trajectory(theta_rad, tof_days)

    # --- Surface Δv ---
    from numpy import nan, log10, isfinite
    DeltaV_masked = np.where(isfinite(DeltaV), DeltaV, nan)
    if 'log_scale' in display_options:
        z_data = np.where(DeltaV_masked>0, log10(DeltaV_masked), nan)
        zaxis_type = "linear"
        color_title = "log₁₀(Δv)"
    else:
        z_data = DeltaV_masked
        zaxis_type = "linear"
        color_title = "Δv (m/s)"
    customdata_str = np.where(np.isnan(DeltaV_masked), "n/a", DeltaV_masked.astype(int).astype(str))

    fig_surface = go.Figure()
    fig_surface.add_trace(go.Surface(
        x=Theta, y=ToF, z=z_data, customdata=customdata_str, colorscale='Viridis',
        name='Surface Δv',
        colorbar=dict(title=dict(text=color_title, side="right"), len=0.8),
        hovertemplate='Angle: %{x:.1f}°<br>Temps: %{y:.1f}j<br>Δv: %{customdata} m/s<extra></extra>'
    ))
    if 'selected_point' in display_options and np.isfinite(deltav_current):
        fig_surface.add_trace(go.Scatter3d(
            x=[theta_rad], y=[tof_days], z=[log10(deltav_current) if 'log_scale' in display_options else deltav_current],
            mode='markers', marker=dict(size=10, color='red', symbol='circle'),
            name=f'Sélection: Δv={deltav_current:.0f} m/s'
        ))
    fig_surface.update_layout(
        title=f'Surface Δv(θ, temps) - C++ {cpp_version}<br>'
              f'<sub>{computation_time:.2f}s, {total_computations/computation_time:.0f} calc/s</sub>',
        scene=dict(
            xaxis_title='Angle θ (rad)',
            yaxis_title='Temps (jours)',
            zaxis_title=color_title,
            zaxis_type=zaxis_type,
            camera=dict(eye=dict(x=1.2, y=1.2, z=0.8))
        ),
        height=500
    )

    # --- Trajectoire 3D ---
    fig_traj = go.Figure()
    if 'mars_orbit' in display_options:
        theta_mars = np.linspace(0, 2*np.pi, 200)
        mars_orbit = np.array([R_Mars*np.cos(theta_mars), R_Mars*np.sin(theta_mars), np.zeros_like(theta_mars)]).T
        fig_traj.add_trace(go.Scatter3d(
            x=mars_orbit[:,0], y=mars_orbit[:,1], z=mars_orbit[:,2],
            mode='lines', line=dict(color='orange', width=2, dash='dash'), name='Orbite Mars'
        ))
    if traj.size>0:
        fig_traj.add_trace(go.Scatter3d(
            x=traj[:,0], y=traj[:,1], z=traj[:,2],
            mode='lines', line=dict(color='blue', width=2), name='Trajectoire fusée'
        ))
    # Soleil, départ, arrivée
    fig_traj.add_trace(go.Scatter3d(x=[0], y=[0], z=[0], mode='markers', marker=dict(size=6, color='yellow'), name='Soleil'))
    fig_traj.add_trace(go.Scatter3d(x=[r1[0]], y=[r1[1]], z=[r1[2]], mode='markers', marker=dict(size=5, color='green'), name='Départ'))
    fig_traj.add_trace(go.Scatter3d(x=[r2[0]], y=[r2[1]], z=[r2[2]], mode='markers', marker=dict(size=5, color='red'), name='Arrivée'))
    fig_traj.update_layout(scene=dict(
        xaxis_title='X (m)', yaxis_title='Y (m)', zaxis_title='Z (m)',
        camera=dict(eye=dict(x=1.5, y=1.5, z=1))
    ), height=500)

    current_values = f"Δv actuel: {deltav_current:.2f} m/s | θ={theta_deg:.1f}°, ToF={tof_days:.1f} jours"
    detailed_info = f"Point final r2: {r2}"

    return fig_surface, fig_traj, current_values, detailed_info

# === RUN SERVER ===
def open_browser():
    webbrowser.open_new("http://127.0.0.1:8050/")

threading.Timer(1.0, open_browser).start()
# Lancer l'app
# Version qui fonctionne avec toutes les versions de Dash
try:
    app.run(debug=False, host='127.0.0.1', port=8050)
except AttributeError:
    app.run_server(debug=False, host='127.0.0.1', port=8050)

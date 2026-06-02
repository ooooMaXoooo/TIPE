import numpy as np
import matplotlib.pyplot as plt

import cmath
import ImportData

try:
    import TIPE_SimuOrbit as cpp  # version Release (rapide)
    cpp_version = "Release"
except ImportError:
    try:
        import TIPE_SimuOrbit_d as cpp  # fallback vers Debug (pour dev)
        cpp_version = "Debug"
    except ImportError:
        raise ImportError("\x1B[38;5;202m Neither TIPE_SimuOrbit.pyd nor TIPE_SimuOrbit_d.pyd found \033[0m")

from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading


dissimilarity_threshold = 100  # à modifier à la main pour coller au c++
family_gap              = 0.4  # espace vertical (en individus) entre deux familles dans une barre
outlier_iqr_factor      = 1.5  # facteur IQR pour la détection des outliers (boîte à moustaches)
outlier_dark_factor     = 0.30 # facteur d'assombrissement des couleurs outliers (0 = noir, 1 = normal)


# ── Cache ──────────────────────────────────────────────────────────────────────

_rocket_cache = {}
_cache_lock   = threading.Lock()


def _load_rocket_data(dossier, gen, idx):
    key = (dossier, gen, idx)
    with _cache_lock:
        if key in _rocket_cache:
            return _rocket_cache[key]

    filepath = dossier + "/RocketsData/gen_" + str(gen) + "/ind_" + str(idx) + ".rck"
    (
        _,
        _, _,
        thetaMin, thetaMax, thetaMean,
        NbTours,
        RMin, RMax, RMean,
        _, _,
        _, _,
        _,
        _, _, _,
        nbImpulsions,
        _
    ) = ImportData.importRocketData(filepath)

    data = dict(
        rMin=RMin, rMax=RMax, rMean=RMean,
        thetaMin=thetaMin, thetaMax=thetaMax, thetaMean=thetaMean,
        nbImpuls=nbImpulsions, nbTurns=NbTours
    )
    with _cache_lock:
        _rocket_cache[key] = data
    return data


def preload_generation(dossier, gen, all_clusters):
    indices = [idx for cluster in all_clusters for idx in cluster]
    with ThreadPoolExecutor() as ex:
        futures = {ex.submit(_load_rocket_data, dossier, gen, idx): idx for idx in indices}
        for f in as_completed(futures):
            f.result()


def _distance_from_cache(dossier, gen1, idx1, gen2, idx2):
    d1 = _load_rocket_data(dossier, gen1, idx1)
    d2 = _load_rocket_data(dossier, gen2, idx2)
    return cpp.genetic.distance_rocket(
        rMin1=d1["rMin"],         rMax1=d1["rMax"],         rMean1=d1["rMean"],
        thetaMin1=d1["thetaMin"], thetaMax1=d1["thetaMax"], thetaMean1=d1["thetaMean"],
        nbImpuls1=d1["nbImpuls"], nbTurns1=d1["nbTurns"],
        rMin2=d2["rMin"],         rMax2=d2["rMax"],         rMean2=d2["rMean"],
        thetaMin2=d2["thetaMin"], thetaMax2=d2["thetaMax"], thetaMean2=d2["thetaMean"],
        nbImpuls2=d2["nbImpuls"], nbTurns2=d2["nbTurns"]
    )


# ── Palette basée sur les longueurs d'onde ─────────────────────────────────────

def wavelength_to_rgb(wavelength, gamma=0.8):
    """
    Convertit une longueur d'onde (nm, 380–700) en couleur RGB dans [0, 1].
    Algorithme de Dan Bruton.
    """
    w = float(wavelength)
    if   380 <= w < 440:  r, g, b = -(w - 440) / (440 - 380), 0.0, 1.0
    elif 440 <= w < 490:  r, g, b = 0.0, (w - 440) / (490 - 440), 1.0
    elif 490 <= w < 510:  r, g, b = 0.0, 1.0, -(w - 510) / (510 - 490)
    elif 510 <= w < 580:  r, g, b = (w - 510) / (580 - 510), 1.0, 0.0
    elif 580 <= w < 645:  r, g, b = 1.0, -(w - 645) / (645 - 580), 0.0
    elif 645 <= w <= 700: r, g, b = 1.0, 0.0, 0.0
    else:                 r, g, b = 0.0, 0.0, 0.0

    if   380 <= w < 420: factor = 0.3 + 0.7 * (w - 380) / (420 - 380)
    elif 680 < w <= 700: factor = 0.3 + 0.7 * (700 - w) / (700 - 680)
    else:                factor = 1.0

    return ((r * factor) ** gamma,
            (g * factor) ** gamma,
            (b * factor) ** gamma)


def classical_mds_1d(D):
    """MDS classique 1D à partir d'une matrice de distances D (n×n)."""
    n  = D.shape[0]
    D2 = D ** 2
    H  = np.eye(n) - np.ones((n, n)) / n
    B  = -0.5 * H @ D2 @ H
    B  = (B + B.T) / 2
    eigenvalues, eigenvectors = np.linalg.eigh(B)
    idx   = np.argmax(eigenvalues)
    coord = eigenvectors[:, idx] * np.sqrt(max(eigenvalues[idx], 0.0))
    return coord


def compute_family_palette(all_color_ids, inter_family_distances,
                           iqr_factor=1.5, dark_factor=0.30):
    """
    Attribue une couleur à chaque famille via MDS 1D + longueurs d'onde.

    Les familles outliers (éloignées du groupe principal selon la règle IQR)
    reçoivent une version assombrie de leur couleur naturelle.
    Les familles majoritaires voient leurs positions renormalisées sur [0, 1]
    pour maximiser la diversité de couleurs perçues.

    Retourne (palette_dict, positions_dict) où positions_dict sert au tri.
    """
    n = len(all_color_ids)

    if n == 1:
        return {all_color_ids[0]: wavelength_to_rgb(550)}, {all_color_ids[0]: 0.5}

    id_to_local = {cid: i for i, cid in enumerate(all_color_ids)}

    # ── Matrice de distances ────────────────────────────────────────────────
    D = np.full((n, n), np.inf)
    np.fill_diagonal(D, 0.0)

    for (cid_a, cid_b), dist in inter_family_distances.items():
        i, j = id_to_local[cid_a], id_to_local[cid_b]
        D[i, j] = dist
        D[j, i] = dist

    # Floyd-Warshall
    for k in range(n):
        for i in range(n):
            for j in range(n):
                candidate = D[i, k] + D[k, j]
                if candidate < D[i, j]:
                    D[i, j] = candidate

    finite_vals = D[np.isfinite(D)]
    max_finite  = float(np.max(finite_vals)) if finite_vals.size > 0 else 1.0
    D           = np.where(np.isinf(D), max_finite * 2.0, D)

    # ── MDS classique 1D ────────────────────────────────────────────────────
    raw_positions = classical_mds_1d(D)

    # Normalisation globale en [0, 1] (sert de référence pour le tri)
    lo, hi        = raw_positions.min(), raw_positions.max()
    positions     = (raw_positions - lo) / (hi - lo) if hi > lo else np.full(n, 0.5)

    # ── Détection des outliers par règle IQR (boîte à moustaches) ──────────
    Q1, Q3    = np.percentile(positions, [25, 75])
    IQR       = Q3 - Q1
    low_fence = Q1 - iqr_factor * IQR
    high_fence= Q3 + iqr_factor * IQR
    is_outlier= (positions < low_fence) | (positions > high_fence)

    nb_outliers = int(np.sum(is_outlier))
    if nb_outliers > 0:
        print(f"  → {nb_outliers} famille(s) outlier(s) détectée(s) sur {n} "
              f"(IQR×{iqr_factor} : [{low_fence:.3f}, {high_fence:.3f}])")

    # ── Renormalisation des majoritaires sur [0, 1] ─────────────────────────
    # Maximise la diversité de couleurs pour les familles qui comptent
    majority_pos = positions[~is_outlier]
    if majority_pos.size > 1:
        maj_lo, maj_hi = majority_pos.min(), majority_pos.max()
        if maj_hi > maj_lo:
            renorm = np.where(
                is_outlier,
                positions,                                      # outliers : position brute
                (positions - maj_lo) / (maj_hi - maj_lo)       # majoritaires : renormalisés
            )
        else:
            renorm = positions
    else:
        renorm = positions

    # ── Conversion en longueurs d'onde puis RGB ─────────────────────────────
    WL_MIN, WL_MAX = 420.0, 680.0

    palette_dict  = {}
    positions_dict = {}

    for cid in all_color_ids:
        local_idx = id_to_local[cid]
        pos       = renorm[local_idx]
        rgb       = wavelength_to_rgb(WL_MIN + pos * (WL_MAX - WL_MIN))

        if is_outlier[local_idx]:
            rgb = tuple(c * dark_factor for c in rgb)

        palette_dict[cid]   = rgb
        positions_dict[cid] = float(positions[local_idx])  # tri sur la position brute

    return palette_dict, positions_dict


# ── Propagation ────────────────────────────────────────────────────────────────

def _compute_row(args):
    i, cluster_i, clusters_n, dossier, gen_n1, gen_n = args
    row = np.empty(len(clusters_n))
    for j, cluster_j in enumerate(clusters_n):
        total = sum(
            _distance_from_cache(dossier, gen_n1, idx_i, gen_n, idx_j)
            for idx_i in cluster_i
            for idx_j in cluster_j
        )
        row[j] = total / (len(cluster_i) * len(cluster_j))
    return i, row


def propagate_colors(clusters_n, clusters_n1, dossier, gen_n, gen_n1,
                     threshold, next_color, color_map_n, pbar=None):
    """
    Retourne (color_map_n1, next_color, inter_family_distances).
    inter_family_distances : dict {(cid_a, cid_b): distance_minimale_observée}
    """
    Cn  = len(clusters_n)
    Cn1 = len(clusters_n1)

    D = np.full((Cn1, Cn), np.inf)

    args_list = [
        (i, clusters_n1[i], clusters_n, dossier, gen_n1, gen_n)
        for i in range(Cn1)
    ]

    with ThreadPoolExecutor() as ex:
        futures = {ex.submit(_compute_row, args): args[0] for args in args_list}
        for future in as_completed(futures):
            i, row = future.result()
            D[i] = row
            if pbar is not None:
                nb_calls = sum(len(clusters_n1[i]) * len(cj) for cj in clusters_n)
                pbar.update(nb_calls)

    pairs = sorted(
        [(D[i, j], i, j) for i in range(Cn1) for j in range(Cn)],
        key=lambda x: x[0]
    )

    color_map_n1 = {}
    matched_n1   = set()
    matched_n    = set()

    for dist, i, j in pairs:
        if dist > threshold:
            break
        if i in matched_n1 or j in matched_n:
            continue
        color_map_n1[i] = color_map_n[j]
        matched_n1.add(i)
        matched_n.add(j)

    for i in range(Cn1):
        if i not in matched_n1:
            color_map_n1[i] = next_color
            next_color += 1

    # Distances inter-familles depuis la matrice D
    inter_family_distances = {}
    for i in range(Cn1):
        for j in range(Cn):
            cid_i = color_map_n1[i]
            cid_j = color_map_n[j]
            if cid_i == cid_j:
                continue
            key = (min(cid_i, cid_j), max(cid_i, cid_j))
            if key not in inter_family_distances or D[i, j] < inter_family_distances[key]:
                inter_family_distances[key] = D[i, j]

    return color_map_n1, next_color, inter_family_distances


# ── Histogramme HAC ────────────────────────────────────────────────────────────

def histo_clusters(fig, ax, gen_final, dossier, gen_start=1, gen_step=1):
    """
    Affiche l'histogramme des clusters pour les générations sélectionnées.

    gen_final  : dernière génération à afficher
    dossier    : répertoire de la simulation
    gen_start  : première génération à afficher (1-indexé, défaut=1)
    gen_step   : pas entre deux générations affichées (défaut=1)
    """
    filename_base = dossier + "/Stats/HAC/gen_"

    selected_gens = list(range(gen_start, gen_final + 1, gen_step))
    Nb_selected   = len(selected_gens)

    AllClusters   = []
    ClustersSizes = []

    for gen in selected_gens:
        (N, clusters) = ImportData.importClusters(filename_base + str(gen) + ".cluster")
        AllClusters.append(clusters)
        ClustersSizes.append(N)

    print("Chargement des données...")
    for i, gen in enumerate(tqdm(selected_gens, desc="Lecture fichiers", unit="gen")):
        preload_generation(dossier, gen, AllClusters[i])

    total_calls = sum(
        len(ci) * len(cj)
        for i in range(1, Nb_selected)
        for ci in AllClusters[i]
        for cj in AllClusters[i - 1]
    )

    color_maps             = []
    next_color             = 0
    all_inter_family_dists = {}

    initial_map = {j: j for j in range(ClustersSizes[0])}
    next_color  = ClustersSizes[0]
    color_maps.append(initial_map)

    with tqdm(total=total_calls, desc="Propagation des couleurs", unit="dist") as pbar:
        for i in range(1, Nb_selected):
            color_map, next_color, new_dists = propagate_colors(
                clusters_n  = AllClusters[i - 1],
                clusters_n1 = AllClusters[i],
                dossier     = dossier,
                gen_n       = selected_gens[i - 1],
                gen_n1      = selected_gens[i],
                threshold   = dissimilarity_threshold,
                next_color  = next_color,
                color_map_n = color_maps[i - 1],
                pbar        = pbar
            )
            color_maps.append(color_map)
            for key, dist in new_dists.items():
                if key not in all_inter_family_dists or dist < all_inter_family_dists[key]:
                    all_inter_family_dists[key] = dist

    all_color_ids = sorted(set(cid for cm in color_maps for cid in cm.values()))
    print(f"Nombre total de familles distinctes : {len(all_color_ids)}")

    palette_dict, positions_dict = compute_family_palette(
        all_color_ids, all_inter_family_dists,
        iqr_factor  = outlier_iqr_factor,
        dark_factor = outlier_dark_factor
    )

    # Affichage — clusters triés par position MDS croissante
    bar_width = gen_step * 0.9

    for i, gen in enumerate(selected_gens):
        bottom = 0
        sorted_clusters = sorted(
            range(ClustersSizes[i]),
            key=lambda j: positions_dict[color_maps[i][j]]
        )
        for k, cluster_j in enumerate(sorted_clusters):
            if k > 0:
                ax.bar(gen, family_gap, bottom=bottom,
                       color="white", width=bar_width, edgecolor="none")
                bottom += family_gap

            size  = len(AllClusters[i][cluster_j])
            color = palette_dict[color_maps[i][cluster_j]]
            ax.bar(gen, size, bottom=bottom, color=color,
                   width=bar_width, edgecolor="black", linewidth=0.3)
            bottom += size

    ax.set_xlabel("Génération")
    ax.set_ylabel("Taille des clusters")
    ax.set_title("Évolution des clusters par génération (HAC — average linkage)")
    plt.tight_layout()




# ── Graphe score ────────────────────────────────────────────────────────────

def plotScores (fig, ax, dossier, end_gen, start_gen=1, gen_step=1) :
    filename_base = dossier + "Stats/Simple/gen_"

    
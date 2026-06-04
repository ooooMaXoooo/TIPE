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
family_gap              = 0    # espace vertical (en individus) entre deux familles dans une barre
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


# ── Barycentres ────────────────────────────────────────────────────────────────

_FIELDS = ["rMin", "rMax", "rMean", "thetaMin", "thetaMax", "thetaMean", "nbImpuls", "nbTurns"]

def compute_cluster_barycenter(dossier, gen, cluster):
    """
    Calcule le barycentre d'un cluster : moyenne des features de chaque individu.
    Toutes les données doivent être déjà dans le cache (appeler preload_generation avant).
    """
    sums = {f: 0.0 for f in _FIELDS}
    for idx in cluster:
        data = _load_rocket_data(dossier, gen, idx)
        for f in _FIELDS:
            sums[f] += data[f]
    n = len(cluster)
    return {f: sums[f] / n for f in _FIELDS}


def _distance_barycenters(b1, b2):
    """Distance entre deux barycentres via la fonction C++."""
    return cpp.genetic.distance_rocket(
        rMin1=b1["rMin"],         rMax1=b1["rMax"],         rMean1=b1["rMean"],
        thetaMin1=b1["thetaMin"], thetaMax1=b1["thetaMax"], thetaMean1=b1["thetaMean"],
        nbImpuls1=int(b1["nbImpuls"]), nbTurns1=b1["nbTurns"],
        rMin2=b2["rMin"],         rMax2=b2["rMax"],         rMean2=b2["rMean"],
        thetaMin2=b2["thetaMin"], thetaMax2=b2["thetaMax"], thetaMean2=b2["thetaMean"],
        nbImpuls2=int(b2["nbImpuls"]), nbTurns2=b2["nbTurns"]
    )








def get_cluster_representative(dossier, gen, cluster):
    """
    Retourne l'indice de l'individu le plus proche du barycentre du cluster.
    Si le cluster est un singleton, retourne directement son unique membre.
    """
    if len(cluster) == 1:
        return cluster[0]

    bary = compute_cluster_barycenter(dossier, gen, cluster)

    best_idx  = cluster[0]
    best_dist = np.inf
    for idx in cluster:
        data = _load_rocket_data(dossier, gen, idx)
        dist = _distance_barycenters(data, bary)
        if dist < best_dist:
            best_dist = dist
            best_idx  = idx

    return best_idx


def get_generation_representatives(dossier, gen_idx, gen_number,
                                   clusters, color_map, palette_dict):
    """
    Pour une génération donnée, retourne la liste des représentants de chaque cluster
    avec leur couleur associée.

    Paramètres
    ----------
    dossier      : répertoire de la simulation
    gen_idx      : index dans AllClusters / color_maps (0-based)
    gen_number   : numéro de génération réel (pour les chemins de fichiers)
    clusters     : AllClusters[gen_idx]
    color_map    : color_maps[gen_idx]
    palette_dict : dict color_id -> (r, g, b)

    Retourne
    --------
    Liste de dict :
        { "idx": int, "color": (r,g,b), "color_id": int, "cluster_size": int }
    """
    result = []
    for cluster_j, cluster in enumerate(clusters):
        color_id  = color_map[cluster_j]
        color     = palette_dict[color_id]
        rep_idx   = get_cluster_representative(dossier, gen_number, cluster)
        result.append({
            "idx":          rep_idx,
            "color":        color,
            "color_id":     color_id,
            "cluster_size": len(cluster)
        })
    return result







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
    Retourne (palette_dict, positions_dict).
    """
    n = len(all_color_ids)
    if n == 1:
        return {all_color_ids[0]: wavelength_to_rgb(550)}, {all_color_ids[0]: 0.5}

    id_to_local = {cid: i for i, cid in enumerate(all_color_ids)}

    D = np.full((n, n), np.inf)
    np.fill_diagonal(D, 0.0)
    for (cid_a, cid_b), dist in inter_family_distances.items():
        i, j = id_to_local[cid_a], id_to_local[cid_b]
        D[i, j] = dist
        D[j, i] = dist

    for k in tqdm(range(n), desc="Floyd-Warshall", unit="iter", leave=False):
        D = np.minimum(D, D[:, k:k+1] + D[k:k+1, :])

    finite_vals = D[np.isfinite(D)]
    max_finite  = float(np.max(finite_vals)) if finite_vals.size > 0 else 1.0
    D           = np.where(np.isinf(D), max_finite * 2.0, D)

    raw_positions = classical_mds_1d(D)
    lo, hi        = raw_positions.min(), raw_positions.max()
    positions     = (raw_positions - lo) / (hi - lo) if hi > lo else np.full(n, 0.5)

    Q1, Q3     = np.percentile(positions, [25, 75])
    IQR        = Q3 - Q1
    low_fence  = Q1 - iqr_factor * IQR
    high_fence = Q3 + iqr_factor * IQR
    is_outlier = (positions < low_fence) | (positions > high_fence)

    nb_outliers = int(np.sum(is_outlier))
    if nb_outliers > 0:
        print(f"  → {nb_outliers} famille(s) outlier(s) sur {n} "
              f"(IQR×{iqr_factor} : [{low_fence:.3f}, {high_fence:.3f}])")

    majority_pos = positions[~is_outlier]
    if majority_pos.size > 1:
        maj_lo, maj_hi = majority_pos.min(), majority_pos.max()
        if maj_hi > maj_lo:
            renorm = np.where(
                is_outlier,
                positions,
                (positions - maj_lo) / (maj_hi - maj_lo)
            )
        else:
            renorm = positions
    else:
        renorm = positions

    WL_MIN, WL_MAX = 420.0, 680.0
    palette_dict   = {}
    positions_dict = {}
    for cid in all_color_ids:
        local_idx = id_to_local[cid]
        pos       = renorm[local_idx]
        rgb       = wavelength_to_rgb(WL_MIN + pos * (WL_MAX - WL_MIN))
        if is_outlier[local_idx]:
            rgb = tuple(c * dark_factor for c in rgb)
        palette_dict[cid]   = rgb
        positions_dict[cid] = float(positions[local_idx])

    return palette_dict, positions_dict

# ── Propagation (centroid linkage) ─────────────────────────────────────────────

def propagate_colors(clusters_n, clusters_n1, dossier, gen_n, gen_n1,
                     threshold, next_color, color_map_n, pbar=None):
    """
    Propage les couleurs de gen_n vers gen_n1 par matching glouton sur la
    distance entre barycentres (centroid linkage).

    Retourne (color_map_n1, next_color, inter_family_distances).
    """
    Cn  = len(clusters_n)
    Cn1 = len(clusters_n1)

    # Calcul des barycentres — O(|cluster| × nb_features) par cluster
    bary_n  = [compute_cluster_barycenter(dossier, gen_n,  c) for c in clusters_n]
    bary_n1 = [compute_cluster_barycenter(dossier, gen_n1, c) for c in clusters_n1]

    # Matrice de distances barycentre-à-barycentre — O(Cn × Cn1) appels C++
    D = np.array([
        [_distance_barycenters(bary_n1[i], bary_n[j]) for j in range(Cn)]
        for i in range(Cn1)
    ])

    if pbar is not None:
        # On compte Cn × Cn1 "distances" calculées (une par paire de barycentres)
        pbar.update(Cn * Cn1)

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

    # Distances inter-familles pour le MDS global
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

def chargeData_histo_clusters(dossier, gen_start, gen_final, gen_step):
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

    # Nombre de paires de barycentres calculées (une par paire de clusters)
    total_calls = sum(
        len(AllClusters[i]) * len(AllClusters[i - 1])
        for i in range(1, Nb_selected)
    )

    color_maps             = []
    next_color             = 0
    all_inter_family_dists = {}

    initial_map = {j: j for j in range(ClustersSizes[0])}
    next_color  = ClustersSizes[0]
    color_maps.append(initial_map)

    with tqdm(total=total_calls, desc="Propagation des couleurs", unit="paires") as pbar:
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

    return (selected_gens, ClustersSizes, positions_dict, color_maps, AllClusters, palette_dict)


def affiche_histo_clusters(fig, ax, gen_step, selected_gens, ClustersSizes,
                           positions_dict, color_maps, AllClusters, palette_dict):
    bar_width = gen_step

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

            size = len(AllClusters[i][cluster_j])
            color = palette_dict[color_maps[i][cluster_j]]
            ax.bar(gen, size, bottom=bottom, color=color,
                   width=bar_width, edgecolor="black", linewidth=0)
            bottom += size

    ax.set_xlabel("Génération")
    ax.set_ylabel("Taille des clusters")
    ax.set_title("Évolution des clusters par génération (HAC — average linkage)")
    plt.tight_layout()


def histo_clusters(fig, ax, gen_final, dossier, gen_start=1, gen_step=1):
    (selected_gens, ClustersSizes, positions_dict,
     color_maps, AllClusters, palette_dict) = chargeData_histo_clusters(
        dossier, gen_start, gen_final, gen_step)
    affiche_histo_clusters(fig, ax, gen_step, selected_gens, ClustersSizes,
                           positions_dict, color_maps, AllClusters, palette_dict)







# ── Streamgraph HAC ────────────────────────────────────────────────────────────

def _distinct_palette(all_ids, positions_dict):
    """
    Une couleur DISTINCTE par famille (garantit autant de couleurs que de
    familles), espacees par rang le long de 'turbo'. L'ordre suit la position
    MDS, donc familles proches -> teintes proches (continuite de l'aspect), et
    la couleur reste attachee a l'id de famille -> continuite au fil des gen.
    """
    order = sorted(all_ids, key=lambda fid: positions_dict.get(fid, 0.5))
    n = max(len(order), 1)
    cmap = plt.get_cmap("turbo")
    return {fid: cmap((k + 0.5) / n) for k, fid in enumerate(order)}



def affiche_streamgraph_clusters(fig, ax, selected_gens, color_maps, AllClusters,
                                 positions_dict, palette_dict=None, baseline="sym", verbose=False):
    """Dessine le streamgraph et verifie la coherence familles <-> couleurs."""
    # familles ordonnees par position MDS : familles proches -> bandes adjacentes
    all_ids = sorted({cid for cm in color_maps for cid in cm.values()},
                     key=lambda fid: positions_dict.get(fid, 0.5))
    if palette_dict is None:
        palette_dict = _distinct_palette(all_ids, positions_dict)
    row = {fid: k for k, fid in enumerate(all_ids)}

    # matrice effectifs [famille, generation]
    M = np.zeros((len(all_ids), len(selected_gens)))
    for i in range(len(selected_gens)):
        for j, cluster in enumerate(AllClusters[i]):
            M[row[color_maps[i][j]], i] += len(cluster)

    colors = [palette_dict[fid] for fid in all_ids]

    # ── Verification : autant de couleurs que de familles ──
    n_fam = len(all_ids)
    n_col = len({tuple(round(c, 4) for c in palette_dict[fid][:3]) for fid in all_ids})
    n_transient = int((np.sum(M > 0, axis=1) == 1).sum())

    if verbose :
        print(f"[stats2] {n_fam} familles distinctes, {n_col} couleurs distinctes, "
            f"{len(colors)} bandes tracees, sur {len(selected_gens)} generations.")
        if n_col != n_fam:
            print(f"[stats2] ATTENTION : {n_fam - n_col} collision(s) de couleur "
                f"(familles de teinte quasi identique).")
        else:
            print("[stats2] OK : 1 couleur par famille.")
        if n_transient:
            print(f"[stats2] {n_transient}/{n_fam} familles n'apparaissent que sur 1 generation "
                f"-> seuil de matching peut-etre trop strict (augmente threshold).")

    ax.stackplot(selected_gens, M, colors=colors, baseline=baseline)
    ax.set_xlabel("Generation")
    ax.set_ylabel("Effectif des familles")
    #ax.set_yscale('log')

    ax.set_title("Abondance des familles au fil des generations (couleurs continues)")
    plt.tight_layout()
    return all_ids, M




# ── Graphe score ────────────────────────────────────────────────────────────────

def toHex(dec):
    digits = "0123456789ABCDEF"
    x    = dec % 16
    rest = dec // 16
    if rest == 0:
        return digits[x]
    return toHex(rest) + digits[x]

def hexToDec(hex):
    return int(hex, 16)

def hexGradient(hexStart, hexFinal, N):
    res      = []
    decStart = hexToDec(hexStart)
    decFinal = hexToDec(hexFinal)
    for i in range(N):
        n = int(decStart + ((decFinal - decStart) * i / (N - 1)))
        res.append(toHex(n))
    return res


def plotKindsCount(fig, ax, dossier, end_gen, start_gen=1, gen_step=1):
    kindsGap      = 0
    filename_base = dossier + "/Stats/Simple/gen_"

    Invalid_colors = [("#" + hex + "0000") for hex in hexGradient("79", "FF", 8)]
    Neutral_colors = [("#FF" + hex + "00") for hex in hexGradient("99", "EE", 3)]
    Valid_color    = "#007500"
    bar_width      = gen_step

    for gen in range(start_gen, end_gen + 1, gen_step):
        bottom = 0
        (_, Invalids, Neutrals, Valids, _, _, _) = ImportData.importStats(
            filename_base + str(gen))

        for i in range(len(Invalids)):
            if i > 0:
                ax.bar(gen, kindsGap, bottom=bottom, color="white",
                       width=bar_width, edgecolor="none")
                bottom += kindsGap
            size = Invalids[i]
            ax.bar(gen, size, bottom=bottom, color=Invalid_colors[i],
                   width=bar_width, edgecolor="black", linewidth=0.3)
            bottom += size

        for i in range(len(Neutrals)):
            ax.bar(gen, kindsGap, bottom=bottom, color="white",
                   width=bar_width, edgecolor="none")
            bottom += kindsGap
            size = Neutrals[i]
            ax.bar(gen, size, bottom=bottom, color=Neutral_colors[i],
                   width=bar_width, edgecolor="black", linewidth=0.3)
            bottom += size

        ax.bar(gen, kindsGap, bottom=bottom, color="white",
               width=bar_width, edgecolor="none")
        bottom += kindsGap
        ax.bar(gen, Valids, bottom=bottom, color=Valid_color,
               width=bar_width, edgecolor="black", linewidth=0)


def plotScores(fig, ax, dossier, end_gen, start_gen=1, gen_step=1):
    filename_base = dossier + "/Stats/Simple/gen_"
    X             = range(start_gen, end_gen + 1, gen_step)
    Maxs, Mins, Means = [], [], []

    for gen in range(start_gen, end_gen + 1, gen_step):
        (_, _, _, _, maxScore, minScore, meanScore) = ImportData.importStats(
            filename_base + str(gen))
        Maxs.append(maxScore)
        Mins.append(minScore)
        Means.append(meanScore)

    ax.plot(X, Maxs,  ':', color="black")
    ax.plot(X, Mins,  ':', color="black")
    ax.plot(X, Means, '-', color="blue", label='Score moyen')
    ax.set_xlabel("Génération")
    ax.set_ylabel("Score")
    ax.set_title("Évolution du score")
    plt.tight_layout()

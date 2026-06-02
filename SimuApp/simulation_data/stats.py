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

family_gap = 0.4  # espace vertical (en individus) entre deux familles dans une barre


# ── Cache ──────────────────────────────────────────────────────────────────────

_rocket_cache = {}
_cache_lock   = threading.Lock()


def _load_rocket_data(dossier, gen, idx):
    """Charge les données d'une fusée depuis le cache ou depuis le disque."""
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
    """Pré-charge en parallèle toutes les fusées d'une génération."""
    indices = [idx for cluster in all_clusters for idx in cluster]
    with ThreadPoolExecutor() as ex:
        futures = {ex.submit(_load_rocket_data, dossier, gen, idx): idx for idx in indices}
        for f in as_completed(futures):
            f.result()


def _distance_from_cache(dossier, gen1, idx1, gen2, idx2):
    """Calcule la distance en utilisant uniquement le cache (pas de lecture disque)."""
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

    return color_map_n1, next_color


# ── Histogramme ────────────────────────────────────────────────────────────────

def histo_clusters(gen_final, dossier, gen_start=1, gen_step=1):
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

    color_maps  = []
    next_color  = 0
    initial_map = {j: j for j in range(ClustersSizes[0])}
    next_color  = ClustersSizes[0]
    color_maps.append(initial_map)

    with tqdm(total=total_calls, desc="Propagation des couleurs", unit="dist") as pbar:
        for i in range(1, Nb_selected):
            color_map, next_color = propagate_colors(
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

    # Palette : répartir uniformément toutes les familles sur le spectre HSV
    all_color_ids = sorted(set(cid for cm in color_maps for cid in cm.values()))
    nb_families   = len(all_color_ids)
    id_to_idx     = {cid: i for i, cid in enumerate(all_color_ids)}

    cmap    = plt.colormaps["hsv"].resampled(nb_families + 1)
    palette = [cmap(i / nb_families) for i in range(nb_families)]

    # Affichage — clusters triés par index de palette (rouge en bas),
    # avec un espace de family_gap entre chaque famille
    fig, ax = plt.subplots()
    bar_width = gen_step * 0.9

    for i, gen in enumerate(selected_gens):
        bottom = 0
        sorted_clusters = sorted(
            range(ClustersSizes[i]),
            key=lambda j: id_to_idx[color_maps[i][j]]
        )
        for k, cluster_j in enumerate(sorted_clusters):
            if k > 0:
                # espace blanc entre deux familles
                ax.bar(gen, family_gap, bottom=bottom,
                       color="white", width=bar_width, edgecolor="none")
                bottom += family_gap

            size  = len(AllClusters[i][cluster_j])
            color = palette[id_to_idx[color_maps[i][cluster_j]]]
            ax.bar(gen, size, bottom=bottom, color=color,
                   width=bar_width, edgecolor="black", linewidth=0.3)
            bottom += size

    ax.set_xlabel("Génération")
    ax.set_ylabel("Taille des clusters")
    ax.set_title("Évolution des clusters par génération (HAC — average linkage)")
    plt.tight_layout()


# ── Main ───────────────────────────────────────────────────────────────────────

dossier   = "simu_01_06_2026_23_58_44"

gen_start = 1
gen_final = 20
gen_step  = 1

histo_clusters(gen_final, dossier, gen_start, gen_step)

plt.show()
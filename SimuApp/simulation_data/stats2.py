"""
stats2.py - Variante de stats.py : remplace l'histogramme en barres par un
streamgraph (aires empilees) de l'abondance des familles au fil des generations.

Reutilise TOUTES les optimisations de stats.py : chargement multithread
(preload_generation), cache, et propagation des couleurs entre generations
(propagate_colors) -> les couleurs sont continues au fil des generations.

Verifie aussi qu'il y a exactement autant de couleurs que de familles, et
diagnostique le nombre de familles "transitoires" (presentes sur une seule
generation), qui expliquent un nombre de couleurs surprenant quand le seuil de
matching est trop strict.

Usage :
    python stats2.py <dossier> <gen_final> [gen_start] [gen_step] [threshold]
"""

import numpy as np
import matplotlib.pyplot as plt

import stats


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
    ax.set_title("Abondance des familles au fil des generations (couleurs continues)")
    plt.tight_layout()
    return all_ids, M


def streamgraph_clusters(fig, ax, gen_final, dossier, gen_start=1, gen_step=1,
                         threshold=None, baseline="sym", use_stats_palette=False):
    """
    Charge (multithread + propagation de stats.py) puis trace le streamgraph.
    use_stats_palette=True -> palette longueur d'onde de stats.py (gradient mais
    risque de collisions) ; defaut -> palette distincte (1 couleur par famille).
    """
    if threshold is not None:
        stats.dissimilarity_threshold = threshold

    (selected_gens, _ClustersSizes, positions_dict, color_maps,
     AllClusters, palette_dict) = stats.chargeData_histo_clusters(
        dossier, gen_start, gen_final, gen_step)

    return affiche_streamgraph_clusters(
        fig, ax, selected_gens, color_maps, AllClusters, positions_dict,
        palette_dict=(palette_dict if use_stats_palette else None), baseline=baseline)


if __name__ == "__main__":
    import sys
    dossier   = sys.argv[1] if len(sys.argv) > 1 else "simu"
    gen_final = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    gen_start = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    gen_step  = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    threshold = float(sys.argv[5]) if len(sys.argv) > 5 else None

    fig, ax = plt.subplots(figsize=(12, 6))
    streamgraph_clusters(fig, ax, gen_final, dossier, gen_start, gen_step, threshold, use_stats_palette=True)
    plt.show()

"""
projection_3d.py - Visualise les familles de trajectoires (HAC) et l'impact de
l'algo genetique au fil des generations.

4 visuels (option --what) :
  embedding   : nuage 3D d'UNE generation, 1 point = 1 trajectoire, couleur = famille
  animation   : nuage 3D partage avec slider par generation -> on voit la population converger
  convergence : diversite (qui chute) vs score (qui monte) au fil des generations
  streamgraph : abondance des familles au fil des generations (Muller plot)

Fichier autonome : ne depend que de numpy + matplotlib (+ plotly optionnel pour
les rendus interactifs, + tqdm pour les barres de progression).

La distance entre trajectoires est la L1 ponderee du clustering (coefs de
RocketData.cpp::distance), calculee de facon vectorisee (bit-identique au C++).

Usage :
    python projection_3d.py <dossier> <gen>                 # embedding d'une generation
    python projection_3d.py <dossier> <gen_debut> <gen_fin> # animation+convergence+streamgraph
    python projection_3d.py <dossier> 1 200 --what animation [--step N] [--out f]
"""

import os
import sys
import argparse

import numpy as np
from tqdm import tqdm


# ── Distance du clustering : L1 ponderee (coefs de RocketData.cpp::distance) ──
# features : rMin, rMax, rMean, thetaSpan(=thetaMax-thetaMin), thetaMean, nbTurns, nbImpuls
WEIGHTS = np.array([100., 250., 25., 100., 25., 500., 10.])


# ── Lecture des fichiers (parsers inlines, pur Python) ───────────────────────

def _read_clusters(path):
    with open(path, encoding="utf-8") as f:
        lines = [x.strip() for x in f if x.strip()]
    return [[int(v) for v in lines[1 + j].split(";")] for j in range(int(lines[0]))]


def _features(path):
    """Vecteur de 7 features d'un individu (voir WEIGHTS) depuis un .rck."""
    with open(path, encoding="utf-8") as f:
        L = [x.strip() for x in f if x.strip()]
    # lignes : 3-5 theta(Min,Max,Mean), 6 nbTours, 7-9 R(Min,Max,Mean), 16 nbImpuls
    return [float(L[7]), float(L[8]), float(L[9]),
            float(L[4]) - float(L[3]), float(L[5]), float(L[6]), float(L[16])]


def _read_stats(base):
    """(score_max, score_min, score_mean) depuis <base>.stats, ou None si absent."""
    path = base + ".stats"
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as f:
        L = [x.strip() for x in f if x.strip()]
    return float(L[1]), float(L[2]), float(L[3])


def _load_gen(dossier, gen):
    """Retourne F (n,7), cluster_of_each (n,), clusters (liste d'indices)."""
    clusters = _read_clusters(os.path.join(dossier, "Stats", "HAC", f"gen_{gen}.cluster"))
    F, cl = [], []
    for cid, members in enumerate(clusters):
        for idx in members:
            F.append(_features(os.path.join(dossier, "RocketsData", f"gen_{gen}", f"ind_{idx}.rck")))
            cl.append(cid)
    return np.array(F, dtype=float), np.array(cl), clusters


# ── Distances vectorisees + MDS ───────────────────────────────────────────────

def dist_matrix(F):
    """Matrice n x n des distances (L1 ponderee), vectorisee."""
    D = np.zeros((F.shape[0], F.shape[0]))
    for k in range(F.shape[1]):
        D += WEIGHTS[k] * np.abs(F[:, k][:, None] - F[:, k][None, :])
    return D


def cross_dist(A, B):
    """Matrice |A| x |B| des distances (L1 ponderee) entre deux jeux de points."""
    D = np.zeros((A.shape[0], B.shape[0]))
    for k in range(A.shape[1]):
        D += WEIGHTS[k] * np.abs(A[:, k][:, None] - B[:, k][None, :])
    return D


def classical_mds(D, k=3):
    """MDS classique : matrice de distances n x n -> coordonnees (n, k)."""
    n = D.shape[0]
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ (D ** 2) @ H
    vals, vecs = np.linalg.eigh((B + B.T) / 2)
    top = np.argsort(vals)[::-1][:k]
    return vecs[:, top] * np.sqrt(np.clip(vals[top], 0, None))


def landmark_mds(F_all, F_land, k=3):
    """
    Landmark MDS : place exactement les pivots F_land (MDS classique), puis
    positionne tous les points F_all relativement aux pivots. O(n x L), evite la
    matrice n x n. Ici les pivots sont les barycentres de familles.
    """
    L = F_land.shape[0]
    DL2 = cross_dist(F_land, F_land) ** 2
    mu = DL2.mean(axis=0)
    H = np.eye(L) - np.ones((L, L)) / L
    B = -0.5 * H @ DL2 @ H
    vals, vecs = np.linalg.eigh((B + B.T) / 2)
    top = np.argsort(vals)[::-1][:k]
    lam = np.clip(vals[top], 1e-9, None)
    sharp = (vecs[:, top] / np.sqrt(lam))           # (L, k)
    Dx2 = cross_dist(F_all, F_land) ** 2            # (n, L)
    return -0.5 * (Dx2 - mu[None, :]) @ sharp        # (n, k)


# ── Couleurs + propagation des familles entre generations ────────────────────

def palette(n):
    import matplotlib.pyplot as plt
    if n <= 10:
        return list(plt.get_cmap("tab10").colors)[:n]
    if n <= 20:
        return list(plt.get_cmap("tab20").colors)[:n]
    cmap = plt.get_cmap("hsv")
    return [cmap(i / n) for i in range(n)]


def propagate_families(bary_per_gen, threshold):
    """
    Donne a chaque cluster un id de famille stable entre generations, par matching
    glouton des barycentres (meme logique que stats.py::propagate_colors).
    Retourne (maps, nb_familles) ou maps[i][cluster] = id_famille.
    """
    maps = [list(range(len(bary_per_gen[0])))]
    nxt = len(bary_per_gen[0])
    for i in range(1, len(bary_per_gen)):
        D = cross_dist(bary_per_gen[i], bary_per_gen[i - 1])
        pairs = sorted((D[a, b], a, b) for a in range(D.shape[0]) for b in range(D.shape[1]))
        match, used_a, used_b = {}, set(), set()
        for d, a, b in pairs:
            if d > threshold:
                break
            if a in used_a or b in used_b:
                continue
            match[a] = maps[i - 1][b]
            used_a.add(a)
            used_b.add(b)
        cm = []
        for a in range(D.shape[0]):
            if a in match:
                cm.append(match[a])
            else:
                cm.append(nxt)
                nxt += 1
        maps.append(cm)
    return maps, nxt


# ── #embedding : nuage 3D d'une generation ───────────────────────────────────

def project(dossier, gen, backend="auto", out=None, show=True, embed=True):
    print(f"Lecture generation {gen}...")
    F, cl, clusters = _load_gen(dossier, gen)
    sizes = [len(c) for c in clusters]
    print(f"  {len(F)} trajectoires, {len(clusters)} familles.")

    coords = classical_mds(dist_matrix(F))
    cols = palette(len(clusters))
    title = f"Familles de trajectoires - gen {gen} - {len(F)} traj., {len(clusters)} familles"

    if backend == "auto":
        backend = "plotly" if _has_plotly() else "matplotlib"
    if backend == "plotly":
        out = out or f"familles_g{gen}.html"
        _scatter_plotly(coords, cl, sizes, cols, title, out, embed)
    else:
        out = out or f"familles_g{gen}.png"
        _scatter_matplotlib(coords, cl, sizes, cols, title, out, show)
    print(f"Ecrit : {out}")


def _has_plotly():
    try:
        import plotly  # noqa: F401
        return True
    except ImportError:
        return False


def _scatter_matplotlib(coords, group, sizes, cols, title, out, show):
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")
    for g in range(len(sizes)):
        m = group == g
        ax.scatter(coords[m, 0], coords[m, 1], coords[m, 2], s=12, color=cols[g],
                   label=f"Famille {g} (n={sizes[g]})")
    ax.set_xlabel("MDS 1"); ax.set_ylabel("MDS 2"); ax.set_zlabel("MDS 3")
    ax.set_title(title)
    if len(sizes) <= 15:
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    if show:
        plt.show()


def _scatter_plotly(coords, group, sizes, cols, title, out, embed):
    import plotly.graph_objects as go
    fig = go.Figure()
    for g in range(len(sizes)):
        m = group == g
        r, b, gr = cols[g][:3]
        fig.add_trace(go.Scatter3d(
            x=coords[m, 0], y=coords[m, 1], z=coords[m, 2], mode="markers",
            name=f"Famille {g} (n={sizes[g]})",
            marker=dict(size=4, color=f"rgb({int(r*255)},{int(b*255)},{int(gr*255)})")))
    fig.update_layout(title=title, scene=dict(
        xaxis_title="MDS 1", yaxis_title="MDS 2", zaxis_title="MDS 3"))
    fig.write_html(out, include_plotlyjs=(True if embed else "cdn"))


# ── Chargement multi-generations (anime / convergence / streamgraph) ──────────

def _load_gens(dossier, gens):
    data = []
    for gen in tqdm(gens, desc="Lecture generations", unit="gen"):
        F, cl, clusters = _load_gen(dossier, gen)
        bary = np.array([F[cl == c].mean(axis=0) for c in range(len(clusters))])
        data.append(dict(gen=gen, F=F, cl=cl, clusters=clusters, bary=bary,
                         sizes=[len(c) for c in clusters]))
    return data


# ── #animation : nuage 3D + slider par generation (landmark MDS) ─────────────

def animation(dossier, gens, threshold=100., per_gen_cap=None,
              out="evolution_animation.html", embed=True, seed=0):
    data = _load_gens(dossier, gens)
    maps, nfam = propagate_families([d["bary"] for d in data], threshold)

    rng = np.random.default_rng(seed)
    F_all, pt_gen, pt_fam = [], [], []
    for gi, d in enumerate(data):
        fam = np.array([maps[gi][c] for c in d["cl"]])
        idx = np.arange(len(d["F"]))
        if per_gen_cap and len(idx) > per_gen_cap:
            idx = rng.choice(idx, per_gen_cap, replace=False)
        F_all.append(d["F"][idx]); pt_gen.append(np.full(len(idx), d["gen"])); pt_fam.append(fam[idx])
    F_all = np.vstack(F_all); pt_gen = np.concatenate(pt_gen); pt_fam = np.concatenate(pt_fam)

    print(f"Landmark MDS sur {len(F_all)} individus ({nfam} familles, {sum(len(d['clusters']) for d in data)} pivots)...")
    F_land = np.vstack([d["bary"] for d in data])
    coords = landmark_mds(F_all, F_land)

    import plotly.graph_objects as go
    cols = palette(nfam)
    rgb = {f: "rgb({},{},{})".format(*[int(255 * c) for c in cols[f][:3]]) for f in range(nfam)}
    ranges = [[float(coords[:, k].min()), float(coords[:, k].max())] for k in range(3)]

    def traces(gen):
        ts = []
        for f in range(nfam):
            m = (pt_gen == gen) & (pt_fam == f)
            ts.append(go.Scatter3d(
                x=coords[m, 0], y=coords[m, 1], z=coords[m, 2], mode="markers",
                name=f"Famille {f}", legendgroup=str(f),
                marker=dict(size=3, color=rgb[f])))
        return ts

    frames = [go.Frame(data=traces(g), name=str(g)) for g in gens]
    fig = go.Figure(data=traces(gens[0]), frames=frames)
    fig.update_layout(
        title="Impact de la GA : convergence de la population (slider = generation)",
        scene=dict(xaxis=dict(title="MDS 1", range=ranges[0]),
                   yaxis=dict(title="MDS 2", range=ranges[1]),
                   zaxis=dict(title="MDS 3", range=ranges[2])),
        updatemenus=[dict(type="buttons", showactive=False, x=0.1, y=0, buttons=[
            dict(label="Play", method="animate",
                 args=[None, dict(frame=dict(duration=600, redraw=True), fromcurrent=True)]),
            dict(label="Pause", method="animate",
                 args=[[None], dict(frame=dict(duration=0, redraw=False), mode="immediate")])])],
        sliders=[dict(active=0, x=0.0, y=0.0, len=1.0, pad=dict(t=50, b=10),
                      currentvalue=dict(prefix="Generation "), steps=[
            dict(label=str(g), method="animate",
                 args=[[str(g)], dict(mode="immediate", frame=dict(duration=0, redraw=True),
                                      transition=dict(duration=0))])
            for g in gens])])
    fig.write_html(out, include_plotlyjs=(True if embed else "cdn"))
    print(f"Ecrit : {out}")


# ── #convergence : diversite vs score ────────────────────────────────────────

def convergence(dossier, gens, out="evolution_convergence.png", show=True):
    div, smax, smean = [], [], []
    for gen in tqdm(gens, desc="Convergence", unit="gen"):
        F, _, _ = _load_gen(dossier, gen)
        div.append(float(cross_dist(F, F.mean(axis=0)[None, :]).mean()))
        s = _read_stats(os.path.join(dossier, "Stats", "Simple", f"gen_{gen}"))
        smax.append(s[0] if s else np.nan)
        smean.append(s[2] if s else np.nan)

    import matplotlib.pyplot as plt
    fig, ax1 = plt.subplots(figsize=(11, 6))
    ax1.plot(gens, div, "-", color="tab:blue", label="Diversite (dist. moy. au barycentre)")
    ax1.set_xlabel("Generation"); ax1.set_ylabel("Diversite", color="tab:blue")
    ax1.tick_params(axis="y", labelcolor="tab:blue")
    if not all(np.isnan(smean)):
        ax2 = ax1.twinx()
        ax2.plot(gens, smax, "--", color="tab:green", label="Score max")
        ax2.plot(gens, smean, "-", color="tab:olive", label="Score moyen")
        ax2.set_ylabel("Score", color="tab:green"); ax2.tick_params(axis="y", labelcolor="tab:green")
        lines = ax1.get_lines() + ax2.get_lines()
        ax1.legend(lines, [l.get_label() for l in lines], loc="center right")
    else:
        ax1.legend(loc="best")
    ax1.set_title("Impact de la GA : la diversite chute pendant que le score monte")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"Ecrit : {out}")
    if show:
        plt.show()


# ── #streamgraph : abondance des familles (Muller plot) ──────────────────────

def streamgraph(dossier, gens, threshold=100., out="evolution_streamgraph.png", show=True):
    data = _load_gens(dossier, gens)
    maps, nfam = propagate_families([d["bary"] for d in data], threshold)

    M = np.zeros((nfam, len(gens)))
    for gi, d in enumerate(data):
        for c, fid in enumerate(maps[gi]):
            M[fid, gi] += d["sizes"][c]

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.stackplot(gens, M, colors=palette(nfam), baseline="sym")
    ax.set_xlabel("Generation"); ax.set_ylabel("Effectif des familles")
    ax.set_title("Impact de la GA : abondance des familles au fil des generations")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"Ecrit : {out}")
    if show:
        plt.show()


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description="Familles de trajectoires (HAC) + impact de la GA.")
    p.add_argument("dossier")
    p.add_argument("gen_start", type=int)
    p.add_argument("gen_final", type=int, nargs="?", default=None,
                   help="si fourni : visuels d'evolution sur [gen_start, gen_final]")
    p.add_argument("--what", choices=["embedding", "animation", "convergence", "streamgraph", "all"],
                   default=None)
    p.add_argument("--step", type=int, default=1)
    p.add_argument("--backend", choices=["auto", "plotly", "matplotlib"], default="auto")
    p.add_argument("--out", default=None)
    p.add_argument("--threshold", type=float, default=100., help="seuil de matching des familles")
    p.add_argument("--per-gen-cap", type=int, default=None, help="animation : max d'individus par generation")
    p.add_argument("--no-show", action="store_true")
    p.add_argument("--cdn", action="store_true", help="plotly via CDN (HTML leger)")
    a = p.parse_args()
    show = not a.no_show
    what = a.what or ("embedding" if a.gen_final is None else "all")

    if what == "embedding":
        project(a.dossier, a.gen_start, a.backend, a.out, show, not a.cdn)
        return

    # animation / convergence / streamgraph : il FAUT une plage de generations
    if a.gen_final is None:
        sys.exit(f"Erreur : '--what {what}' demande une PLAGE de generations (un debut ET une fin).\n"
                 f"Exemple : python projection_3d.py {a.dossier} 1 200 --what {what}")
    gens = list(range(a.gen_start, a.gen_final + 1, a.step))
    if len(gens) < 2:
        sys.exit(f"Erreur : plage trop courte ({len(gens)} generation). Donne au moins 2 generations "
                 f"(verifie gen_debut < gen_fin et --step).")

    if what in ("animation", "all"):
        animation(a.dossier, gens, a.threshold, a.per_gen_cap,
                  out=(a.out if what == "animation" and a.out else "evolution_animation.html"),
                  embed=not a.cdn)
    if what in ("convergence", "all"):
        convergence(a.dossier, gens,
                    out=(a.out if what == "convergence" and a.out else "evolution_convergence.png"), show=show)
    if what in ("streamgraph", "all"):
        streamgraph(a.dossier, gens, a.threshold,
                    out=(a.out if what == "streamgraph" and a.out else "evolution_streamgraph.png"), show=show)


if __name__ == "__main__":
    main()

    
# use : python3.12-64.exe -u .\projection_3d.py <dossier> <end_gen> --what all
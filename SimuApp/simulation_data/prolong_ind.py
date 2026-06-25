try:
    import TIPE_SimuOrbit as cpp  # version Release (rapide)
    cpp_version = "Release"
except ImportError:
    try:
        import TIPE_SimuOrbit_d as cpp  # fallback vers Debug (pour dev)
        cpp_version = "Debug"
    except ImportError:
        raise ImportError("\x1B[38;5;202m Neither TIPE_SimuOrbit.pyd nor TIPE_SimuOrbit_d.pyd found \033[0m")

import Donnees_astres
import ImportData
import ATS
import numpy as np
import matplotlib.pyplot as plt
import cmath
import stats

#plt.style.use('dark_background')


def simulate_individual (
        dossier, gen, idx,
        fig, ax,
        verbose=False,
        is_best=False,
        trajColor="#229A28",
        plotBackground=True,
        col_start='blue', col_final='red',
        alpha=2,
        col_sun = 'yellow',
        col_startStartRing="skyblue",
        col_startFinalRing="blue",
        col_finalStartRing="lightcoral",
        col_finalFinalRing="red",
        col_rocketRing="#0A680C",
        col_final_rocket_position='black',
        linewidthStart=2, linewidthFinal=2,
        linewidthRocket=2,
        linestyleStart = '-', linestyleFinal='-',
        linestyleRocket='-'
        ) :
    filepath = dossier + "/RocketsData/gen_" + str(gen) + "/ind_" + str(idx) + ".rck"

    if is_best :
        filepath = dossier + "/Bests/RocketsData/gen_" + str(gen) + ".rck"


    (_, _, dt,
    _,
    planete_depart, planete_arrivee,
    _,
    pos_init_depart, pos_init_arrivee,
    pos_final_depart, pos_final_arrivee,
    _,
    _,
    est_etat_lie, r_min, r_max) = ImportData.lire_donnees(dossier + "/Bests/Physics/gen_" + str(gen) + ".phys")

    if not is_best :
        est_etat_lie, r_min, r_max = False, -1, -1

    x_depart, y_depart, x_final, y_final = [], [], [], []
    ImportData.init_liste(dossier + "/start_planet.txt", x_depart, y_depart)
    ImportData.init_liste(dossier + "/final_planet.txt", x_final, y_final)

    (
        _, 
        tof, max_time,
        _, _, _,
        _,
        _, _, _,
        startPlanet, finalPlanet,
        _, _,
        _,
        imp0, impulsions_, dates,
        nbImpulsions,
        genome
   ) = ImportData.importRocketData(filepath)
    
    impulsions = [(1e-15, imp0)]
    for i in range(nbImpulsions-1):
        impulsions.append((dates[i], impulsions_[i]))

    X_rocket, Y_rocket = cpp.genetic.GetIndividualTrajectory__configD32_2_2(
        dt_seconds = dt,
        max_time_days = max_time,
        simulation_duration_days= alpha*tof,
        tof_days=tof,
        mass_kg = 1e3,
        startPlanet = startPlanet,
        finalPlanet = finalPlanet,
        genome = genome
    )

    pos_init_depart, pos_final_depart, pos_init_arrivee, pos_final_arrivee, x_depart, y_depart, x_final, y_final, X_rocket, Y_rocket = ATS.conversion_data(pos_init_depart, pos_final_depart, pos_init_arrivee, pos_final_arrivee, x_depart, y_depart, x_final, y_final, X_rocket, Y_rocket)

    fig, ax = ATS.affiche_orbite(
        fig, ax,
        pos_init_depart, pos_final_depart, pos_init_arrivee, pos_final_arrivee,
        x_depart, y_depart, x_final, y_final,
        X_rocket, Y_rocket,
        planete_depart, planete_arrivee, 
        est_etat_lie, r_min, r_max, 
        impulsions, dt, trajColor=trajColor,
        plotBackground=plotBackground,
        col_start=col_start, col_final=col_final,
        col_sun = col_sun,
        col_startStartRing=col_startStartRing,
        col_startFinalRing=col_startFinalRing,
        col_finalStartRing=col_finalStartRing,
        col_finalFinalRing=col_finalFinalRing,
        col_rocketRing=col_rocketRing,
        col_final_rocket_position=col_final_rocket_position,
        linewidthStart=linewidthStart, linewidthFinal=linewidthFinal,
        linewidthRocket=linewidthRocket,
        linestyleStart = linestyleStart, linestyleFinal=linestyleFinal,
        linestyleRocket=linestyleRocket
    )


    ax.set_aspect('equal')
    ax.set_title(f"Trajectoire des objets du système - gen {gen}")
    ax.set_xlabel("x (km)")
    ax.set_ylabel("y (km)")
    ax.grid(False)

    return fig, ax




def plot_generation_families(fig, ax, dossier, gen_number, selected_gens,
                              AllClusters, color_maps, palette_dict, positions_dict):
    """
    Affiche un histogramme des familles présentes à une génération donnée.
    Abscisse : indice local de la famille (trié par position MDS).
    Ordonnée : taille du cluster.
    Couleur   : couleur de la famille dans la palette globale.

    Retourne la liste ordonnée des color_ids affichés, pour pouvoir les
    passer directement à plot_families_trajectories.
    """
    gen_idx = selected_gens.index(gen_number)
    clusters  = AllClusters[gen_idx]
    color_map = color_maps[gen_idx]

    # Trier les clusters par position MDS pour que les couleurs proches soient adjacentes
    sorted_clusters = sorted(
        range(len(clusters)),
        key=lambda j: positions_dict.get(color_map[j], 0.5)
    )

    color_ids_ordered = [color_map[j] for j in sorted_clusters]
    sizes             = [np.log(len(clusters[j]) + 1) for j in sorted_clusters]
    colors            = [palette_dict[cid] for cid in color_ids_ordered]
    x                 = np.arange(len(sorted_clusters))

    ax.bar(x, sizes, color=colors, width=0.6, edgecolor="black", linewidth=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(
        [str(cid) for cid in color_ids_ordered],
        fontsize=7, rotation=45, ha="right"
    )
    ax.set_xlabel("Identifiant de famille (color_id)")
    ax.set_ylabel("Taille du cluster (échelle log)")
    ax.set_title(f"Familles présentes à la génération {gen_number} "
                 f"({len(sorted_clusters)} familles)")

    plt.tight_layout()
    return color_ids_ordered


def plot_families_trajectories(fig, ax, dossier, selected_gens, AllClusters, color_maps,
                                palette_dict, family_color_ids,
                                gen_indices=None, activateLegend=False, limits=3):
    """
    Affiche sur un unique ax les trajectoires représentatives des familles
    sélectionnées, sur potentiellement plusieurs générations.

    Paramètres
    ----------
    fig, ax           : figure et axes matplotlib sur lesquels tracer
    dossier           : répertoire de la simulation
    selected_gens     : retour de chargeData_histo_clusters
    AllClusters       : retour de chargeData_histo_clusters
    color_maps        : retour de chargeData_histo_clusters
    palette_dict      : retour de chargeData_histo_clusters
    family_color_ids  : liste de color_ids à afficher (ex: [3, 7, 12])
    gen_indices       : indices (0-based) dans selected_gens à inspecter.
                        Si None, toutes les générations sont inspectées.
    dt_seconds        : pas de temps pour la simulation (défaut 3600s)
    """
    #family_color_ids = set(family_color_ids)

    if gen_indices is None:
        gen_indices = range(len(selected_gens))

    k=0


    for i in range(len(family_color_ids)) :
        gen_idx = gen_indices[i]

        gen_number = selected_gens[gen_idx]
        clusters   = AllClusters[gen_idx]
        color_map  = color_maps[gen_idx]

        representative = [
            rep for rep in stats.get_generation_representatives(
                dossier, gen_idx, gen_number, clusters, color_map, palette_dict
            )
            if rep["color_id"] == family_color_ids[i]
        ][0]

        fid   = representative["color_id"]
        color = palette_dict[fid]
        label = f"Famille {fid} — gen {gen_number} ({representative['cluster_size']} ind.)"

        n_lines_before = len(ax.get_lines())
        simulate_individual(
            dossier, gen_number, representative["idx"],
            fig, ax,
            trajColor=color,
            plotBackground=True if i==0 else False,
            col_start="#000000", col_final="#000000")
        k+=1

        line = ax.get_lines()[-1]
        line.set_label(label)

    if activateLegend :
        ax.legend(fontsize=7)
    ax.set_aspect("equal")
    ax.set_title(
        f"Trajectoires des familles {sorted(family_color_ids)}",
        fontsize=10
    )

    troisCentMillions = limits * 1.5e8
    ax.set_xlim(-troisCentMillions, troisCentMillions)
    ax.set_ylim(-troisCentMillions, troisCentMillions)
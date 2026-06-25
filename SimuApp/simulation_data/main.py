from ATS import affiche_fichier, affiche_couple_generations
import matplotlib.pyplot as plt
import os
import prolong_ind

plt.close("all")

dossier = "simu_07_06_2026_15_43_01"
generation = 31
verbose = True

# fig, ax = affiche_fichier(dossier, generation, window_title="Terre_jupiter", verbose=verbose)


fig, ax = plt.subplots()

prolong_ind.simulate_individual(
    dossier=dossier, gen=generation, idx=0,
    fig=fig, ax=ax,
    verbose=verbose,
    is_best=True,
    trajColor="#0000FF",
    plotBackground=True,
    col_start="#0000FF00", col_final='#000000',
    alpha=500/478.75,
    # alpha=275/236.333333,
    # alpha=1,
    col_sun = "#E2E200",
    col_startStartRing="#076807",
    col_startFinalRing="#0000FF00",
    col_finalStartRing="#FF005D00",
    col_finalFinalRing="#FF0000",
    col_rocketRing="#02879B",
    col_final_rocket_position="#00000000",
    linewidthStart=2, linewidthFinal=2,
    linewidthRocket=2,
    linestyleStart = '-', linestyleFinal='--',
    linestyleRocket='-'
)

# prolong_ind.simulate_individual(
#     dossier=dossier, gen=generation, idx=0,
#     fig=fig, ax=ax,
#     verbose=verbose,
#     is_best=True,
#     trajColor="#0CDA51",
#     plotBackground=True,
#     col_start='#2383C2', col_final='#FF0038',
#     alpha=1.2,
#     col_sun = '#FFDE21',
#     col_startStartRing="#2383C2",
#     col_startFinalRing="#2383C2",
#     col_finalStartRing="#FF0038",
#     col_finalFinalRing="#FF0038",
#     col_rocketRing="#0CDA5100",
#     col_final_rocket_position='black',
#     linewidthStart=4, linewidthFinal=4,
#     linewidthRocket=4,
#     linestyleStart = '--', linestyleFinal='--',
#     linestyleRocket='-'
# )

plt.show()
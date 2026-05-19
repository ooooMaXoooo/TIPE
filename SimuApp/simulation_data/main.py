from ATS import affiche_fichier
import matplotlib.pyplot as plt

is_super = False

dossier = "simu_18_05_2026_15_57_44__cool"
generation = 41

if is_super :
    dossier += "__super"

fig, ax = affiche_fichier(dossier, generation)

plt.close("all")
plt.show()
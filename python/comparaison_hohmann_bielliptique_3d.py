import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ----------------------------
# Constantes
# ----------------------------
G = 6.67e-11
M_soleil = 1.989e30
sqrt_u = np.sqrt(G * M_soleil)
UA = 150e9

# ----------------------------
# Fonctions Δv
# ----------------------------
def delta_v(r1, r2):
    dv1 = sqrt_u * np.abs(np.sqrt(2/r1 - 2/(r1+r2)) - np.sqrt(1/r1))
    dv2 = sqrt_u * np.abs(np.sqrt(1/r2) - np.sqrt(2/r2 - 2/(r2+r1)))
    return dv1 + dv2

def delta_u(r1, r2, r3):
    du1 = sqrt_u * np.abs(np.sqrt(2/r1 - 2/(r1+r3)) - np.sqrt(1/r1))
    du2 = sqrt_u * np.abs(np.sqrt(2/r3 - 2/(r3+r2)) - np.sqrt(2/r3 - 2/(r3+r1)))
    du3 = sqrt_u * np.abs(np.sqrt(1/r2) - np.sqrt(2/r2 - 2/(r2+r3)))
    return du1 + du2 + du3

# ----------------------------
# Plages de r1, r2, r3
# ----------------------------
nb_points = 15
R1 = np.linspace(0.72*UA, 30.11*UA, nb_points)
R2 = np.linspace(1*UA, 30.11*UA, nb_points)
R3 = np.linspace(1*UA, 50*UA, nb_points)

# ----------------------------
# Calcul des données
# ----------------------------
r1_list, r2_list, r3_list, diff_list = [], [], [], []

for r1 in R1:
    for r2 in R2:
        if r2 <= r1:  # condition r2 > r1
            continue
        if np.abs(r2 - r1) < 0.28*UA:  # condition distance minimale
            continue
        for r3 in R3:
            if r3 <= r2:  # condition r3 > r2
                continue
            dv = delta_v(r1, r2)
            du = delta_u(r1, r2, r3)
            diff = dv - du
            r1_list.append(r1/UA)
            r2_list.append(r2/UA)
            r3_list.append(r3/UA)
            diff_list.append(diff)

# Conversion en arrays
r1_arr = np.array(r1_list)
r2_arr = np.array(r2_list)
r3_arr = np.array(r3_list)
diff_arr = np.array(diff_list)

# ----------------------------
# Filtrer les points Bi-elliptique moins coûteux
# ----------------------------
mask_bi_less_costly = diff_arr > 0
r1_bi = r1_arr[mask_bi_less_costly]
r2_bi = r2_arr[mask_bi_less_costly]
r3_bi = r3_arr[mask_bi_less_costly]
diff_bi = diff_arr[mask_bi_less_costly]

# ----------------------------
# Analyse des rapports critiques
# ----------------------------
ratio_r2_r1 = r2_bi / r1_bi
ratio_r3_r2 = r3_bi / r2_bi

print("=== Analyse des conditions empiriques avec conditions supplémentaires ===")
print(f"Nombre de points Bi-elliptique moins coûteux : {len(r1_bi)}")
print(f"r1 [UA] min/max : {r1_bi.min():.2f} / {r1_bi.max():.2f}")
print(f"r2 [UA] min/max : {r2_bi.min():.2f} / {r2_bi.max():.2f}")
print(f"r3 [UA] min/max : {r3_bi.min():.2f} / {r3_bi.max():.2f}")
print(f"Ratio r2/r1 min/max : {ratio_r2_r1.min():.2f} / {ratio_r2_r1.max():.2f}")
print(f"Ratio r3/r2 min/max : {ratio_r3_r2.min():.2f} / {ratio_r3_r2.max():.2f}")

# ----------------------------
# Scatter 3D pour visualiser
# ----------------------------
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111, projection='3d')
colors = ['blue' if d < 0 else 'red' for d in diff_arr]
ax.scatter(r1_arr, r2_arr, r3_arr, c=colors, marker='o')
ax.set_xlabel("r1 [UA]")
ax.set_ylabel("r2 [UA]")
ax.set_zlabel("r3 [UA]")
ax.set_title("Hohmann vs Bi-elliptique : bleu = Bi-elliptique moins coûteux, rouge = Hohmann moins coûteux")
plt.show()

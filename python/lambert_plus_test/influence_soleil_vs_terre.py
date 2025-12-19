from scipy.constants import G
import numpy as np
import matplotlib.pyplot as plt

# Constantes
G_km3s2 = G / 1e9  # [km^3/s^2]
M_sun = 1.9885e30  # [kg]
M_earth = 5.972e24  # [kg]
D = 1.496e8  # Distance Terre-Soleil [km]

# Rayon de la sphère d'influence (SOI) de la Terre
R_SOI_earth = D * (M_earth / M_sun)**(2/5)

# Distances au centre de la Terre (logarithmique)
r = np.logspace(np.log10(0.01 * R_SOI_earth), np.log10(10 * R_SOI_earth), 500)

# Accélérations gravitationnelles
a_earth = G_km3s2 * M_earth / r**2
a_sun = G_km3s2 * M_sun / (D - r)**2

# Influence en pourcentage
influence_earth_pct = 100 * a_earth / (a_earth + a_sun)
influence_sun_pct = 100 * a_sun / (a_earth + a_sun)

# Trouver le point où les deux influences sont égales (~50%-50%)
idx_equilibrium = np.argmin(np.abs(influence_sun_pct - 50))
r_eq = r[idx_equilibrium]

# --- Affichage graphique ---
plt.figure(figsize=(10,6))
plt.plot(r / R_SOI_earth, influence_sun_pct, label="Influence Soleil [%]", color="orange")
plt.plot(r / R_SOI_earth, influence_earth_pct, label="Influence Terre [%]", color="blue")

# Lignes verticales repères
plt.axvline(x=1, color='gray', linestyle='--', label="Sphère d'influence théorique")
plt.axvline(x=r_eq / R_SOI_earth, color='red', linestyle='--', label="Zone d'équilibre (~50/50)")

# Axes et style
plt.xscale('log')
plt.xlabel("Distance (en multiples de R_SOI)")
plt.ylabel("Pourcentage d'influence gravitationnelle")
plt.title("Part relative de l'influence gravitationnelle Terre / Soleil")
plt.legend()
plt.grid(True, which='both')
plt.tight_layout()
plt.show()

# --- Affichage numérique ---
print(f"Rayon SOI de la Terre ≈ {R_SOI_earth:,.0f} km")
print(f"Point d'équilibre Terre/Soleil ≈ {r_eq:,.0f} km ({r_eq / R_SOI_earth:.2f} × R_SOI)")

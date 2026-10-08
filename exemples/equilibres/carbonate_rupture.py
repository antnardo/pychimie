"""Décomposition du carbonate de calcium : rupture d'équilibre à pression fixée.

CaCO₃(s) = CaO(s) + CO₂(g), ΔrH° et ΔrS° par la loi de Hess. Sous P = 1 bar, CO₂ est le
seul gaz : son activité vaut P/P° quelle que soit sa quantité, et Q ne dépend pas de ξ.
Si K° < P/P°, rien ne se décompose ; si K° > P/P°, tout se décompose — G(ξ) est
monotone, minimal sur une borne. Dans un récipient de volume fixé, a(CO₂) = nRT/(VP°)
varie avec ξ et l'équilibre redevient possible, tant qu'il reste du carbonate.
"""

import matplotlib.pyplot as plt
import numpy as np

from pychimie import Reaction, Systeme, balayage

calcaire = Reaction.depuis_tables("CaCO3(s) = CaO(s) + CO2(g)")
print(calcaire, f"; K° = 1 à T = {calcaire.DrH0() / calcaire.DrS0():.0f} K")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
for T in (1000, 1110, 1200):
    s = Systeme(calcaire, {"CaCO3": 1}, T=T, P=1)
    etat = s.etat_final()
    xi = np.linspace(0, 1, 200)
    ax1.plot(xi, s.G(xi) / 1e3, label=f"{T} K : $\\xi_f$ = {etat.xi:g}, K° = {etat.K0:.2f}")
ax1.set(xlabel=r"$\xi$ (mol)", ylabel=r"$G(\xi) - G(0)$ (kJ)", title="$P$ = 1 bar")
ax1.legend()

T = np.linspace(950, 1250, 121)
for n in (1.0, 0.1):
    s = Systeme(calcaire, {"CaCO3": n}, V_gaz=10, T=1000)
    b = balayage(lambda T, s=s: s.avec(T=T), T)
    ax2.plot(T, b.taux_avancement, label=f"{n:g} mol de CaCO₃, V = 10 L")
    rupture = T[~b.equilibre]
    if rupture.size:
        print(f"{n:g} mol, V = 10 L : rupture d'équilibre dès {rupture[0]:.0f} K")
b = balayage(lambda T: Systeme(calcaire, {"CaCO3": 1}, T=T, P=1), T)
ax2.plot(T, b.taux_avancement, "k--", label="$P$ = 1 bar")
ax2.set(xlabel="$T$ (K)", ylabel=r"$\xi_f/\xi_\mathrm{max}$")
ax2.legend()
for ax in (ax1, ax2):
    ax.grid(True)
fig.tight_layout()
plt.show()

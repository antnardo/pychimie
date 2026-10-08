"""Synthèse de l'ammoniac : K°(T) et déplacements d'équilibre, avec les tables.

ΔrH°, ΔrS° et ΔrC°p viennent des tables par la loi de Hess. On compare K°(T) en
approximation d'Ellingham et avec Kirchhoff, puis on perturbe l'équilibre d'un mélange
stœchiométrique (700 K, 200 bar) : hausse de T, de P, ajout d'argon à P puis à V fixés,
ajout de N₂. Les courbes donnent le taux d'avancement en fonction de T, de P et de la
quantité d'argon.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import R

from pychimie import Reaction, Systeme, balayage, perturber

kirchhoff = Reaction.depuis_tables("N2(g) + 3 H2(g) = 2 NH3(g)", kirchhoff=True)
ellingham = kirchhoff.ellingham()
print(kirchhoff)

s = Systeme(kirchhoff, {"N2": 1, "H2": 3}, T=700, P=200)
etat = s.etat_final()
print(f"équilibre : ξ = {etat.xi:.4f} mol, x(NH3) = {etat.fraction('NH3'):.4f}")
# Le même équilibre dans un récipient fermé : le volume du mélange à 200 bar, en litres.
V = sum(etat.quantites) * R * 700 / 200e5 * 1e3
a_volume_fixe = Systeme(kirchhoff, {"N2": 1, "H2": 3}, T=700, V_gaz=V)
perturbations = {
    "T : 700 → 800 K": perturber(s, T=800),
    "P : 200 → 300 bar": perturber(s, P=300),
    "+1 mol Ar à P fixée": perturber(s, ajout={"Ar(g)": 1}),
    "+1 mol Ar à V fixé": perturber(a_volume_fixe, ajout={"Ar(g)": 1}),
    "+0,1 mol N2": perturber(s, ajout={"N2": 0.1}),
}
for nom, p in perturbations.items():
    print(f"{nom:22} Q/K° = {p.Q / p.K0:8.4f}  sens {p.sens:+d}  δξ = {p.deplacement:+.5f} mol")

fig, axes = plt.subplots(1, 3, figsize=(13, 4))
T = np.linspace(300, 1200, 200)
axes[0].semilogy(T, kirchhoff.K0(T), label="Kirchhoff ($\\Delta_rC_p^\\circ$ constant)")
axes[0].semilogy(T, ellingham.K0(T), "--", label="Ellingham")
axes[0].set(xlabel="$T$ (K)", ylabel="$K^\\circ$")
axes[0].legend()

for P in (1, 50, 200):
    b = balayage(lambda T, P=P: s.avec(T=T, P=P), np.linspace(400, 1000, 60))
    axes[1].plot(b.valeurs, b.taux_avancement, label=f"$P = {P}$ bar")
axes[1].set(xlabel="$T$ (K)", ylabel=r"$\xi_\mathrm{éq}/\xi_\mathrm{max}$")
axes[1].legend()

n_Ar = np.linspace(0, 20, 41)
b = balayage(lambda n: s.avec(inertes={"Ar(g)": n}), n_Ar)
axes[2].plot(n_Ar, b.taux_avancement, label="$P$ fixée")
b = balayage(lambda n: a_volume_fixe.avec(inertes={"Ar(g)": n}), n_Ar)
axes[2].plot(n_Ar, b.taux_avancement, "--", label="$V$ fixé")
axes[2].set(xlabel="$n(\\mathrm{Ar})$ (mol)", ylabel=r"$\xi_\mathrm{éq}/\xi_\mathrm{max}$")
axes[2].legend()
for ax in axes:
    ax.grid(True)
fig.tight_layout()
plt.show()

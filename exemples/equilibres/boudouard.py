"""Équilibre de Boudouard C(s) + CO₂(g) = 2 CO(g) : Q(ξ), ξ_éq, et rupture d'équilibre.

Reprend ``Optimisation_Boudouard.py`` (2 mol de C, 1 mol de CO₂ et de CO, 2 bar,
1200 K), puis refait le calcul avec seulement 0,1 mol de carbone : Q(ξ_max) reste sous
K°, le carbone disparaît avant l'équilibre et G(ξ) est minimal sur la borne. ΔrH° et
ΔrS° sont ceux des tables NBS, que la loi de Hess redonne.
"""

import matplotlib.pyplot as plt

from pychimie import Reaction, Systeme

boudouard = Reaction("C(s) + CO2(g) = 2 CO(g)", DrH0=172.5e3, DrS0=175.8)
print("Hess :", Reaction.depuis_tables("C(s) + CO2(g) = 2 CO(g)"))
T, P = 1200, 2  # K, bar

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
s = Systeme(boudouard, {"C": 2, "CO2": 1, "CO": 1}, T=T, P=P)
xi, xi_eq = s.avancements(2000), s.xi_eq()
print(f"ξ ∈ [{s.xi_min:g}, {s.xi_max:g}] mol, K° = {s.K0():.4g}, ξ_éq = {xi_eq:.6f} mol")
ax1.semilogy(xi, s.Q(xi), label="$Q(\\xi)$")
ax1.axhline(s.K0(), color="C1", label=f"$K^\\circ$ à $T = {T - 273}$ °C")
ax1.axvline(xi_eq, color="C2", ls="--", label=f"$\\xi_\\mathrm{{éq}} = {xi_eq:.2f}$ mol")
ax1.set(xlabel=r"Avancement $\xi$ (mol)", ylabel=f"Quotient réactionnel pour $P = {P}$ bar")
ax1.legend()
ax1.grid(True)

peu = Systeme(boudouard, {"C": 0.1, "CO2": 1}, T=T, P=P)
etat = peu.etat_final()
print(
    f"0,1 mol de C : ξ_f = {etat.xi:g} mol, équilibre : {etat.equilibre}, "
    f"disparus : {', '.join(etat.disparus)}, Q = {etat.Q:.3g} < K° = {etat.K0:.3g}"
)
xi = peu.avancements(500, debut=0, fin=peu.xi_max)
ax2.plot(xi, peu.G(xi) / 1e3, label=r"$G(\xi) - G(0)$, 0,1 mol de C")
ax2.axvline(peu.xi_max, color="C3", ls="--", label=r"$\xi_\mathrm{max}$ : C(s) disparaît")
ax2.set(xlabel=r"Avancement $\xi$ (mol)", ylabel="kJ")
ax2.legend()
ax2.grid(True)
fig.tight_layout()
plt.show()

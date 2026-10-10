"""Synthèse de l'ammoniac : Q(ξ), K°(T) et ξ_éq sous 100 bar à 400 K.

Reprend ``Optimisation_NH3_1_P.py`` : 1 mol de chaque gaz, ΔrH° et ΔrS° du cours. Q(ξ)
croît de 0 (plus de NH₃) à l'infini (plus de H₂) et coupe K° à l'équilibre ; la seconde
figure montre ΔrG(ξ) = RT ln(Q/K°), qui s'annule au même ξ, et G(ξ), minimal là.
"""

import matplotlib.pyplot as plt

from pychimie import Reaction, Systeme

nh3 = Reaction("N2(g) + 3 H2(g) = 2 NH3(g)", DrH0=-92e3, DrS0=-175.6)
T, P = 400, 100  # K, bar
s = Systeme(nh3, {"N2": 1, "H2": 1, "NH3": 1}, T=T, P=P)
xi = s.avancements(2000)
xi_eq = s.xi_eq()
print(f"ξ ∈ [{s.xi_min:g}, {s.xi_max:.4f}] mol, K° = {s.K0():.4g}, ξ_éq = {xi_eq:.6f} mol")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
ax1.semilogy(xi, s.Q(xi), label="$Q(\\xi)$")
ax1.axhline(s.K0(), color="C1", label=f"$K^\\circ$ à $T = {T - 273}$ °C")
ax1.axvline(xi_eq, color="C2", ls="--", label=f"$\\xi_\\mathrm{{éq}} = {xi_eq:.2f}$ mol")
ax1.set(xlabel=r"Avancement $\xi$ (mol)", ylabel=f"Quotient réactionnel pour $P = {P}$ bar")
ax1.legend()
ax1.grid(True)

ax2.plot(xi, s.DrG(xi) / 1e3, label=r"$\Delta_rG(\xi)$ (kJ/mol)")
ax2.plot(xi, s.G(xi) / 1e3, label=r"$G(\xi) - G(0)$ (kJ)")
ax2.axhline(0, color="k", lw=0.5)
ax2.axvline(xi_eq, color="C2", ls="--")
ax2.set(xlabel=r"Avancement $\xi$ (mol)", ylim=(-60, 60))
ax2.legend()
ax2.grid(True)
fig.tight_layout()
plt.show()

"""Synthèse de l'ammoniac : influence des proportions initiales, sous 50 bar à 450 K.

Reprend ``Optimisation_NH3_3_r.py`` : 1 mol de N₂ + H₂, x = n(N₂)/(n(N₂) + n(H₂)). La
fraction molaire de NH₃ à l'équilibre est maximale pour les proportions
stœchiométriques, x = 1/4, alors que le taux d'avancement r = ξ_éq/ξ_max croît encore.
"""

import matplotlib.pyplot as plt
import numpy as np

from pychimie import Reaction, Systeme, balayage

nh3 = Reaction("N2(g) + 3 H2(g) = 2 NH3(g)", DrH0=-92e3, DrS0=-175.6)
s = Systeme(nh3, {"N2": 1, "H2": 1}, T=450, P=50)

x_N2 = np.linspace(0.01, 1, 400, endpoint=False)
b = balayage(lambda x: s.avec(initial={"N2": x, "H2": 1 - x}), x_N2)
print(f"x(NH3) maximal pour x(N2) = {x_N2[np.argmax(b.fraction('NH3'))]:.3f}")

fig, ax = plt.subplots()
ax.axvline(1 / 4, color="k", ls="--")
ax.plot(x_N2, b.fraction("NH3"), label=r"$x_\mathrm{éq}(\mathrm{NH_3})$")
ax.plot(x_N2, b.xi_max, label=r"$\xi_\mathrm{max}$")
ax.plot(x_N2, b.xi, label=r"$\xi_\mathrm{éq}$")
ax.plot(x_N2, b.taux_avancement, label=r"$r = \xi_\mathrm{éq}/\xi_\mathrm{max}$")
ax.axis([0, 1, 0, 1.05])
ax.set(xlabel=r"Fraction molaire initiale de $\mathrm{N_2}$")
ax.legend()
plt.show()

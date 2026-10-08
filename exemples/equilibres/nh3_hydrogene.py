"""Synthèse de l'ammoniac : Q(ξ) depuis l'état initial pour dix quantités de H₂.

Reprend ``Optimisation_NH3_2_ni.py`` : 1 mol de N₂ et de NH₃, n(H₂) de 10⁻³ à 10 mol,
sous P = P°. Moins il y a de H₂, plus ξ_max est petit et plus Q y diverge tôt.
"""

import matplotlib.pyplot as plt
import numpy as np

from pychimie import Reaction, Systeme

nh3 = Reaction("N2(g) + 3 H2(g) = 2 NH3(g)", DrH0=-92e3, DrS0=-175.6)

fig, ax = plt.subplots()
for n_H2 in np.logspace(-3, 1, 10):
    s = Systeme(nh3, {"N2": 1, "H2": n_H2, "NH3": 1})
    xi = s.avancements(debut=0)
    ax.loglog(xi, s.Q(xi), label=f"$n(\\mathrm{{H_2}}) = {n_H2:.1e}$ mol")
ax.axis([1e-5, 10, 1e-2, 1e14])
ax.set(xlabel=r"Avancement $\xi$ (mol)", ylabel=r"Quotient réactionnel pour $P = P^\circ$")
ax.legend(fontsize="small")
ax.grid(True)
plt.show()

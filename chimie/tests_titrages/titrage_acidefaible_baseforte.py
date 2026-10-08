#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# print(S.melange['Acide test'].Ca(pH = 12))
import numpy as np
import matplotlib.pyplot as plt

from chimie.titrage import PolyAcide, p as plog

T = PolyAcide([6.3, 10.3])  # carbonates

# Titrages par une base forte
CB = 0.1
plt.figure(figsize=(6, 4))
for CT in [0.01, 0.1, 1]:

    pH_init = T.pH(Ctot=CT)
    pHs = np.linspace(pH_init, 14-plog(CB)-.2, 200)

    n_article = T.n(Ctot=CT, pH=pHs)
    n = T.n_dilue(CA=CT, CB=CB, pH=pHs)
    n2 = min(np.floor(2 * max(n)) / 2, 3)

    lines = plt.plot(n, pHs, label=f"Ctot={CT:.1e} M")
    plt.plot(n_article, pHs, '--', lw=.7, color=lines[0].get_color())

    EPs = []
    for p in np.arange(0, n2 + 0.5, step=0.5):
        EP = T.pH_EP(CA=CT, CB=CB, p=p)
        EPs.append((p, EP))
        if p < 3:
            plt.plot([0, p], [EP, EP], "o-", lw=.5, color=lines[0].get_color())
            plt.text(-0.05, EP, f"{EP:.1f}", fontsize=10, ha='right', va='center')

    print(EPs)

plt.xlim(0, 3)
plt.legend()
plt.ylabel("$p$H")
plt.xlabel(r"$n = V \dfrac{C_\text{burette}}{C_\text{acide faible}V_\text{becher}}$")
plt.grid()
plt.tight_layout()
plt.savefig("acidefaible_baseforte.pdf")
plt.show()

#!python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

from chimie.titrage import PolyAcide

T = PolyAcide([9.3])

# Titrages par une base forte
plt.figure(figsize=(6, 4))
for CT in [0.1]:

    pHs = np.linspace(.5, 13.5, 200)

    n = T.n(Ctot=CT, pH=pHs)
    n1 = np.ceil(2 * min(n)) / 2
    n2 = np.floor(2 * max(n)) / 2

    lines = plt.plot(1-n, pHs, label=f"Ctot={CT:.1e} M")

    EPs = []
    for p in np.arange(n1, n2 + 0.5, step=0.5):
        EP = T.pH_EP(CA=CT, CB=0, p=p)
        EPs.append((p, EP))
        if -1 <= p <= 1:
            plt.plot([0, 1-p], [EP, EP], "o--", lw=.5, color=lines[0].get_color())
            plt.text(-0.05, EP, f"{EP:.1f}", fontsize=10, ha='right', va='center')

    print(EPs)

plt.xlim(0, 2)
plt.legend()
plt.ylabel("$p$H")
plt.xlabel("$1-n$")
plt.grid()
plt.tight_layout()
plt.savefig("basefaible_acidefort.pdf")
plt.show()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# print(S.melange['Acide test'].Ca(pH = 12))
import numpy as np
import matplotlib.pyplot as plt
from timeit import timeit

from chimie.titrage import PolyAcide

T = PolyAcide([6.3, 10.3])  # carbonates

pH = np.linspace(0, 14, 1000)

plt.figure(figsize=(5, 5))
for CT in [0.00001, 0.0001, 0.001, 0.01, 0.1, 1]:
    n = T.n(Ctot=CT, pH=pH)
    plt.plot(n, pH, label=f"CT={CT:.1e}M")

plt.xlim(-2, 4)
plt.legend()
plt.ylabel("$p$H")
plt.xlabel(r"$n=\frac{C_B}{C_A}$ ($n<0$ : acide fort)")
plt.grid()
plt.tight_layout()
plt.savefig("acidefaible_melangebaseforte.pdf")
plt.show()

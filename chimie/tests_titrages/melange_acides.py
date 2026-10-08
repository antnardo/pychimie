#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

from chimie.titrage import PolyAcide, Melange

H2CO3 = PolyAcide([6.3, 10.3])  # carbonates
H3PO4 = PolyAcide([2.1, 7.2, 12.7])  # phosphates

M = Melange()
M.add_polyacide(H2CO3, 0.1)
M.add_polyacide(H3PO4, 0.1)

# Mélange avec une base forte
plt.figure(figsize=(6, 4))

pHs = np.linspace(.5, 13.5, 200)
n = M.n(pH=pHs)
plt.plot(n, pHs)

plt.xlim(-3, 3)
plt.ylabel("$p$H")
plt.xlabel("$n$")
plt.grid()
plt.tight_layout()
plt.show()

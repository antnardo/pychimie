#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

from chimie.titrage import PolyAcide, Melange, p


HCL = PolyAcide([-3.5])
NAOH = PolyAcide([15.7])

CA = 0.1
CB = 0.1
M = Melange()
M.add_polyacide(HCL, CA)
M.add_polyacide(NAOH, CB)

y = CB/CA

def pH_acide(n):
    return p(CB) + p((1-n)/(y+n))

def pH_basique(n):
    return 14 - p(CB) - p((n-1)/(y+n))

n1 = np.linspace(0, .999, 1000)
n2 = np.linspace(1.001, 2, 1000)

plt.figure(figsize=(6, 4))

pHs = np.linspace(.5, 13.5, 200)
n = M.n(pH=pHs)
plt.plot(2*n, pHs)

plt.plot(n1, pH_acide(n1))
plt.plot(n2, pH_basique(n2))


plt.xlim(0, 2)
plt.ylim(.5, 13.5)
plt.ylabel("$p$H")
plt.xlabel("$n$")
plt.grid()
plt.tight_layout()
plt.show()

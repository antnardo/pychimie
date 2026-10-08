#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

from chimie.titrage import p, AcideFort

CA = 0.1
CB = 0.4
y = CB/CA

HCL = AcideFort()

# Calcul direct
def pH_acide(n):
    return p(CB) + p((1-n)/(y+n))

def pH_basique(n):
    return 14 - p(CB) - p((n-1)/(y+n))

pH = np.linspace(.5, 13.5, 1000)

n1 = np.linspace(0, .999, 1000)
n2 = np.linspace(1.001, 2, 1000)
pHA = pH_acide(n1)
pHB = pH_basique(n2)

n = HCL.n_dilue(CA, CB, pH=pH)
n_nc = HCL.n(CA, pH=pH)

plt.figure()
i_acide = np.argwhere(pH < 7)
i_basique = np.argwhere(pH > 7)
nc1 = n_nc[i_acide]
nc2 = n_nc[i_basique]
line = plt.plot(nc1, pH[i_acide]-p(1+nc1/y), lw=3, label='correction')
plt.plot(nc2, pH[i_basique]+p(1+nc2/y), lw=3, color=line[0].get_color())

plt.plot(n, pH, label="méthode corrigée")

plt.plot(n_nc, pH, label="article faux")

line = plt.plot(n1, pHA, '--', label="ATTENDU")
plt.plot(n2, pHB, '--', color=line[0].get_color())

plt.xlim(0, 2)
plt.ylim(.5, 13.5)
plt.ylabel("$p$H")
plt.xlabel("$n$")
plt.grid()
plt.legend()
plt.tight_layout()

plt.show()

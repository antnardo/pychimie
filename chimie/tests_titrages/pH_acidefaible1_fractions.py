#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# print(S.melange['Acide test'].Ca(pH = 12))
import numpy as np
import matplotlib.pyplot as plt
from timeit import timeit

from chimie.titrage import PolyAcide, p

pKa1 = 6.3
pKa2 = 10.3
T = PolyAcide([pKa1, pKa2])  # carbonates

# Fractions molaires
pHs = np.linspace(1, 13, 200)
x0, x1, x2 = T.fractions(pH=pHs)

plt.figure(figsize=(6, 4))
plt.plot(pHs, x0, label="H2CO3")
plt.plot(pHs, x1, label="HCO3-")
plt.plot(pHs, x2, label="CO32-")
plt.xlabel("pH")
plt.ylabel("fractions molaires")
plt.grid()
plt.xlim(1, 13)
plt.ylim(0, 1)
plt.legend()
plt.tight_layout()
plt.show(block=False)
plt.savefig("fractions_H2CO3.pdf")

# pH de la solution

CT = T.Ctot(pH=pHs)
CT2 = np.logspace(-9, 4, 200)
pH2 = T.pH_inverse(Ctot=CT2)
pH3 = T.pH(Ctot=CT2)

pCA = p(CT)

plt.figure()
plt.plot(pCA, pHs, label="forward")
plt.plot(p(CT2), pH2, label="inverse")
plt.plot(p(CT2), pH3, label="forward interpolé")

plt.plot(pCA, .5*(pCA+pKa1), 'r--', lw=.5, label=r"\frac12(pC_A+pKa_1)")
plt.plot(pCA, pCA, 'r--', lw=.5)
plt.plot(pCA, 7*np.ones_like(pCA), 'r--', lw=.5)

plt.legend()
plt.xlim(0, 9)
plt.ylim(0, 9)
plt.ylabel("pH")
plt.xlabel("$C_t$ ($n=0$)")
plt.grid()
plt.tight_layout()
plt.show(block=False)
plt.savefig("pH_H2CO3.pdf")

# t1 = timeit("T.Ctot(pH=pHs)", number=100, globals=globals())
# t2 = timeit("T.pH_inverse(Ctot=CT2)", number=100, globals=globals())
# t3 = timeit("T.pH(Ctot=CT2)", number=100, globals=globals())
# print(t1, t2, t3)
# 0.006292375008342788 15.741244045013445 2.2830512060027104

plt.show()

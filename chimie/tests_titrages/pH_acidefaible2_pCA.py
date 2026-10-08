#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# print(S.melange['Acide test'].Ca(pH = 12))
import numpy as np
import matplotlib.pyplot as plt
from timeit import timeit

from chimie.titrage import PolyAcide, p

pKa1 = 4
T = PolyAcide([pKa1])  # carbonates

pHs = np.linspace(1, 13, 200)
CT = T.Ctot(pH=pHs)
pCA = p(CT)

plt.figure(figsize=(5, 5))
plt.plot(pHs, pCA)
plt.xlim(0, 8)
plt.ylim(0, 8)
plt.xlabel("pH")
plt.ylabel("$pC_A$")
plt.grid()
plt.tight_layout()
plt.show(block=False)
plt.savefig("pCA1.pdf")

plt.figure(figsize=(5, 5))
plt.plot(pCA, pHs)

plt.plot(pCA, .5*(pCA+pKa1), 'r--', lw=.5, label=r"\frac12(pC_A+pKa)")
plt.plot(pCA, pCA, 'r--', lw=.5)
plt.plot(pCA, 7*np.ones_like(pCA), 'r--', lw=.5)
plt.xlim(0, 8)
plt.ylim(0, 8)
plt.ylabel("pH")
plt.xlabel("$pC_A$ ($n=0$)")
plt.grid()
plt.tight_layout()
plt.show(block=False)
plt.savefig("pCA2.pdf")


plt.show()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# print(S.melange['Acide test'].Ca(pH = 12))
import numpy as np
import matplotlib.pyplot as plt

from chimie.titrage import PolyAcide, p as plog

T = PolyAcide([6.3, 10.3])  # carbonates

# Titrages par une base forte
CB = 0.1
CTS = [0.01, 0.1, 1]

plt.xlim(0, 3)
plt.ylim(0, 14)
plt.legend()
plt.ylabel("$p$H")
plt.xlabel(r"$n = V \dfrac{C_\text{burette}}{C_\text{acide faible}V_\text{becher}}$")
plt.grid()
plt.tight_layout()
plt.savefig("acidefaible_baseforte_tangentes.pdf")

plt.figure()
CT = 0.1
pH_init = T.pH(Ctot=CT)
pHs = np.linspace(pH_init, 14-plog(CB)-.2, 400)
n = T.n_dilue(CA=CT, CB=CB, pH=pHs)
# plt.plot(n, pHs)
pentes = 1/T.beta_dilue(Ctot=CT, CB=CB, pH=pHs)
for i in range(0, len(n), 1):
    n0 = n[i]
    pH0 = pHs[i]
    p = pentes[i]
    nx = np.array([n0-.5, n0+.5])
    t = T.tangente(nx, n0, pH0, p)
    plt.plot(nx, t, color='k', lw=0.3)
    plt.plot([n0], [pH0], color='k', ms=1, marker='o')
plt.xlim(0, 3)
plt.ylim(0, 14)
plt.ylabel("$p$H")
plt.xlabel("$n$")
plt.grid()
plt.tight_layout()
plt.savefig("verif_pentes_dilue.pdf")

plt.show()

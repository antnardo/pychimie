#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# print(S.melange['Acide test'].Ca(pH = 12))
import numpy as np
import matplotlib.pyplot as plt

from chimie.titrage import PolyAcide, p as plog

pKa1 = 4
pKa2 = 10
T = PolyAcide([pKa1, pKa2])  # carbonates

# Titrages par une base forte
CB = 0.1
CTS = [0.01, 0.1, 1]

# Calculs classiques
def pH1(n):
    return pKa1 - plog(n/(1-n))

def pH2(n):
    return pKa2 - plog((n-1)/(2-n))

def pH3(n, CT):
    return 14 - plog((n-2)/(CB+n*CT)*CT*CB)

pHeq1 = .5*(pKa1+pKa2)

n1 = np.linspace(0, 1, 50)
n2 = np.linspace(1, 2, 50)
pH1_th = pH1(n1)
pH2_th = pH2(n2)

plt.figure(figsize=(6, 4))
ax = plt.gca()

plt.plot([0, 3], [pKa1, pKa1], 'r--')
plt.plot([0, 3], [pKa2, pKa2], 'r--')

plt.plot(n1, pH1_th, 'k--')
plt.plot(n2, pH2_th, 'k--')

plt.plot([1], [pHeq1], 'sk')

for CT in CTS:
    pH_init = T.pH(Ctot=CT)
    pHs = np.linspace(pH_init, 14-plog(CB)-.2, 200)

    n_article = T.n(Ctot=CT, pH=pHs)
    n = T.n_dilue(CA=CT, CB=CB, pH=pHs)
    ni = max(np.ceil(2 * min(n)) / 2, -1)
    nf = min(np.floor(2 * max(n)) / 2, 3)

    lines = plt.plot(n, pHs, label=f"Ctot={CT:.1e} M")
    c = lines[0].get_color()

    pHeq0 = .5*(pKa1+plog(CT))
    pHeq2 = 7+.5*pKa2+.5*np.log(CT*CB/(CT+CB))
    n3 = np.linspace(2, nf, 50)
    pH3_th = pH3(n3, CT)
    plt.plot([0, 2], [pHeq0, pHeq2], 's', color=c)
    plt.plot(n3, pH3_th, '--', color=c)

    plt.plot(n_article, pHs, '--', lw=.7, color=c)

    method = T.draw_tangentes(ax=ax, Ctot=CT, CB=CB, n0=1, delta=0.4,
                              draw_kwargs={'color': c}, draw_points_kwargs={'color': c})
    print(f"CT={CT:.1e}M: n={method[0]}, pH={method[1]}, pente={method[6]}")

    method = T.draw_tangentes(ax=ax, Ctot=CT, CB=CB, n0=2, delta=0.4,
                              draw_kwargs={'color': c}, draw_points_kwargs={'color': c})
    print(f"CT={CT:.1e}M: n={method[0]}, pH={method[1]}, pente={method[6]}")

    EPs = []
    for p in np.arange(ni, nf + 0.5, step=0.5):
        EP = T.pH_EP(CA=CT, CB=CB, p=p)
        EPs.append((p, EP))
        if p < 3:
            plt.plot([0, p], [EP, EP], "o-", lw=.5, color=lines[0].get_color())
            plt.text(-0.05, EP, f"{EP:.1f}", fontsize=10, ha='right', va='center')

    # print(EPs)

plt.xlim(0, 3)
plt.ylim(0, 14)
plt.legend()
plt.ylabel("$p$H")
plt.xlabel(r"$n = V \dfrac{C_\text{burette}}{C_\text{acide faible}V_\text{becher}}$")
plt.grid()
plt.tight_layout()
plt.savefig("acidefaible_baseforte_tangentes.pdf")

plt.show()

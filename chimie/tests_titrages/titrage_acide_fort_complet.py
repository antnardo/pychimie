#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt


def p(c):
    return -np.log10(c)

def c(p):
    return 10**(-p)

PKE = 14
KE = c(PKE)
TOL = KE / 100
LN10 = np.log(10)

def delta(h):
    return KE/h - h

def n(h, CT):
    return 1 + delta(h)/CT

def n_bis(h, CT, y):
    return (1 + delta(h)/CT)/(1 - delta(h)/y/CT)

CA = 0.1
CB = 0.4
y = CB/CA

def pH_acide(n):
    return p(CB) + p((1-n)/(y+n))

def pH_basique(n):
    return 14 - p(CB) - p((n-1)/(y+n))

pH = np.linspace(.5, 13.5, 1000)
h = c(pH)
n_calc = n(h, CA)
n_calc2 = n_bis(h, CA, y)

nt = np.array([0, 2])
n1 = np.linspace(0, .999, 1000)
n2 = np.linspace(1.001, 2, 1000)
pHA = pH_acide(n1)
pHB = pH_basique(n2)

def pente_acide(n):
    return (y+1)/(1-n)/(y+n)/LN10

def pente_basique(n):
    return (y+1)/(n-1)/(y+n)/LN10

def tangente_acide(n, n1):
    return pente_acide(n1)*(n-n1) + pH_acide(n1)

def tangente_basique(n, n2):
    return pente_basique(n2)*(n-n2) + pH_basique(n2)

def get_tangentes(p):
    n = np.linspace(0, .9999, 1000)
    p1 = pente_acide(n)
    n1 = np.interp(p, p1, n)
    n = np.linspace(2, 1.0001, 1000)
    p2 = pente_basique(n)
    n2 = np.interp(p, p2, n)
    return n1, n2

def methode(n1, n2):
    t1 = tangente_acide(nt, n1)
    t2 = tangente_basique(nt, n2)
    te = (t1+t2)/2
    return t1, t2, te, te.mean()

def draw_tangentes(ax, n1, n2):
    t1, t2, te, tE = methode(n1, n2)
    lines = ax.plot(nt, t1, '--', lw=.5)
    ax.plot(nt, t2, '--', color=lines[0].get_color(), lw=.5)
    ax.plot(nt, te, '--', color=lines[0].get_color(), lw=.5)
    ax.plot([1], [tE], 'o', color=lines[0].get_color())
    return tE

plt.figure()
p1 = pente_acide(n1)
p2 = pente_basique(n2)
plt.plot(n1, p1)
plt.plot(n2, p2)

plt.xlim(0, 2)
plt.ylabel("pente")
plt.xlabel("$n$")
plt.grid()
plt.tight_layout()

pmin = max(min(p1), min(p2))
pmax = 15
d = pmax-pmin
ps = np.linspace(pmin+d*.1, pmin + d, 5)

plt.figure()
ax = plt.gca()
plt.plot(n_calc, pH+p(1+n_calc/y), lw=3, label='correction basique')
plt.plot(n_calc, pH-p(1+n_calc/y), lw=3, label='correction acide')
plt.plot(n_calc, pH, label="article faux")
plt.plot(n_calc2, pH, label="correction complete")

line = plt.plot(n1, pHA, '--', label="ATTENDU")
plt.plot(n2, pHB, '--', color=line[0].get_color())

for pp in ps:
    n1, n2 = get_tangentes(pp)
    te = draw_tangentes(ax, n1, n2)
    print(pp, n1, n2, te)

plt.xlim(0, 2)
plt.ylim(.5, 13.5)
plt.ylabel("$p$H")
plt.xlabel("$n$")
plt.grid()
plt.legend()
plt.tight_layout()


equivalences = []
p_equiv = np.linspace(1, 15, 1000)
n1_equiv = []
n2_equiv = []
for pp in p_equiv:
    n1, n2 = get_tangentes(pp)
    t1, t2, te, tE = methode(n1, n2)
    equivalences.append(tE)
    n1_equiv.append(n1)
    n2_equiv.append(n2)

equivalences = np.array(equivalences)
n1_equiv = np.array(n1_equiv)
n2_equiv = np.array(n2_equiv)

delta1 = -p_equiv*(2-n1_equiv-n2_equiv)/2
delta2 = -np.log10((y+n1_equiv)/(y+n2_equiv))

plt.figure()
plt.plot(p_equiv, delta1+delta2, label="th", lw=2)
plt.plot(p_equiv, delta1, '--')
plt.plot(p_equiv, delta2/2, '--')
plt.plot(p_equiv, 7-equivalences, label="7-$p$H")
plt.plot(p_equiv, 1-n1_equiv, label="position avant")
plt.plot(p_equiv, n2_equiv-1, label="position après")
plt.legend()
plt.xlabel("pente")
plt.ylabel("écart au $p$H équivalence")
plt.title("Méthode des tangeantes approchée")
plt.grid()
plt.tight_layout()

# à partir de ces valeurs, pour un point en particulier
pH_demiequivalence = np.interp(0.5, n_calc, pH)
pH_equivalence = np.interp(1, n_calc, pH)
print(pH_demiequivalence, pH_equivalence)

plt.show()


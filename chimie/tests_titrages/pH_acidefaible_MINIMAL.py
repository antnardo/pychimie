#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from itertools import accumulate
from operator import mul


def p(c):
    return -np.log10(c)

def c(p):
    return 10**(-p)

PKE = 14
KE = c(PKE)
TOL = KE / 100


pKa1 = 6.3
pKa2 = 10.3

pKa = [pKa1, pKa2]
N = len(pKa)
Ka = [c(pK) for pK in pKa]
k = list(accumulate(Ka, mul))
N = len(Ka)

def x0(h):
    return 1/(1+sum(k[j] / h**(j+1) for j in range(N)))

def x(h, j):
    return x0(h) * k[j-1] / h**j

def Y1(h):
    return sum([j*x(h, j) for j in range(1, N+1)])

def delta(h):
    return KE/h - h

def n(h, CT):
    return Y1(h) + delta(h)/CT

pH = np.linspace(1, 13, 200)
h = c(pH)

CT = 0.1
n_calc = n(h, CT)

plt.figure()
plt.plot(n_calc, pH)
plt.show()

# à partir de ces valeurs, pour un point en particulier
pH_demiequivalence = np.interp(0.5, n_calc, pH)
pH_equivalence = np.interp(1, n_calc, pH)
print(pH_demiequivalence, pH_equivalence)

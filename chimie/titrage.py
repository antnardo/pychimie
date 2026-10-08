"""
11/2024

Kalka,H. PolyproticAcids and Beyond—An Algebraic Approach.
Chemistry 2021,3,454–508.
https:// doi.org/10.3390/chemistry3020034

@author: elessar
"""
import numpy as np
from scipy.optimize import brentq
from itertools import accumulate
from functools import cache
import operator
from collections.abc import Iterable
from matplotlib.axes import Axes


def p(c):
    return -np.log10(c)

def c(p):
    return 10**(-p)

PKE = 14
KE = c(PKE)
TOL = KE / 100
PHTOL = 0.01
LN10 = np.log(10)


class PolyAcide:
    draw_kwargs_default = dict(lw=.5, color="red", ls='--')
    draw_points_kwargs_default = dict(color="red", ls='', marker='x')

    def __init__(self, pKas) -> None:
        """liste des pKa, rerangés par ordre croissant
        de sorte que l'indice 0 correspond à la constante de l'équation
        H_NA = H+ + H_(N-1)A-
        etc
        les k_j sont les produits de ces constantes, donc les constantes des équations
        H_NA = j H+ + H_(N-j)A^-j
        """
        self.K = [c(pKa) for pKa in sorted(pKas)]
        self.N = len(self.K)
        self.k = [1] + list(accumulate(self.K, operator.mul))

    @staticmethod
    def _h(*args, **kwargs):
        if len(args) == 1:
            return args[0]  # h par defaut
        h = kwargs.get("h", None)
        if h is not None:
            return h
        pH = kwargs.get("pH", None)
        if pH is not None:
            return c(pH)
        pOH = kwargs.get("pOH", None)
        if pOH is not None:
            return c(PKE - pOH)
        raise ValueError

    @staticmethod
    @cache
    def _pH(pHmin, pHmax, npoints):
        return np.linspace(pHmin, pHmax, npoints)

    @staticmethod
    def w(h):
        """variable intermédiaire liée à la conservation de la charge"""
        return KE / h - h

    def a(self, *args, **kwargs):
        """a1, a2, a3... a_N normalisés par x0"""
        h = self._h(*args, **kwargs)
        return [self.k[j] / h**j for j in range(1, self.N + 1)]

    def x0(self, *args, **kwargs):
        """fraction molaire de l'acide le plus fort
        par CM et par les LAM :
        c_k ="""
        return 1 / (1 + sum(self.a(*args, **kwargs)))

    def x(self, j, *args, **kwargs):
        """fraction molaire de l'acide j"""
        h = self._h(*args, **kwargs)
        x0 = self.x0(h)
        if j == 0:
            return x0
        return x0 * self.k[j] / h**j  # self.a(h)[j-1]

    def fractions(self, *args, **kwargs):
        h = self._h(*args, **kwargs)
        return [self.x(j, h) for j in range(self.N + 1)]

    def Y(self, k, *args, **kwargs):
        """Moments des fractions d'ordre k"""
        h = self._h(*args, **kwargs)
        x0 = self.x0(h)
        a = self.a(h)
        return x0 * sum(j**k * a[j-1] for j in range(1, self.N + 1))

    def n(self, Ctot, *args, **kwargs):
        h = self._h(*args, **kwargs)
        w = self.w(h)
        Y = self.Y(1, h)
        return Y + w / Ctot

    def n_dilue(self, CA, CB, *args, **kwargs):
        h = self._h(*args, **kwargs)
        w = self.w(h)
        Y = self.Y(1, h)
        return (Y + w/CA) / (1 - w/CB)

    def Cb(self, Ctot, *args, **kwargs):
        return Ctot * self.n(Ctot, *args, **kwargs)

    def Ctot(self, Ca=0, Cb=0, **kwargs):
        h = self._h(**kwargs)
        w = self.w(h)
        Y = self.Y(1, h)
        return (Cb - Ca - w) / Y

    def beta(self, Ctot, *args, **kwargs):
        """dn/dpH"""
        h = self._h(*args, **kwargs)
        w = self.w(h)
        Y1, Y2 = [self.Y(k, h) for k in [1, 2]]
        return LN10 * (Y2 - Y1**2 + (w + 2 * h) / Ctot)

    def beta_dilue(self, Ctot, CB, *args, **kwargs):
        h = self._h(*args, **kwargs)
        w = self.w(h)
        Y1 = self.Y(h, 1)
        dilution = 1 / (1-w/CB)
        return (self.beta(Ctot, h) + LN10*self.n(Ctot, h)*(w + 2*h)/CB) * dilution

    def beta2(self, Ctot, *args, **kwargs):
        """d^2n/dpH^2"""
        h = self._h(*args, **kwargs)
        w = self.w(h)
        Y1, Y2, Y3 = [self.Y(k, h) for k in [1, 2, 3]]
        return LN10**2 * (Y3 - 3 * Y1 * Y2 + 2 * Y1**3 + w / Ctot)

    def methode_tangentes(self, Ctot, CB=None, n0=1, delta=0.5, npoints=200, dpH=0.1, dilue=True):
        """recherche autour de n0 à environ delta=0.5 autour
        n0 est a priori entier (sinon pas de garantie)
        """
        if CB is None:
            dilue = False
        # Définition des points de calcul
        n1 = n0 - delta
        n2_approx = n0 + delta
        if dilue:
            pH1 = self.pH_EP(Ctot, p=n1, CB=CB)
            pH2_approx = self.pH_EP(Ctot, p=n2_approx, CB=CB)
        else:
            pH1 = self.pH(Ctot, Cb=n1*Ctot)
            pH2_approx = self.pH(Ctot, Cb=n2_approx*Ctot)
        pH = self._pH(pH1-dpH, pH2_approx+dpH, npoints)
        # Calculs de la courbe
        if dilue:
            p = 1/self.beta_dilue(Ctot, CB, pH=pH)
            n = self.n_dilue(Ctot, CB, pH=pH)
        else:
            p = 1/self.beta(Ctot, pH=pH)
            n = self.n(Ctot, pH=pH)
        # Calcul de la pente
        pT = np.interp(n1, n, p)
        # On cherche le second point précisément, ce n'est croissant qu'après
        # Ce serait plus efficace de chercher où la dérivée devient positive...
        iafter = np.concatenate(np.argwhere(n > n0))[::-1]
        n2 = np.interp(pT, p[iafter], n[iafter])
        if dilue:
            pH2 = self.pH_EP(Ctot, p=n2, CB=CB)
        else:
            pH2 = self.pH(Ctot, Cb=n2*Ctot)
        # On se met au milieu
        nmoy = (n1+n2)/2
        pHmoy = (pH1+pH2)/2
        ecart = nmoy + (pH-pHmoy)/pT - n
        imax = np.argmax(ecart)  # on évite l'éventuel croisement suivant
        pH_methode = np.interp(0, ecart[:imax], pH[:imax])
        n_methode = np.interp(pH_methode, pH, n)
        return n_methode, pH_methode, n1, n2, pH1, pH2, pT

    def draw_tangentes(self, ax: Axes, *args, **kwargs):
        draw_kwargs = dict(self.draw_kwargs_default)
        draw_kwargs.update(kwargs.pop("draw_kwargs", dict()))
        draw_points_kwargs = dict(self.draw_points_kwargs_default)
        draw_points_kwargs.update(kwargs.pop("draw_points_kwargs", dict()))
        method_return = kwargs.pop("method_points", None)
        if method_return is None:
            n_methode, pH_methode, n1, n2, pH1, pH2, pT = self.methode_tangentes(*args, **kwargs)
        else:
            n_methode, pH_methode, n1, n2, pH1, pH2, pT = method_return
        nx = kwargs.pop("nx", np.array([n1-.1, n2+.1]))

        t1 = self.tangente(nx, n1, pH1, pT)
        t2 = self.tangente(nx, n2, pH2, pT)
        tmoy = (t1+t2)/2
        ax.plot(nx, t1, **draw_kwargs)
        ax.plot(nx, t2, **draw_kwargs)
        ax.plot(nx, tmoy, **draw_kwargs)
        ax.plot([n_methode, n1, n2], [pH_methode, pH1, pH2], **draw_points_kwargs)
        return n_methode, pH_methode, n1, n2, pH1, pH2, pT

    @staticmethod
    def tangente(n, n0, pH0, pente):
        return pente*(n-n0) + pH0

    @staticmethod
    def _pHmin(pH, CB):
        if pH is not None:
            return pH
        return p(CB)+PHTOL

    @staticmethod
    def _pHmax(pH, CB):
        if pH is not None:
            return pH
        return p(np.sqrt(CB**2/4+KE)-CB/2)-PHTOL

    def pH_EP(self, CA, CB, p=1, pHmin=None, pHmax=None, npoints=200):
        """
        renvoit le pH de l'équivalence p pour un titrage (avec dilution)
        interpolation pour trouver les intersections de n dans un certain intervalle
        p = 1 est l'équivalence classique d'un titrage par une base forte
        p = -1 par un acide fort
        pour p demi entier, cela correspond aux zéros de la fonction d^2n/dpH^2
            donc l'interpolation est très bonne à ces endroits (la dérivée seconde y est nulle...)
        """
        pHmin = self._pHmin(pHmin, CB)
        pHmax = self._pHmax(pHmax, CB)
        pH = self._pH(pHmin, pHmax, npoints)
        n = self.n_dilue(CA=CA, CB=CB, pH=pH)  # fonction croissante de pH
        return np.interp(p, n, pH)

    def pH_melange(self, Ctot, CA, CB, pHmin=0, pHmax=14, npoints=200):
        """
        renvoit le pH d'un solution acide fiable + base forte + acide fort
        C donnés sans dilution (les C de la solution)
        """
        pH = self._pH(pHmin, pHmax, npoints)
        n = self.n(Ctot=Ctot, pH=pH)  # fonction croissante de pH
        p = (CB-CA)/Ctot
        return np.interp(p, n, pH)

    def pH(self, Ctot, Ca=0, Cb=0, **kwargs):
        if isinstance(Ctot, Iterable):
            return [self.pH(c, Ca, Cb, **kwargs) for c in Ctot]
        return self.pH_melange(Ctot, Ca, Cb, **kwargs)

    def pH_inverse(self, Ctot, Ca=0, Cb=0):
        """ calcul inverse """
        if isinstance(Ctot, Iterable):
            return [self.pH_inverse(c, Ca, Cb) for c in Ctot]
        n = (Cb - Ca) / Ctot

        def equa(h):
            return sum(
                self.k[j] * h ** (self.N - j) * (h**2 + (n - j) * Ctot * h - KE)
                for j in range(self.N + 1)
            )

        try:
            h = brentq(equa, TOL, 100, xtol=TOL)
        except ValueError:
            h = np.NaN
        return -np.log10(h)

def AcideFort():
    return PolyAcide([-5])

def BaseForte():
    return PolyAcide([PKE+5])

class Melange:
    def __init__(self) -> None:
        self.acides = []
        self.concentrations = []
        self.fractions = []
        self.CT = 0

    def __len__(self):
        return len(self.acides)

    def dilution(self, x):
        for i in range(len(self)):
            self.concentrations[i] /= x

    def add_polyacide(self, P, C):
        self.acides.append(P)
        self.concentrations.append(C)
        self._set_fractions()

    def _set_fractions(self):
        self.CT = sum(self.concentrations)
        self.fractions = [x/self.CT for x in self.concentrations]

    @staticmethod
    def _h(*args, **kwargs):
        return PolyAcide._h(*args, **kwargs)

    @staticmethod
    def w(h):
        """variable intermédiaire liée à la conservation de la charge"""
        return PolyAcide.w(h)

    def n(self, *args, **kwargs):
        h = self._h(*args, **kwargs)
        w = self.w(h)
        Y = sum([x*A.Y(1, h) for x, A in zip(self.fractions, self.acides)])
        return Y + w / self.CT


H2CO3 = PolyAcide([6.37, 10.30])

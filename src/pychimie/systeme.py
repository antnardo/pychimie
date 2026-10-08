"""Un système qui évolue par une réaction : avancement, Q(ξ), ΔrG(ξ), G(ξ), état final.

Un :class:`Systeme` est une réaction, des quantités initiales (mol), d'éventuels
constituants inactifs, et des conditions : T (K), P (bar) imposée, V (L) de la
solution, ou V_gaz (L) si le gaz occupe un volume fixé au lieu d'être à pression
fixée. Tout le reste en découle : les quantités nᵢ(ξ) = nᵢ,₀ + νᵢξ, les activités,
Q(ξ), ΔrG(ξ) = ΔrG° + RT ln Q, et l'état final.

**État final.** ΔrG(ξ) = dG/dξ croît avec ξ (G est convexe pour des mélanges idéaux) :
l'état final est le minimum de G sur [ξ_min, ξ_max]. Ou bien ΔrG s'annule à
l'intérieur — c'est l'équilibre, Q = K° —, ou bien il garde son signe jusqu'à une
borne, où disparaît un constituant **seul dans sa phase** (solide ou liquide pur, gaz
seul à pression fixée) dont l'activité ne tend pas vers zéro : c'est la rupture
d'équilibre, et la transformation est totale. L'ancien ``Reaction`` ne testait que le
réactif limitant (``is_rupture_min`` n'était jamais appelé) ; ici la même règle vaut aux
deux bornes, parce qu'elle découle de la seule monotonie de ΔrG.

**Précision.** On cherche la racine depuis la borne la plus proche, en paramétrant
ξ = borne ∓ u et en calculant les quantités comme m + νu, où m vaut exactement zéro
pour les constituants qui s'annulent à la borne. Avec K° = 10²⁴ (2 SO₂ + O₂ à 298 K),
il reste 10⁻¹² mol de SO₂, que ξ = 1 − 10⁻¹² ne distinguerait pas de zéro : ces
quantités-là sont exactes à la précision de la machine, et u se cherche en logarithme.

**G(ξ).** Plutôt que d'intégrer ΔrG (l'ancien ``compute_G`` appelait ``quad`` et
plantait sur le tuple qu'il rend), G se calcule en forme close :

    G(ξ) − G(0) = ξ ΔrG° + RT [Φ(ξ) − Φ(0)],   Φ = Σᵢ nᵢ ln aᵢ − Σ_dilués nᵢ,

dont la dérivée est exactement ΔrG(ξ) : pour une activité en fraction molaire,
Σⱼ nⱼ ∂ln xⱼ/∂nᵢ = 0 (Gibbs-Duhem) ; pour une activité proportionnelle à nᵢ (soluté,
gaz à volume fixé), le terme −nᵢ compense la dérivée de nᵢ ln nᵢ. À volume de gaz fixé,
la fonction minimale est l'énergie libre F ; G(ξ) rend alors F(ξ) − F(0), de même
dérivée.

Les quantités sont en mol, P en bar (P° = 1 bar), V en L (c° = 1 mol/L).
"""

import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.constants import R
from scipy.optimize import brentq, minimize_scalar

from pychimie._calcul import Flottants, sortie
from pychimie.equation import Constituant, indice
from pychimie.reaction import Reaction

__all__ = ["EtatFinal", "RuptureEquilibre", "Systeme"]

P_STANDARD = 1e5  # Pa ; les pressions de pychimie sont en bar, P° = 1 bar.
C_STANDARD = 1.0  # mol/L
_EPS = np.finfo(float).eps

Tableau = NDArray[np.float64]


class RuptureEquilibre(ValueError):  # noqa: N818 — nom français, sans suffixe « Error »
    """L'état final n'est pas un équilibre : une phase a disparu avant que Q n'atteigne K°."""


@dataclass(frozen=True, slots=True)
class EtatFinal:
    """L'état final d'un système : équilibre (Q = K°), ou rupture d'équilibre.

    ``constituants`` et ``quantites`` vont ensemble : ceux de la réaction, puis les
    inactifs. Quand l'état n'est pas un équilibre, ``disparus`` nomme les constituants
    de la réaction absents de l'état final.
    """

    xi: float
    xi_min: float
    xi_max: float
    constituants: tuple[Constituant, ...]
    quantites: tuple[float, ...]
    equilibre: bool
    disparus: tuple[str, ...]
    Q: float
    K0: float

    def n(self, nom: str) -> float:
        """La quantité finale de ``nom`` (mol)."""
        return self.quantites[indice(self.constituants, nom)]

    def fraction(self, nom: str) -> float:
        """La fraction molaire finale de ``nom`` dans sa phase (gaz, mélange, corps pur)."""
        c = self.constituants[indice(self.constituants, nom)]
        if c.melange is None:
            raise ValueError(f"{c.nom} est un soluté : pas de fraction molaire, voir n()")
        total = sum(
            q
            for k, q in zip(self.constituants, self.quantites, strict=True)
            if k.melange == c.melange
        )
        return self.n(nom) / total if total > 0 else math.nan

    @property
    def taux_avancement(self) -> float:
        """τ = ξ_f / ξ_max (nan si aucun réactif n'est présent)."""
        return self.xi / self.xi_max if self.xi_max != 0 else math.nan

    @property
    def sens(self) -> int:
        """+1 si la réaction a avancé dans le sens direct, −1 dans le sens inverse, 0 sinon."""
        return int(np.sign(self.xi))


@dataclass(frozen=True, slots=True, init=False)
class Systeme:
    """Une réaction, des quantités initiales et des conditions (voir le module)."""

    reaction: Reaction
    initial: tuple[float, ...]
    inertes: tuple[tuple[Constituant, float], ...]
    T: float | None
    P: float
    V: float | None
    V_gaz: float | None

    def __init__(
        self,
        reaction: Reaction,
        initial: Mapping[str, float],
        *,
        T: float | None = None,
        P: float = 1.0,
        V: float | None = None,
        V_gaz: float | None = None,
        inertes: Mapping[str, float] | None = None,
    ) -> None:
        """``Systeme(nh3, {"N2": 1, "H2": 3}, T=700, P=200)``.

        Un constituant absent de ``initial`` vaut 0 mol. Les inactifs s'écrivent avec
        leur phase : ``inertes={"Ar(g)": 1}``. ``P`` est la pression imposée ;
        ``V_gaz``, si on le donne, fixe le volume du gaz à la place (P en découle).
        """
        constituants = reaction.equation.constituants
        quantites = [0.0] * len(constituants)
        vus: set[int] = set()
        for nom, n in initial.items():
            i = indice(constituants, nom)
            if i in vus:
                raise ValueError(f"{constituants[i].nom} est donné deux fois")
            vus.add(i)
            quantites[i] = _quantite(nom, n)
        liste_inertes = []
        for nom, n in (inertes or {}).items():
            c = Constituant.lire(nom)
            if c.nom in reaction.equation.noms:
                raise ValueError(f"{c.nom} participe à la réaction : le mettre dans initial")
            liste_inertes.append((c, _quantite(nom, n)))
        for nom, valeur in (("T", T), ("P", P), ("V", V), ("V_gaz", V_gaz)):
            if valeur is not None and valeur <= 0:
                raise ValueError(f"{nom} doit être strictement positif")
        tous = (*constituants, *(c for c, _ in liste_inertes))
        if V is None and any(c.etat == "aq" for c in tous):
            raise ValueError("des solutés : donner le volume V de la solution (L)")
        if V_gaz is not None and T is None:
            raise ValueError("à volume de gaz fixé, l'activité nRT/(VP°) demande T")
        object.__setattr__(self, "reaction", reaction)
        object.__setattr__(self, "initial", tuple(quantites))
        object.__setattr__(self, "inertes", tuple(liste_inertes))
        object.__setattr__(self, "T", None if T is None else float(T))
        object.__setattr__(self, "P", float(P))
        object.__setattr__(self, "V", None if V is None else float(V))
        object.__setattr__(self, "V_gaz", None if V_gaz is None else float(V_gaz))

    def avec(self, **changements: Any) -> "Systeme":
        """Le même système, avec d'autres conditions ou d'autres quantités initiales.

        ``initial`` et ``inertes`` remplacent la composition entière :
        ``s.avec(initial={"N2": x, "H2": 1 - x})``, ``s.avec(T=800)``.
        """
        actuel: dict[str, Any] = {
            "initial": dict(zip(self.reaction.equation.noms, self.initial, strict=True)),
            "T": self.T,
            "P": self.P,
            "V": self.V,
            "V_gaz": self.V_gaz,
            "inertes": {c.nom: n for c, n in self.inertes},
        }
        inconnus = set(changements) - set(actuel)
        if inconnus:
            raise TypeError(f"paramètres inconnus : {', '.join(sorted(inconnus))}")
        actuel.update(changements)
        return Systeme(self.reaction, **actuel)

    # --- bornes et grandeurs ------------------------------------------------------------

    @property
    def constituants(self) -> tuple[Constituant, ...]:
        """Ceux de la réaction, puis les inactifs."""
        return (*self.reaction.equation.constituants, *(c for c, _ in self.inertes))

    @property
    def xi_max(self) -> float:
        """L'avancement où disparaît le réactif limitant."""
        nu, n0 = self._nu[: self._r], self._n0[: self._r]
        return float(np.min(n0[nu < 0] / -nu[nu < 0])) + 0.0

    @property
    def xi_min(self) -> float:
        """L'avancement (négatif ou nul) où disparaît le produit limitant."""
        nu, n0 = self._nu[: self._r], self._n0[: self._r]
        return float(np.max(-n0[nu > 0] / nu[nu > 0])) + 0.0

    @property
    def K0(self) -> float:
        """K° à la température du système."""
        return float(self.reaction.K0(self.T))

    @property
    def DrG0(self) -> float:
        """ΔrG° à la température du système (J/mol)."""
        return float(self.reaction.DrG0(self._T_requis("ΔrG°")))

    def avancements(
        self, n: int = 2000, *, debut: float | None = None, fin: float | None = None
    ) -> Tableau:
        """``n`` avancements entre ξ_min et ξ_max, bornes exclues (Q y vaut 0 ou ∞).

        Comme l'ancien ``liste_avancement`` : un pas (fin − début)/n est retiré à chaque
        bout. ``debut=0`` part de l'état initial.
        """
        a = self.xi_min if debut is None else debut
        b = self.xi_max if fin is None else fin
        delta = (b - a) / n
        return np.linspace(a + delta, b - delta, n)

    def n(self, nom: str, xi: ArrayLike) -> Flottants:
        """nᵢ(ξ) = nᵢ,₀ + νᵢξ (mol)."""
        x, forme = _grille(xi)
        return sortie(self._quantites(x)[indice(self.constituants, nom)].reshape(forme))

    def activite(self, nom: str, xi: ArrayLike) -> Flottants:
        """aᵢ(ξ)."""
        x, forme = _grille(xi)
        ln_a = self._ln_activites(self._quantites(x))[indice(self.constituants, nom)]
        return sortie(np.exp(ln_a).reshape(forme))

    def pression(self, xi: ArrayLike) -> Flottants:
        """La pression (bar) : P imposée, ou nRT/V_gaz à volume de gaz fixé."""
        x, forme = _grille(xi)
        if self.V_gaz is None:
            return sortie(np.full(forme, self.P))
        gaz = [i for i, c in enumerate(self.constituants) if c.etat == "g"]
        n_gaz = self._quantites(x)[gaz].sum(axis=0)
        T = self._T_requis("la pression")
        return sortie((n_gaz * R * T / (self.V_gaz * 1e-3) / P_STANDARD).reshape(forme))

    def Q(self, xi: ArrayLike) -> Flottants:
        """Le quotient réactionnel Q(ξ) = Π aᵢ^νᵢ."""
        x, forme = _grille(xi)
        with np.errstate(over="ignore"):
            return sortie(np.exp(self._ln_Q(self._quantites(x))).reshape(forme))

    def DrG(self, xi: ArrayLike) -> Flottants:
        """ΔrG(ξ) = ΔrG° + RT ln Q(ξ) (J/mol) : négatif, la réaction avance."""
        T = self._T_requis("ΔrG(ξ)")
        x, forme = _grille(xi)
        return sortie((self.DrG0 + R * T * self._ln_Q(self._quantites(x))).reshape(forme))

    def G(self, xi: ArrayLike) -> Flottants:
        """G(ξ) − G(0) (J), à T et P fixées ; F(ξ) − F(0) à volume de gaz fixé."""
        T = self._T_requis("G(ξ)")
        x, forme = _grille(xi)
        phi = self._phi(self._quantites(x)) - self._phi(self._n0[:, None])
        return sortie((x * self.DrG0 + R * T * phi).reshape(forme))

    def minimum_G(self) -> tuple[float, float]:
        """(ξ, G) au minimum de G, cherché directement, sans passer par Q = K°.

        C'est l'autre lecture de l'état final : le minimum de G est l'équilibre, ou la
        borne d'une rupture. La précision est celle d'une minimisation, ~10⁻⁸ en ξ.
        """
        a, b = self.xi_min, self.xi_max
        if b <= a:
            return 0.0, 0.0
        res = minimize_scalar(
            lambda x: float(self.G(x)),
            bounds=(a, b),
            method="bounded",
            options={"xatol": 1e-12 * (b - a), "maxiter": 2000},
        )
        return float(res.x), float(res.fun)

    # --- état final ---------------------------------------------------------------------

    def etat_final(self) -> EtatFinal:
        """L'équilibre, ou la rupture d'équilibre (voir le module)."""
        ln_K = float(self.reaction.ln_K0(self.T))
        a, b = self.xi_min, self.xi_max
        if b <= a:  # un réactif et un produit absents : rien ne peut avancer
            return self._etat(0.0, self._n0, ln_K, equilibre=False)
        milieu = 0.5 * (a + b)
        ecart = self._ln_Q(self._quantites(np.array([milieu])))[0] - ln_K
        if ecart == 0:
            return self._etat(
                milieu, self._quantites(np.array([milieu]))[:, 0], ln_K, equilibre=True
            )
        # cote = +1 : la racine, si elle existe, est entre le milieu et ξ_max.
        cote = 1.0 if ecart < 0 else -1.0
        borne = b if cote > 0 else a
        m = self._ancrage(borne)
        d = -cote * self._nu  # quantités m + d·u à ξ = borne − cote·u
        U = 0.5 * (b - a)

        def g(u: float) -> float:
            return float(self._ln_Q((m + d * u)[:, None])[0]) - ln_K

        g0 = g(0.0)
        if cote * g0 <= 0:  # ΔrG garde son signe jusqu'à la borne : rupture
            return self._etat(borne, m, ln_K, equilibre=g0 == 0)

        def h(s: float) -> float:
            return g(U * math.exp(s))

        s_bas = -25.0
        while cote * h(s_bas) <= 0:
            s_bas -= 25.0
            if U * math.exp(s_bas) == 0:  # K° hors des flottants : équilibre sur la borne
                return self._etat(borne, m, ln_K, equilibre=True)
        s = 0.0 if cote * h(0.0) >= 0 else brentq(h, s_bas, 0.0, xtol=1e-15, rtol=4 * _EPS)
        u = U * math.exp(s)
        return self._etat(borne - cote * u, m + d * u, ln_K, equilibre=True)

    def xi_eq(self) -> float:
        """ξ à l'équilibre ; :class:`RuptureEquilibre` si l'état final n'en est pas un."""
        etat = self.etat_final()
        if etat.equilibre:
            return etat.xi
        if etat.xi_min == etat.xi_max:
            raise RuptureEquilibre(
                f"aucune évolution possible : {', '.join(etat.disparus)} absents"
            )
        raise RuptureEquilibre(
            f"pas d'équilibre : {', '.join(etat.disparus)} disparaît avant que Q n'atteigne K° "
            f"(Q = {etat.Q:.4g}, K° = {etat.K0:.4g}) ; rupture d'équilibre, l'état final est "
            f"ξ = {etat.xi:.6g} mol (voir etat_final())"
        )

    # --- calcul -------------------------------------------------------------------------

    @property
    def _r(self) -> int:
        return len(self.initial)

    @property
    def _nu(self) -> Tableau:
        return np.array([*self.reaction.equation.nu, *(0.0 for _ in self.inertes)])

    @property
    def _n0(self) -> Tableau:
        return np.array([*self.initial, *(n for _, n in self.inertes)])

    def _T_requis(self, quoi: str) -> float:
        if self.T is None:
            raise ValueError(f"{quoi} demande la température : Systeme(..., T=...)")
        return self.T

    def _quantites(self, xi: Tableau) -> Tableau:
        """nᵢ(ξ), une colonne par avancement ; un résidu d'arrondi négatif devient 0."""
        n0, nu = self._n0[:, None], self._nu[:, None]
        n = n0 + nu * xi
        arrondi = (n < 0) & (n > -8 * _EPS * (n0 + np.abs(nu * xi)))
        return np.where(arrondi, 0.0, n)

    def _ancrage(self, borne: float) -> Tableau:
        """Les quantités à une borne, exactement nulles pour ceux qui s'y annulent."""
        n = self._n0 + self._nu * borne
        echelle = self._n0 + np.abs(self._nu * borne)
        return np.where(np.abs(n) <= 8 * _EPS * echelle, 0.0, n)

    def _melanges(self) -> list[tuple[list[int], float]]:
        """Les phases à fraction molaire, avec ln du facteur d'activité (P/P° pour le gaz)."""
        groupes: dict[str, list[int]] = {}
        for i, c in enumerate(self.constituants):
            if c.melange is not None and not (c.etat == "g" and self.V_gaz is not None):
                groupes.setdefault(c.melange, []).append(i)
        return [(idx, math.log(self.P) if cle == "g" else 0.0) for cle, idx in groupes.items()]

    def _dilues(self) -> list[tuple[int, float]]:
        """Solutés (a = c/c°) et gaz à volume fixé (a = nRT/VP°), avec ln du facteur."""
        dilues = []
        for i, c in enumerate(self.constituants):
            if c.etat == "aq":
                V = (
                    self.V if self.V is not None else math.nan
                )  # V est exigé dès qu'il y a un soluté
                dilues.append((i, -math.log(V * C_STANDARD)))
            elif c.etat == "g" and self.V_gaz is not None:
                facteur = R * self._T_requis("le gaz") / (self.V_gaz * 1e-3 * P_STANDARD)
                dilues.append((i, math.log(facteur)))
        return dilues

    def _ln_activites(self, n: Tableau) -> Tableau:
        """ln aᵢ pour chaque colonne de quantités.

        Une phase entièrement vide (à une borne) a des fractions 0/0 : leur limite est
        le rapport des pentes νᵢ/Σν des constituants qui s'y annulent ensemble — 1 pour
        un constituant seul dans sa phase, d'où la rupture d'équilibre.
        """
        ln_a = np.zeros_like(n)
        with np.errstate(divide="ignore", invalid="ignore"):
            for idx, ln_facteur in self._melanges():
                bloc = n[idx]
                total = bloc.sum(axis=0)
                ln_x = np.log(bloc) - np.log(total)
                vides = total == 0
                if np.any(vides):
                    pentes = self._nu[idx]
                    ln_x[:, vides] = np.log(pentes / pentes.sum())[:, None]
                ln_a[idx] = ln_x + ln_facteur
            for i, ln_facteur in self._dilues():
                ln_a[i] = np.log(n[i]) + ln_facteur
        return ln_a

    def _ln_Q(self, n: Tableau) -> Tableau:
        r = self._r
        return (self._nu[:r, None] * self._ln_activites(n)[:r]).sum(axis=0)

    def _phi(self, n: Tableau) -> Tableau:
        """Φ = Σ nᵢ ln aᵢ − Σ_dilués nᵢ, avec 0·ln 0 = 0."""
        with np.errstate(invalid="ignore"):
            termes = np.where(n > 0, n * self._ln_activites(n), 0.0)
        for i, _ in self._dilues():
            termes[i] -= n[i]
        return termes.sum(axis=0)

    def _etat(self, xi: float, n: Tableau, ln_K: float, *, equilibre: bool) -> EtatFinal:
        r = self._r
        noms = self.reaction.equation.noms
        disparus = (
            () if equilibre else tuple(nom for nom, q in zip(noms, n[:r], strict=True) if q == 0)
        )
        with np.errstate(over="ignore", invalid="ignore"):
            Q = float(np.exp(self._ln_Q(n[:, None])[0]))
            K0 = float(np.exp(ln_K))
        return EtatFinal(
            xi=float(xi) + 0.0,
            xi_min=self.xi_min,
            xi_max=self.xi_max,
            constituants=self.constituants,
            quantites=tuple(float(q) for q in n),
            equilibre=bool(equilibre),
            disparus=disparus,
            Q=Q,
            K0=K0,
        )


def _quantite(nom: str, n: float) -> float:
    if n < 0:
        raise ValueError(f"quantité négative pour {nom} : {n}")
    return float(n)


def _grille(xi: ArrayLike) -> tuple[Tableau, tuple[int, ...]]:
    x = np.asarray(xi, dtype=float)
    return x.reshape(-1), x.shape

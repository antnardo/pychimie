"""Déplacements d'équilibre : perturber un équilibre, balayer un paramètre.

Le raisonnement de cours tient en trois temps : un système est à l'équilibre (Q = K°) ;
on change une condition — T, P, un ajout de constituant actif ou inactif — sans laisser
la composition changer ; on compare le nouveau Q au nouveau K°. :func:`perturber` fait
exactement cela et rend les trois temps : l'équilibre de départ, Q et K° juste après la
perturbation (d'où le sens d'évolution), et le nouvel état final, dont l'avancement
*compté depuis l'état perturbé* est le déplacement δξ.

Les lois de modération en découlent et se vérifient numériquement :

- une hausse de T déplace dans le sens endothermique (ΔrH° > 0 : sens direct) ;
- une hausse de P à T fixée, dans le sens qui diminue la quantité de gaz
  (signe de −Σν_gaz) ;
- un gaz inactif ajouté à T et P fixées agit comme une baisse de pression ; à T et
  V fixés, il ne déplace rien (Q ne dépend que des nᵢRT/V) ;
- un réactif ajouté à T et P fixées ne déplace pas toujours dans le sens direct :
  N₂ ajouté à un mélange où x(N₂) > 1/2 déplace la synthèse de NH₃ dans le sens
  inverse.

:func:`balayage` refait un même calcul pour une suite de valeurs d'un paramètre —
température, pression, quantité ajoutée, proportions initiales — et rend ξ_f, le taux
d'avancement et les fractions molaires finales en tableaux, prêts à tracer. Il se
donne comme :func:`perturber`, ``balayage(s, T=…)``, ``balayage(s, P=…)``,
``balayage(s, ajout={"Ar(g)": …})``, ou par une fabrique quelconque de systèmes.
"""

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from pychimie.equation import Constituant
from pychimie.systeme import EtatFinal, Systeme

__all__ = ["Balayage", "Perturbation", "balayage", "perturber"]


@dataclass(frozen=True, slots=True)
class Perturbation:
    """Un équilibre, la perturbation, et l'état final qui la suit."""

    avant: EtatFinal
    systeme: Systeme
    Q: float
    K0: float
    apres: EtatFinal

    @property
    def sens(self) -> int:
        """+1 si Q < K° juste après la perturbation (sens direct), −1 si Q > K°, 0 sinon."""
        return (self.Q < self.K0) - (self.Q > self.K0)

    @property
    def deplacement(self) -> float:
        """δξ, compté depuis l'état perturbé ; son signe est ``sens`` (sauf blocage)."""
        return self.apres.xi


def _ajouter(
    systeme: Systeme, initial: dict[str, float], inertes: dict[str, float], nom: str, dn: float
) -> None:
    """Ajoute ``dn`` mol de ``nom`` — un constituant de la réaction, ou un inactif écrit
    avec sa phase (``"Ar(g)"``) — à une composition."""
    equation = systeme.reaction.equation
    try:
        initial[equation.noms[equation.indice(nom)]] += dn
    except KeyError:
        if "(" not in nom:  # ni un constituant de la réaction, ni un inactif avec sa phase
            raise
        cle = Constituant.lire(nom).nom
        inertes[cle] = inertes.get(cle, 0.0) + dn


def perturber(
    systeme: Systeme, *, ajout: Mapping[str, float] | None = None, **conditions: Any
) -> Perturbation:
    """Perturbe l'état final de ``systeme`` : nouvelles ``conditions`` (T, P, V, V_gaz),
    ``ajout`` de constituants (``{"N2": 0.1}``, ``{"Ar(g)": 1}`` pour un inactif)."""
    avant = systeme.etat_final()
    r = len(systeme.initial)
    initial = dict(zip(systeme.reaction.equation.noms, avant.quantites[:r], strict=True))
    inertes = {c.nom: n for c, n in zip(avant.constituants[r:], avant.quantites[r:], strict=True)}
    for nom, dn in (ajout or {}).items():
        _ajouter(systeme, initial, inertes, nom, dn)
    perturbe = systeme.avec(initial=initial, inertes=inertes, **conditions)
    return Perturbation(
        avant=avant,
        systeme=perturbe,
        Q=float(perturbe.Q(0.0)),
        K0=float(perturbe.K0()),
        apres=perturbe.etat_final(),
    )


@dataclass(frozen=True, slots=True)
class Balayage:
    """Les états finaux d'une suite de systèmes, un par valeur du paramètre."""

    valeurs: NDArray[np.float64]
    etats: tuple[EtatFinal, ...]
    parametre: str = ""
    """``"T"``, ``"P"``, ``"V"``, ``"V_gaz"``, le nom du constituant ajouté, ou ``""``
    pour une fabrique quelconque."""

    @property
    def xi(self) -> NDArray[np.float64]:
        return np.array([e.xi for e in self.etats])

    @property
    def xi_max(self) -> NDArray[np.float64]:
        return np.array([e.xi_max for e in self.etats])

    @property
    def taux_avancement(self) -> NDArray[np.float64]:
        return np.array([e.taux_avancement for e in self.etats])

    @property
    def equilibre(self) -> NDArray[np.bool_]:
        """Faux là où l'état final est une rupture d'équilibre."""
        return np.array([e.equilibre for e in self.etats])

    def n(self, nom: str) -> NDArray[np.float64]:
        return np.array([e.n(nom) for e in self.etats])

    def fraction(self, nom: str) -> NDArray[np.float64]:
        return np.array([e.fraction(nom) for e in self.etats])


def balayage(
    systeme: Systeme | Callable[[float], Systeme],
    valeurs: ArrayLike | None = None,
    *,
    ajout: Mapping[str, ArrayLike] | None = None,
    **conditions: ArrayLike,
) -> Balayage:
    """L'état final le long d'un paramètre : les mêmes trois que :func:`perturber`.

    - ``balayage(s, T=np.linspace(400, 900, 50))``, ``balayage(s, P=[1, 10, 100])``,
      de même ``V`` et ``V_gaz`` : une condition remplacée ;
    - ``balayage(s, ajout={"Ar(g)": [0, 1, 5]})``, ``ajout={"N2": …}`` : une quantité
      ajoutée à la composition **initiale** (:func:`perturber` l'ajoute, lui, à
      l'équilibre) ;
    - ``balayage(lambda x: s.avec(initial={"N2": x, "H2": 1 - x}), x)`` : une fabrique
      de systèmes, pour tout autre paramètre (proportions initiales…).

    Un seul paramètre à la fois.
    """
    if callable(systeme):
        if valeurs is None or ajout or conditions:
            raise TypeError(
                "avec une fabrique, donner les valeurs seules : balayage(fabrique, valeurs)"
            )
        return _balayer(systeme, valeurs, "")
    parametres = {**conditions, **(ajout or {})}
    if valeurs is not None or len(parametres) != 1:
        raise TypeError(
            "un seul paramètre : balayage(s, T=…), balayage(s, P=…), balayage(s, V=…), "
            "balayage(s, V_gaz=…) ou balayage(s, ajout={nom: …})"
        )
    ((nom, grille),) = parametres.items()
    if nom in conditions:
        return _balayer(lambda v: systeme.avec(**{nom: v}), grille, nom)
    r = len(systeme.initial)
    equation = systeme.reaction.equation

    def avec_ajout(dn: float) -> Systeme:
        initial = dict(zip(equation.noms, systeme.initial[:r], strict=True))
        inertes = {c.nom: n for c, n in systeme.inertes}
        _ajouter(systeme, initial, inertes, nom, dn)
        return systeme.avec(initial=initial, inertes=inertes)

    return _balayer(avec_ajout, grille, nom)


def _balayer(fabrique: Callable[[float], Systeme], valeurs: ArrayLike, parametre: str) -> Balayage:
    v = np.asarray(valeurs, dtype=float).reshape(-1)
    return Balayage(v, tuple(fabrique(float(x)).etat_final() for x in v), parametre)

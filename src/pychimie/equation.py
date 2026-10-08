"""L'équation d'une réaction, telle qu'on l'écrit : ``"N2(g) + 3 H2(g) = 2 NH3(g)"``.

Chaque constituant porte sa phase entre parenthèses, à la fin : ``g``, ``l``, ``s`` ou
``aq``. La phase dit comment se calcule l'activité, et une seule règle couvre tous les
cas : **l'activité d'un constituant d'une phase condensée est sa fraction molaire dans
cette phase**.

- ``s`` ou ``l`` sans suffixe : le constituant est **seul dans sa phase** (corps pur,
  solvant) ; sa fraction molaire vaut 1, son activité aussi, même quand il disparaît —
  d'où les ruptures d'équilibre.
- ``s`` ou ``l`` suivi d'une étiquette (``l1``, ``s2``, ``lorg``…) : les constituants
  de même étiquette forment **un mélange idéal**, comme l'acide, l'alcool, l'ester et
  l'eau d'une estérification écrite ``(l1)``.
- ``g`` : tous les gaz forment un seul mélange idéal ; activité x·P/P°.
- ``aq`` (suffixe accepté, ``aq1``…) : soluté, activité c/c° ; tous les solutés
  partagent le volume de la solution.

Le séparateur des constituants est `` + `` entouré d'espaces, pour que les charges
(``Fe2+(aq)``, ``H+(aq)``) restent dans les formules. Les deux membres sont séparés par
``=``, ``⇌``, ``->`` ou ``→``. Un coefficient peut être entier, décimal ou fractionnaire
(``1/2 O2(g)``).
"""

import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from fractions import Fraction

__all__ = ["Constituant", "Equation", "indice"]

_SEPARATEURS = re.compile(r"\s*(?:=|⇌|->|→)\s*")
_PLUS = re.compile(r"\s+\+\s+")
_TERME = re.compile(
    r"^(?P<coef>\d+/\d+|\d+(?:[.,]\d+)?)?\s*(?P<formule>\S.*?)\((?P<phase>g|aq\w*|s\w*|l\w*)\)$"
)


@dataclass(frozen=True, slots=True)
class Constituant:
    """Un constituant et son nombre stœchiométrique : négatif pour un réactif, nul pour
    un constituant inactif (gaz inerte, solvant d'un mélange)."""

    formule: str
    phase: str
    nu: float = 0.0

    @classmethod
    def lire(cls, texte: str) -> "Constituant":
        """Lit ``"Ar(g)"`` : un constituant inactif, sans coefficient."""
        trouve = _TERME.match(texte.strip())
        if trouve is None or trouve["coef"]:
            raise ValueError(
                f"constituant illisible : {texte!r} (attendu « Ar(g) », avec sa phase)"
            )
        return cls(trouve["formule"], trouve["phase"])

    @property
    def nom(self) -> str:
        """``"NH3(g)"`` : la formule et la phase, ce qui distingue H2O(l) de H2O(g)."""
        return f"{self.formule}({self.phase})"

    @property
    def etat(self) -> str:
        """``"g"``, ``"aq"``, ``"s"`` ou ``"l"``, sans l'étiquette de la phase."""
        return "aq" if self.phase.startswith("aq") else self.phase[0]

    @property
    def melange(self) -> str | None:
        """La phase dont le constituant partage la fraction molaire, ou ``None`` (soluté).

        Un ``s`` ou un ``l`` sans étiquette est seul dans sa phase : la clé est alors
        son propre nom, et aucun autre constituant ne la partage.
        """
        if self.etat == "aq":
            return None
        if self.phase in ("s", "l"):
            return self.nom
        return self.phase


def indice(constituants: Sequence[Constituant], nom: str) -> int:
    """La position du constituant ``nom`` : ``"NH3(g)"``, ou ``"NH3"`` s'il est unique.

    Une formule seule désigne le constituant quand une seule phase la porte ; avec H2O(l)
    et H2O(g) dans la même équation, il faut écrire la phase.
    """
    nom = nom.strip()
    exacts = [i for i, c in enumerate(constituants) if c.nom == nom]
    if exacts:
        return exacts[0]
    formules = [i for i, c in enumerate(constituants) if c.formule == nom]
    if len(formules) == 1:
        return formules[0]
    disponibles = ", ".join(c.nom for c in constituants)
    if formules:
        raise KeyError(f"{nom} est ambigu : préciser la phase parmi {disponibles}")
    raise KeyError(f"{nom} n'est pas un constituant ({disponibles})")


@dataclass(frozen=True, slots=True)
class Equation:
    """Les constituants d'une réaction, réactifs puis produits, avec leurs ν."""

    constituants: tuple[Constituant, ...]

    def __post_init__(self) -> None:
        noms = [c.nom for c in self.constituants]
        doublons = {n for n in noms if noms.count(n) > 1}
        if doublons:
            raise ValueError(f"constituant écrit deux fois : {', '.join(sorted(doublons))}")
        if not any(c.nu < 0 for c in self.constituants):
            raise ValueError("l'équation n'a pas de réactif")
        if not any(c.nu > 0 for c in self.constituants):
            raise ValueError("l'équation n'a pas de produit")
        if any(c.nu == 0 for c in self.constituants):
            raise ValueError("un coefficient stœchiométrique est nul")

    @classmethod
    def lire(cls, texte: str) -> "Equation":
        """Lit ``"C(s) + CO2(g) = 2 CO(g)"``."""
        membres = _SEPARATEURS.split(texte.strip())
        if len(membres) != 2:
            raise ValueError(f"une équation a deux membres séparés par « = » : {texte!r}")
        reactifs, produits = membres
        return cls((*_termes(reactifs, -1), *_termes(produits, +1)))

    @property
    def noms(self) -> tuple[str, ...]:
        return tuple(c.nom for c in self.constituants)

    @property
    def nu(self) -> tuple[float, ...]:
        return tuple(c.nu for c in self.constituants)

    @property
    def delta_nu_gaz(self) -> float:
        """Σν des gaz : son signe dit dans quel sens une hausse de P déplace l'équilibre."""
        return sum(c.nu for c in self.constituants if c.etat == "g")

    def indice(self, nom: str) -> int:
        return indice(self.constituants, nom)

    def __str__(self) -> str:
        def membre(termes: Iterable[Constituant]) -> str:
            return " + ".join(f"{_coef(abs(c.nu))}{c.nom}" for c in termes)

        reactifs = membre(c for c in self.constituants if c.nu < 0)
        produits = membre(c for c in self.constituants if c.nu > 0)
        return f"{reactifs} = {produits}"


def _termes(membre: str, signe: int) -> list[Constituant]:
    termes = []
    for terme in _PLUS.split(membre.strip()):
        trouve = _TERME.match(terme.strip())
        if trouve is None:
            raise ValueError(
                f"constituant illisible : {terme!r} (attendu « 2 NH3(g) », phase g, l, s ou aq)"
            )
        coef = trouve["coef"]
        nu = float(Fraction(coef.replace(",", "."))) if coef else 1.0
        termes.append(Constituant(trouve["formule"], trouve["phase"], signe * nu))
    return termes


def _coef(nu: float) -> str:
    if nu == 1:
        return ""
    fraction = Fraction(nu).limit_denominator(12)
    return f"{fraction} " if float(fraction) == nu else f"{nu:g} "

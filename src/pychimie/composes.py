"""Un composé : une formule, et l'état physique qu'on lui donne dans une équation.

Les équations de réaction écrivent « H2O(l) », « CO2(g) », « Cu^2+(aq) » : la phase
fait partie de l'identité de l'espèce (l'eau liquide et la vapeur n'ont ni la même
enthalpie de formation, ni la même activité). :class:`Compose` sépare cette mention
de la formule, que :class:`~pychimie.formules.Formule` analyse.
"""

import re
from dataclasses import dataclass

from pychimie.elements import TableauPeriodique
from pychimie.formules import Formule

__all__ = ["PHASES", "Compose"]

PHASES = {
    "s": "solide",
    "l": "liquide",
    "g": "gaz",
    "aq": "en solution aqueuse",
}

_PHASE = re.compile(r"^(?P<formule>.+?)\s*\((?P<phase>s|l|g|aq)\)$")


@dataclass(frozen=True, slots=True)
class Compose:
    """Une espèce chimique : formule, phase éventuelle et nom d'usage."""

    formule: Formule
    phase: str | None = None
    nom: str = ""

    def __post_init__(self) -> None:
        if self.phase is not None and self.phase not in PHASES:
            raise ValueError(f"phase inconnue : {self.phase!r} (attendu : {', '.join(PHASES)})")

    @classmethod
    def lire(cls, texte: str, nom: str = "", tableau: TableauPeriodique | None = None) -> "Compose":
        """« H2O(l) » → eau liquide ; sans parenthèse finale, la phase reste inconnue.

        Seules les quatre mentions (s), (l), (g) et (aq) sont des phases : « Ca(OH)2 »
        garde ses parenthèses, qui sont celles de la formule.
        """
        texte = texte.strip()
        m = _PHASE.match(texte)
        formule, phase = (m["formule"], m["phase"]) if m else (texte, None)
        return cls(Formule.lire(formule, tableau), phase, nom)

    def __str__(self) -> str:
        return self.formule.texte + (f"({self.phase})" if self.phase else "")

    @property
    def masse_molaire(self) -> float:
        """Masse molaire en g/mol."""
        return self.formule.masse_molaire()

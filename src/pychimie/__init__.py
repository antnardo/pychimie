"""Chimie pour la classe : éléments, formules, solutions.

>>> from pychimie import element, masse_molaire
>>> element("Fe").masse_molaire
55.845
>>> round(masse_molaire("CuSO4·5H2O"), 2)
249.68
"""

from pychimie.composes import PHASES, Compose
from pychimie.elements import (
    CARACTERES,
    FAMILLES,
    Element,
    TableauPeriodique,
    element,
    tableau_periodique,
)
from pychimie.formules import AnalyseurFormule, Formule, FormuleInvalide, masse_molaire
from pychimie.solutions import RHO_EAU, Solution, masse_a_peser, volume_a_prelever

__all__ = [
    "CARACTERES",
    "FAMILLES",
    "PHASES",
    "RHO_EAU",
    "AnalyseurFormule",
    "Compose",
    "Element",
    "Formule",
    "FormuleInvalide",
    "Solution",
    "TableauPeriodique",
    "element",
    "masse_a_peser",
    "masse_molaire",
    "tableau_periodique",
    "volume_a_prelever",
]

"""Chimie pour la classe : éléments, formules, solutions, tableau périodique en TikZ.

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
from pychimie.periodikz import ARTICLES_WIKIPEDIA, TableauTikz, coordonnees, url_wikipedia
from pychimie.pgfperiodictable import (
    BLANC,
    GRIS,
    Couleur,
    rgb,
    schema_categories,
    schema_degrade,
)
from pychimie.solutions import RHO_EAU, Solution, masse_a_peser, volume_a_prelever

__all__ = [
    "ARTICLES_WIKIPEDIA",
    "BLANC",
    "CARACTERES",
    "FAMILLES",
    "GRIS",
    "PHASES",
    "RHO_EAU",
    "AnalyseurFormule",
    "Compose",
    "Couleur",
    "Element",
    "Formule",
    "FormuleInvalide",
    "Solution",
    "TableauPeriodique",
    "TableauTikz",
    "coordonnees",
    "element",
    "masse_a_peser",
    "masse_molaire",
    "rgb",
    "schema_categories",
    "schema_degrade",
    "tableau_periodique",
    "url_wikipedia",
    "volume_a_prelever",
]

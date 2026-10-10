"""Chimie pour la classe : équilibres, éléments, formules, préparations, tableau périodique.

Thermochimie des équilibres : avancement, Q(ξ), K°(T), G(ξ), rupture, déplacements.

>>> from pychimie import Reaction, Systeme
>>> nh3 = Reaction.depuis_tables("N2(g) + 3 H2(g) = 2 NH3(g)")
>>> s = Systeme(nh3, {"N2": 1, "H2": 3}, T=700, P=200)
>>> round(s.xi_eq(), 4)
0.578

Éléments et formules :

>>> from pychimie import element, masse_molaire
>>> element("Fe").masse_molaire
55.845
>>> round(masse_molaire("CuSO4·5H2O"), 2)
249.68
"""

from pychimie.deplacements import Balayage, Perturbation, balayage, perturber
from pychimie.donnees import DonneesThermo, TableThermo, table_thermo
from pychimie.elements import FAMILLES, Element, TableauPeriodique, element, tableau_periodique
from pychimie.equation import Constituant, Equation
from pychimie.formules import AnalyseurFormule, Formule, FormuleInvalide, masse_molaire
from pychimie.periodikz import TableauTikz, coordonnees, url_wikipedia
from pychimie.pgfperiodictable import (
    BLANC,
    GRIS,
    Couleur,
    rgb,
    schema_categories,
    schema_degrade,
)
from pychimie.preparation import Preparation, masse_a_peser, volume_a_prelever
from pychimie.reaction import Reaction
from pychimie.systeme import EtatFinal, RuptureEquilibre, Systeme

__all__ = [
    "BLANC",
    "FAMILLES",
    "GRIS",
    "AnalyseurFormule",
    "Balayage",
    "Constituant",
    "Couleur",
    "DonneesThermo",
    "Element",
    "Equation",
    "EtatFinal",
    "Formule",
    "FormuleInvalide",
    "Perturbation",
    "Preparation",
    "Reaction",
    "RuptureEquilibre",
    "Systeme",
    "TableThermo",
    "TableauPeriodique",
    "TableauTikz",
    "balayage",
    "coordonnees",
    "element",
    "masse_a_peser",
    "masse_molaire",
    "perturber",
    "rgb",
    "schema_categories",
    "schema_degrade",
    "table_thermo",
    "tableau_periodique",
    "url_wikipedia",
    "volume_a_prelever",
]

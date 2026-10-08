"""Thermochimie des équilibres : avancement, Q(ξ), K°(T), G(ξ), rupture, déplacements.

>>> from pychimie import Reaction, Systeme
>>> nh3 = Reaction.depuis_tables("N2(g) + 3 H2(g) = 2 NH3(g)")
>>> s = Systeme(nh3, {"N2": 1, "H2": 3}, T=700, P=200)
>>> round(s.xi_eq(), 4)
0.578
"""

from pychimie.deplacements import Balayage, Perturbation, balayage, perturber
from pychimie.donnees import DonneesThermo, TableThermo, table_thermo
from pychimie.equation import Constituant, Equation
from pychimie.reaction import Reaction
from pychimie.systeme import EtatFinal, RuptureEquilibre, Systeme

__all__ = [
    "Balayage",
    "Constituant",
    "DonneesThermo",
    "Equation",
    "EtatFinal",
    "Perturbation",
    "Reaction",
    "RuptureEquilibre",
    "Systeme",
    "TableThermo",
    "balayage",
    "perturber",
    "table_thermo",
]

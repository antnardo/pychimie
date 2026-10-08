r"""Schémas de couleurs pour le paquet LaTeX pgf-PeriodicTable.

pgf-PeriodicTable (CTAN, LPPL 1.3) dessine un tableau périodique complet, en
français avec l'option ``fr``, mais ne colore les cases que selon ses propres schémas.
Sa commande ``\pgfPTnewColorScheme{nom}{r/g/b,…}`` en accepte un nouveau : une
couleur RVB (composantes entre 0 et 1) par élément, de Z = 1 à 118. Ce module
l'écrit à partir de n'importe quel attribut de :class:`~pychimie.elements.Element`,
qualitatif ou numérique :

.. code-block:: latex

    \usepackage[fr]{pgf-PeriodicTable}
    \input{schema.tex}  % \pgfPTnewColorScheme{chi}{…}
    \pgfPT[back color scheme=chi]

C'est le complément de :mod:`pychimie.periodikz` quand on préfère la mise en page
soignée et les données du paquet à un gabarit maison.
"""

import re
from collections.abc import Mapping
from typing import Any

from pychimie.elements import TableauPeriodique, tableau_periodique

__all__ = [
    "BLANC",
    "GRIS",
    "Couleur",
    "rgb",
    "schema_categories",
    "schema_degrade",
]

Couleur = tuple[float, float, float]

BLANC: Couleur = (1.0, 1.0, 1.0)
GRIS: Couleur = (0.8, 0.8, 0.8)

_HEX = re.compile(r"#?([0-9a-fA-F]{2})([0-9a-fA-F]{2})([0-9a-fA-F]{2})")


def rgb(couleur: Couleur | str) -> Couleur:
    """Une couleur en triplet de composantes entre 0 et 1, depuis « #rrggbb » ou un triplet."""
    if isinstance(couleur, str):
        m = _HEX.fullmatch(couleur.strip())
        if m is None:
            raise ValueError(f"couleur hexadécimale attendue (#rrggbb) : {couleur!r}")
        r, v, b = (int(c, 16) / 255 for c in m.groups())
        return r, v, b
    if len(couleur) != 3 or not all(0 <= c <= 1 for c in couleur):
        raise ValueError(f"trois composantes entre 0 et 1 attendues : {couleur!r}")
    r, v, b = (float(c) for c in couleur)
    return r, v, b


def _commande(nom: str, couleurs: list[Couleur]) -> str:
    if not re.fullmatch(r"[A-Za-z]+", nom):
        # Le nom finit dans un \csname et dans la clé « back color scheme » : des
        # lettres seules évitent toute surprise de développement.
        raise ValueError(f"le nom d'un schéma ne contient que des lettres : {nom!r}")
    triplets = ",".join("/".join(format(round(c, 3), "g") for c in rvb) for rvb in couleurs)
    return rf"\pgfPTnewColorScheme{{{nom}}}{{{triplets}}}" + "\n"


def schema_categories(
    nom: str,
    attribut: str,
    couleurs: Mapping[Any, Couleur | str],
    *,
    defaut: Couleur | str = BLANC,
    tableau: TableauPeriodique | None = None,
) -> str:
    r"""Un schéma qui donne à chaque valeur de l'attribut sa couleur.

    >>> schema_categories("etats", "etat", {"g": "#ff0000"})[:40]
    '\\pgfPTnewColorScheme{etats}{1/0/0,1/0/0,'
    """
    t = tableau or tableau_periodique()
    table = {valeur: rgb(c) for valeur, c in couleurs.items()}
    fond = rgb(defaut)
    return _commande(nom, [table.get(getattr(e, attribut), fond) for e in t])


def schema_degrade(
    nom: str,
    attribut: str,
    debut: Couleur | str,
    fin: Couleur | str,
    *,
    manquant: Couleur | str = GRIS,
    bornes: tuple[float, float] | None = None,
    tableau: TableauPeriodique | None = None,
) -> str:
    """Un schéma en dégradé linéaire de ``debut`` à ``fin`` selon une propriété numérique.

    La plus petite valeur reçoit ``debut``, la plus grande ``fin`` (ou les ``bornes``
    données, voir :meth:`~pychimie.elements.TableauPeriodique.normaliser`) ; un
    élément sans valeur reçoit ``manquant``.
    """
    t = tableau or tableau_periodique()
    a, b, vide = rgb(debut), rgb(fin), rgb(manquant)
    x = t.normaliser(attribut, bornes)

    def melange(s: float | None) -> Couleur:
        if s is None:
            return vide
        r, v, bl = (ca + s * (cb - ca) for ca, cb in zip(a, b, strict=True))
        return r, v, bl

    return _commande(nom, [melange(x[e.Z]) for e in t])

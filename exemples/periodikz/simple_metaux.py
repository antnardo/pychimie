"""Tableau périodique des familles, avec masses molaires et liens vers Wikipédia.

Produit ``simple_metaux.tex``, que le gabarit ``Tableau_simple_metaux.tex`` inclut :

    python simple_metaux.py && lualatex Tableau_simple_metaux.tex

Les cases sont colorées selon la famille, les symboles selon l'état à 25 °C ; chaque
symbole est un lien vers l'article de Wikipédia en français.
"""

from pathlib import Path

from pychimie import FAMILLES, TableauTikz

ICI = Path(__file__).resolve().parent

tableau = (
    TableauTikz(liens=True)
    .cases("famille")
    .symboles("etat")
    .infos("masse_molaire", "M info", r"\num{{{}}}")
    .infos("nom", "nom")
    .infos("Z", "Z")
)
tableau.legende("famille", ligne=1, horizontale=True, valeurs=FAMILLES)
print(tableau.ecrire(ICI / "simple_metaux.tex"))

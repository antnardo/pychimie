"""Tableau périodique ombré selon l'électronégativité de Pauling.

Produit ``ombrage_chi.tex``, que le gabarit ``Tableau_ombrage_chi.tex`` inclut :

    python ombrage_chi.py && lualatex Tableau_ombrage_chi.tex

:meth:`TableauTikz.ombrage` ramène la propriété à un entier de 0 à 100 et l'écrit en
argument d'un style du gabarit, ``\\dessinecase[chi=73]{1-1}`` : c'est le gabarit qui
en fait une couleur (``chi/.style={fill=red!#1!blue}``). Les éléments sans valeur
reçoivent le style ``sans chi``. Les bornes fixées, (0, 4), rendent l'échelle lisible
dans la légende ; sans elles, 0 et 100 iraient aux valeurs extrêmes de la table.
"""

from pathlib import Path

from pychimie import TableauTikz

ICI = Path(__file__).resolve().parent

tableau = (
    TableauTikz(liens=False)
    .ombrage("electronegativite", "chi", bornes=(0, 4), defaut="sans chi")
    .symboles()
    .infos("electronegativite", "chi info", r"\num{{{}}}")
    .infos("nom", "nom")
    .infos("Z", "Z")
)
print(tableau.ecrire(ICI / "ombrage_chi.tex"))

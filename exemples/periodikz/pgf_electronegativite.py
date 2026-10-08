"""Schémas de couleurs pour pgf-PeriodicTable : électronégativité et familles.

Produit ``schemas_pgfpt.tex``, que ``Tableau_pgf_electronegativite.tex`` inclut :

    python pgf_electronegativite.py && lualatex Tableau_pgf_electronegativite.tex

Le paquet LaTeX pgf-PeriodicTable dessine le tableau et ses données ; pychimie ne
fournit que la couleur de fond de chaque case.
"""

from pathlib import Path

from pychimie import schema_categories, schema_degrade

ICI = Path(__file__).resolve().parent

COULEURS_FAMILLES = {
    "métal alcalin": "#f08080",
    "métal alcalino-terreux": "#ccb833",
    "lanthanide": "#f5d0e0",
    "actinide": "#f5b3b3",
    "métal de transition": "#f5cccc",
    "métal pauvre": "#a6a6a6",
    "métalloïde": "#99b380",
    "non-métal": "#8fd9a8",
    "halogène": "#f5e066",
    "gaz noble": "#b3d9f0",
}

schemas = schema_degrade(
    "chi", "electronegativite", "#d62728", "#1f77b4", manquant="#cccccc"
) + schema_categories("familles", "famille", COULEURS_FAMILLES)
chemin = ICI / "schemas_pgfpt.tex"
chemin.write_text(schemas, encoding="utf-8")
print(chemin)

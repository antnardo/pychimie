# pychimie

Thermochimie (avancement, équilibre, rupture), données atomiques et tableau périodique en TikZ

Code d'origine importé tel quel le 8 octobre 2026 ; voir `../AUDIT.md`.

## Éléments, formules et tableau périodique

### Les 118 éléments

`element()` rend un élément par son numéro atomique, son symbole ou son nom français
(sans tenir compte de la casse ni des accents) ; `tableau_periodique()` rend la table
entière, que l'on filtre ou dont on extrait une colonne.

```python
from pychimie import element, tableau_periodique

fer = element("fer")             # ou element(26), element("Fe")
fer.masse_molaire                # 55.845 (g/mol, CIAAW)
fer.theta_fusion                 # 1538.0 (°C) ; fer.T_fusion en K
fer.famille, fer.configuration   # 'métal de transition', '[Ar] 3d6 4s2'

t = tableau_periodique()
[e.symbole for e in t.filtrer(famille="halogène")]
t.valeurs("electronegativite")   # {'H': 2.2, 'He': None, …}
```

Une grandeur a une seule unité (g/mol, pm, eV, K, kg/m³) et vaut `None` quand elle
est inconnue. La table est figée dans `src/pychimie/donnees/elements.csv` ; ses
sources et leurs licences — masses de la CIAAW, propriétés du paquet mendeleev
(MIT), aucune donnée de Wikipédia — sont détaillées dans
`src/pychimie/donnees/SOURCES.md`.

### Formules brutes et solutions

```python
from pychimie import Formule, Solution, masse_a_peser, masse_molaire

masse_molaire("CuSO4·5H2O")              # 249.68 g/mol
Formule.lire("K4[Fe(CN)6]").hill         # 'C6FeK4N6'
Formule.lire("SO₄²⁻").charge             # -2
masse_a_peser("CuSO4·5H2O", 0.10, 0.1)   # 2.50 g pour 100 mL à 0,10 mol/L
Solution.commerciale("H2SO4", 0.96, 1.84).concentration   # 18.0 mol/L
```

Les hydrates (`·`, `.` ou `*`), les parenthèses et crochets imbriqués, les indices et
exposants Unicode sont lus ; une charge s'écrit `^2-`, `²⁻` ou, pour ±1, d'un signe
seul (`NH4+`). Les équilibres en solution (pH, titrages) relèvent du paquet
`titrages`.

### Tableau périodique en TikZ

`TableauTikz` produit, élément par élément, les commandes qu'un gabarit TikZ met en
forme : cases stylées selon un attribut (`famille=halogène`), ombrées selon une
propriété numérique (`chi=73`), symboles en hyperliens vers Wikipédia, informations
et légende.

```python
from pychimie import FAMILLES, TableauTikz

t = TableauTikz(liens=True).cases("famille").symboles("etat")
t.infos("masse_molaire", "M info", r"\num{{{}}}").infos("nom", "nom").infos("Z", "Z")
t.legende("famille", ligne=1, horizontale=True, valeurs=FAMILLES)
t.ecrire("simple_metaux.tex")
```

Le gabarit `exemples/periodikz/Tableau_simple_metaux.tex` définit les coordonnées
(`<groupe>-<période>`, `<rang>-L` et `<rang>-A` pour le bloc f), les commandes
`\dessinecase`, `\legend`, `\legendhoriz` et les styles. Pour le paquet LaTeX
pgf-PeriodicTable, `schema_categories` et `schema_degrade` écrivent un
`\pgfPTnewColorScheme` (exemple : `exemples/periodikz/pgf_electronegativite.py`).

### Régénérer la table

```bash
uv run python scripts/telecharger_ciaaw.py
uv run --group donnees python scripts/generer_elements.py
```

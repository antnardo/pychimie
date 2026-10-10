# pychimie

Thermochimie (avancement, équilibre, rupture), données atomiques et tableau périodique en TikZ

Les équilibres en solution (pH, titrages pH-métriques, conductimétriques et
potentiométriques) relèvent du paquet `titrages`. L'historique de la reprise de
l'ancien module `chimie` est dans `../AUDIT.md` (8 octobre 2026).

## Équilibres chimiques

Le module de thermochimie suit le cours de CPGE : une réaction, son K°(T), et un système
qui évolue par cette réaction, de l'état initial à l'état final — équilibre ou rupture
d'équilibre. Il remplace l'ancien `chimie.thermochimie.Reaction`.

```python
from pychimie import Reaction, Systeme, perturber

nh3 = Reaction("N2(g) + 3 H2(g) = 2 NH3(g)", DrH0=-92.2e3, DrS0=-198.8)  # J/mol, J/K/mol
s = Systeme(nh3, {"N2": 1, "H2": 3}, T=700, P=200)  # mol, K, bar
s.xi_eq()  # 0.5774 mol
xi = s.avancements()  # ξ entre ξ_min et ξ_max, bornes exclues
s.Q(xi), s.DrG(xi), s.G(xi)
perturber(s, P=300).sens  # +1 : une hausse de P déplace vers NH₃
```

### Réaction et K°(T)

- `Reaction(équation, DrH0=…, DrS0=…)` : les valeurs d'un énoncé. Elles sont
  obligatoires — l'ancien code prenait 1 J/mol et 1 J/K/mol quand on les oubliait.
- `Reaction.depuis_tables(équation)` : ΔrH° et ΔrS° par la loi de Hess, à partir de la
  table de `src/pychimie/donnees/` ; `kirchhoff=True` garde ΔrC°p.
- `Reaction.depuis_K0(équation, K0, T=…, DrH0=…)` : K° à une température, et la loi de
  van 't Hoff si ΔrH° est donné ; sans ΔrH°, K° ne dépend pas de T.
- `r.K0(T)`, `r.ln_K0(T)`, `r.DrG0(T)`, `r.DrH0(T)`, `r.DrS0(T)` acceptent des tableaux.
  L'approximation d'Ellingham est le cas ΔrC°p = 0 (`r.ellingham()`).

L'équation s'écrit avec la phase de chaque constituant : `g`, `l`, `s` ou `aq`. Un `s`
ou un `l` sans étiquette est un corps pur (activité 1) ; `l1`, `s2`… désignent un
mélange idéal (activité égale à la fraction molaire). Le solvant est un constituant
comme un autre : donner sa quantité (55,5 mol pour un litre d'eau).

### Système, état final, rupture

Unités : quantités en mol, `T` en K, `P` en bar (P° = 1 bar), volumes en L
(c° = 1 mol/L). Une pression au-delà de 10⁴ bar ou une température au-dessous de
100 K est refusée, avec la conversion dans le message : `P=200e5` (des pascals) et
`T=25` (des degrés Celsius) ne passent pas en silence.

`Systeme(reaction, initial, T=…, P=…, V=…, V_gaz=…, inertes=…)` : `P` est la pression
imposée (bar), `V` le volume de la solution (L) quand il y a des solutés, `V_gaz` le
volume du gaz s'il est fixé à la place de P, `inertes` les constituants inactifs
(`{"Ar(g)": 1}`).

- `K0()`, `DrG0()` : K° et ΔrG° à la température du système ; `K0(T)`, `DrG0(T)` à
  une autre, comme `Reaction.K0(T)`.
- `xi_min`, `xi_max`, `avancements()`, `n(nom, ξ)`, `activite(nom, ξ)`, `Q(ξ)`,
  `DrG(ξ)`, `G(ξ)` (à une constante près : G(0) = 0), `minimum_G()`.
- `etat_final()` : avancement, quantités, fractions molaires, Q et K°, et, en cas de
  rupture d'équilibre, les constituants disparus. `xi_eq()` lève `RuptureEquilibre`
  quand l'état final n'est pas un équilibre.
- `avec(T=…, P=…, initial=…)` : le même système dans d'autres conditions.

### Déplacements d'équilibre

- `perturber(s, T=…)`, `perturber(s, P=…)`, `perturber(s, ajout={"Ar(g)": 1})` :
  l'équilibre de départ, Q et K° juste après la perturbation, le sens d'évolution et le
  déplacement δξ.
- `balayage(lambda T: s.avec(T=T), valeurs)` : ξ_f, taux d'avancement, quantités et
  fractions finales le long d'un paramètre (T, P, quantité d'inerte, proportions
  initiales).

Les scripts de `exemples/equilibres/` refont ceux du cours (`Optimisation_NH3_*.py`,
`Optimisation_Boudouard.py`) et y ajoutent K°(T) avec Kirchhoff, les déplacements
d'équilibre et la rupture d'équilibre du carbonate de calcium.

### Données

`src/pychimie/donnees/thermo298.csv` donne ΔfH°, S°m et C°p à 298,15 K pour 127
espèces : gaz simples et composés de l'azote, du soufre et du phosphore (NH₃, NO,
NO₂, N₂O₄, SO₂, SO₃, H₂S, H₂SO₄, HNO₃, PCl₃, PCl₅, NH₄Cl, urée…), halogènes et
hydracides, le carbone et ses composés jusqu'à C₂ (CO, CO₂, CH₄, C₂H₂, C₂H₄, C₂H₆,
méthanol, éthanol, éthanal, acide éthanoïque), les métaux usuels avec leurs oxydes,
sulfures et sels (Na, K, Mg, Ca, Al, Si, Ti, Cr, Mn, Fe, Ni, Cu, Zn, Ag, Sn, Pb, Hg),
et 28 espèces en solution aqueuse (H⁺ par convention, OH⁻, Na⁺, Cl⁻, SO₄²⁻, Fe²⁺,
Fe³⁺, MnO₄⁻, CH₃COO⁻…) à l'état standard c° = 1 mol/L. Les valeurs viennent des tables
NBS (Wagman et al., 1982) et, pour FeO, des tables NIST-JANAF (Chase, 1998) ; la
plupart des ions n'ont pas de C°p dans ces tables (`Cp0` vaut `None`, et
`depuis_tables(kirchhoff=True)` le dit). `src/pychimie/donnees/SOURCES.md` donne les
références, les conventions, ce qui a été écarté faute de valeur sûre, et la licence.
Le fichier est écrit par `scripts/ecrire_donnees_thermo.py`.

## Éléments, formules et tableau périodique

### Les 118 éléments

`element()` rend un élément par son numéro atomique, son symbole ou son nom français
(sans tenir compte de la casse ni des accents) ; `tableau_periodique()` rend la table
entière, que l'on filtre ou dont on extrait une colonne.

```python
from pychimie import element, tableau_periodique

fer = element("fer")  # ou element(26), element("Fe")
fer.masse_molaire  # 55.845 (g/mol, CIAAW)
fer.theta_fusion  # 1538.0 (°C) ; fer.T_fusion en K
fer.famille, fer.configuration  # 'métal de transition', '[Ar] 3d6 4s2'

t = tableau_periodique()
[e.symbole for e in t.filtrer(famille="halogène")]
t.valeurs("electronegativite")  # {'H': 2.2, 'He': None, …}
```

Une grandeur a une seule unité (g/mol, pm, eV, K, kg/m³) et vaut `None` quand elle
est inconnue. La table est figée dans `src/pychimie/donnees/elements.csv` ; ses
sources et leurs licences — masses de la CIAAW, propriétés du paquet mendeleev
(MIT), aucune donnée de Wikipédia — sont détaillées dans
`src/pychimie/donnees/SOURCES.md`.

### Formules brutes et préparation de solutions

```python
from pychimie import Formule, Preparation, masse_a_peser, masse_molaire

masse_molaire("CuSO4·5H2O")  # 249.68 g/mol
Formule.lire("K4[Fe(CN)6]").hill  # 'C6FeK4N6'
Formule.lire("SO₄²⁻").charge  # -2
Formule.lire("Cu^2+(aq)").phase  # 'aq'
masse_a_peser("CuSO4·5H2O", 0.10, 0.1)  # 2.50 g pour 100 mL à 0,10 mol/L
Preparation.commerciale("H2SO4", 0.96, 1.84).concentration  # 18.0 mol/L
```

Les hydrates (`·`, `.`, `*` ou la virgule française : `CuSO4, 5 H2O`), les parenthèses et crochets imbriqués, les indices et
exposants Unicode sont lus ; une charge s'écrit `^2-`, `²⁻` ou, pour ±1, d'un signe
seul (`NH4+`) ; une phase finale `(s)`, `(l)`, `(g)` ou `(aq)` est gardée à part
(`Formule.phase`). `Preparation` (masse à peser, volume à prélever, solution du
commerce) n'est pas la `Solution` de `titrages`, qui décrit une solution par ses
couples acide-base : les équilibres en solution (pH, titrages) relèvent de ce paquet.

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

## Licence

Le code — le paquet (`src/`), les tests et les scripts — est sous licence MIT
(`LICENSE`). Les exemples (`exemples/`), qui sont du contenu scientifique et
pédagogique plutôt que des algorithmes, sont sous licence Creative Commons
Attribution - Pas d'utilisation commerciale - Partage dans les mêmes conditions 4.0
International, CC BY-NC-SA 4.0 (`LICENSE-CC-BY-NC-SA-4.0`).

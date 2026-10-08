# Changelog

Le format suit [Keep a Changelog](https://keepachangelog.com/fr/1.1.0/) et le
[versionnage sémantique](https://semver.org/lang/fr/).

## [Unreleased]

### Ajouté

- Table des 118 éléments reconstruite (`src/pychimie/donnees/elements.csv`) à partir
  des masses atomiques standard abrégées de la CIAAW (table 2024 : IUPAC 2021 et
  révisions de Gd, Lu, Zr) et du paquet mendeleev (MIT), avec les noms français ;
  sources et licences dans `donnees/SOURCES.md`, scripts de génération dans
  `scripts/`. mendeleev n'est qu'une dépendance du groupe `donnees`.
- `pychimie.elements` : `Element`, `TableauPeriodique`, `element()`,
  `tableau_periodique()` ; accès par Z, symbole ou nom, `filtrer`, `valeurs`,
  `normaliser` ; températures en °C (`theta_fusion`…), caractère métallique.
- `pychimie.formules` : analyseur de formules brutes (parenthèses imbriquées,
  hydrates, charges, indices et exposants Unicode), masse molaire, fractions
  massiques, notation de Hill.
- `pychimie.composes` : `Compose`, une formule et sa phase (« H2O(l) »).
- `pychimie.solutions` : masse à peser, volume à prélever, solution du commerce
  (titre massique et densité).
- `pychimie.periodikz` : `TableauTikz` (cases, ombrage, symboles en hyperliens
  vers Wikipédia FR, informations, légende).
- `pychimie.pgfperiodictable` : schémas `\pgfPTnewColorScheme` pour le paquet
  LaTeX pgf-PeriodicTable.
- Exemples dans `exemples/periodikz` et `exemples/formules`.

### Corrigé

Par rapport à `chimie/atomes.py`, `chimie/periodikz.py` et `chimie/atomes.csv` :

- masses atomiques périmées (Zn, Ge, Se, Mo, Cd, Yb…), Dy (160,5 au lieu de
  162,50), masses des éléments superlourds présentées comme standard ;
- numéro CAS de He (celui de l'hydrogène atomique), et tous les numéros CAS,
  perdus par `int()` ;
- températures de 1000 °C et plus, perdues à cause du séparateur de milliers ;
  As, dont l'ébullition précédait la fusion (il se sublime) ;
- masses volumiques en g/L, g/cm³ et kg/m³ mêlées et moyennées : une seule
  unité, le kg/m³ ;
- configuration électronique de Nb ; H classé parmi les alcalins ;
- `__gt__` toujours vrai (maximum et ombrage faux), `NotImplementedError`
  rendue au lieu de `NotImplemented`, `__format__` et f-string fautifs,
  demi-largeur d'une plage calculée après la moyenne, `get_dict` à clé fixe,
  échappements invalides dans les expressions régulières, `open()` sans
  encodage ;
- `periodikz` : valeurs `Decimal('NaN')` non détectées, `except KeyError` pour
  `list.remove`, nombre de colonnes de légende arrondi vers le bas et division
  par zéro ;
- exemple `test-Tikz` : `with_style` au lieu de `attr_style`, sortie
  `simple_taux.tex` au lieu de `simple_metaux.tex`.

### Supprimé

- `chimie/atomes.py`, `chimie/atomes.csv` (données d'origine Wikipédia),
  `chimie/periodikz.py`, `chimie/composes.py`, `chimie/composes.csv`,
  `chimie/solutions.py`, `chimie/test-Tikz/` : repris dans `src/pychimie` et
  `exemples/periodikz`.

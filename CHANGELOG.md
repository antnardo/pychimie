# Changelog

Le format suit [Keep a Changelog](https://keepachangelog.com/fr/1.1.0/) et le
[versionnage sémantique](https://semver.org/lang/fr/).

## [0.1.0] - 2026-10-10

### Ajouté

- Licences : MIT pour le code, CC BY-NC-SA 4.0 pour les exemples.
- Équilibres chimiques, en remplacement de `chimie/thermochimie.py` :
  - `Equation` lit une équation écrite en texte, phases comprises ;
  - `Reaction` : K°(T) en approximation d'Ellingham ou avec Kirchhoff (ΔrC°p
    constant), ΔrH°, ΔrS° et ΔrC°p par la loi de Hess (`depuis_tables`), ou K° à une
    température (`depuis_K0`) ;
  - `Systeme` : avancement, Q(ξ), ΔrG(ξ), G(ξ) et son minimum, état final (équilibre
    ou rupture d'équilibre), à pression fixée ou à volume de gaz fixé, avec solutés,
    mélanges liquides et constituants inactifs ;
  - `perturber` et `balayage` : déplacements d'équilibre par T, P, ajout d'un
    constituant actif ou inactif, proportions initiales.
- Table de ΔfH°, S°m et C°p à 298,15 K pour 26 espèces (tables NBS 1982, NIST-JANAF
  1998 pour FeO), écrite par `scripts/ecrire_donnees_thermo.py`, sources et licence
  dans `src/pychimie/donnees/SOURCES.md`.
- `exemples/equilibres/` : les scripts `Optimisation_NH3_*.py` et
  `Optimisation_Boudouard.py` refaits, K°(T) avec Kirchhoff, déplacements
  d'équilibre, rupture d'équilibre du carbonate de calcium.
- Table des 118 éléments reconstruite (`src/pychimie/donnees/elements.csv`) à partir
  des masses atomiques standard abrégées de la CIAAW (table 2024 : IUPAC 2021 et
  révisions de Gd, Lu, Zr) et du paquet mendeleev (MIT), avec les noms français ;
  sources et licences dans `donnees/SOURCES.md`, scripts de génération dans
  `scripts/`. mendeleev n'est qu'une dépendance du groupe `donnees`.
- `pychimie.elements` : `Element`, `TableauPeriodique`, `element()`,
  `tableau_periodique()` ; accès par Z, symbole ou nom, `filtrer`, `valeurs`,
  `normaliser` ; températures en °C (`theta_fusion`…), caractère métallique. Un
  élément inconnu est signalé par l'entrée telle qu'on l'a écrite, avec la règle
  de recherche.
- `pychimie.formules` : analyseur de formules brutes (parenthèses imbriquées,
  hydrates avec `·`, `.`, `*` ou la virgule française, charges, indices et exposants
  Unicode, phase finale `(s)`, `(l)`, `(g)` ou `(aq)` gardée dans `Formule.phase`),
  masse molaire, fractions massiques, notation de Hill.
- `pychimie.preparation` : `Preparation` (masse à peser, volume à prélever, solution
  du commerce : titre massique et densité) — nommée ainsi pour ne pas se confondre
  avec `titrages.Solution`, qui décrit une solution par ses couples acide-base.
- `pychimie.periodikz` : `TableauTikz` (cases, ombrage, symboles en hyperliens
  vers Wikipédia FR, informations, légende).
- `pychimie.pgfperiodictable` : schémas `\pgfPTnewColorScheme` pour le paquet
  LaTeX pgf-PeriodicTable.
- Exemples dans `exemples/periodikz` et `exemples/formules`.

### Corrigé

Par rapport à `chimie/thermochimie.py` :

- `R_GP` en `Decimal` faisait planter `Ko`, `xeq`, `DrG` et `compute_G` : R vient de
  `scipy.constants`, en flottant.
- `compute_G` rangeait dans un tableau le tuple que rend `quad` ; G(ξ) est désormais en
  forme close, de dérivée exactement ΔrG(ξ).
- `is_rupture_min` n'était jamais appelé (`is_rupture` testait deux fois le réactif
  limitant) : l'état final traite les deux bornes par la même règle.
- Les solutés se reconnaissaient par `phase == 'aq'`, malgré les suffixes annoncés ;
  ils se reconnaissent par le préfixe `aq`.
- ΔrH° = 1 J/mol et ΔrS° = 1 J/K/mol par défaut : ces deux grandeurs sont
  obligatoires.
- `xeq` échouait (« bisect ne fonctionne pas ») en cas de rupture d'équilibre ;
  `xi_eq` lève `RuptureEquilibre` et `etat_final` rend l'état de rupture.
- `xeq` s'arrêtait à 2·10⁻¹² mol près ; ξ_éq est à la précision de la machine, et les
  quantités résiduelles restent exactes quand K° vaut 10²⁴ ou 10⁻²⁰.

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

- `exemples/titrages/` et `tests/test_migration_titrages.py`, qui ne concernaient que
  le paquet `titrages` : déposés dans son dépôt (`examples/07_…` à `11_…`,
  `tests/test_migration_chimie.py`), avec la dépendance de développement `titrages`.
- `__all__` ne reprend plus `ARTICLES_WIKIPEDIA`, `CARACTERES`, `PHASES`, `RHO_EAU`
  ni `DonneesThermo`, que rien n'appelait : ils restent dans leurs modules.
- Les sorties TikZ des exemples (`exemples/periodikz/simple_metaux.tex`,
  `schemas_pgfpt.tex`) ne sont plus suivies par git : les scripts les régénèrent.
- `chimie/thermochimie.py`.
- `chimie/atomes.py`, `chimie/atomes.csv` (données d'origine Wikipédia),
  `chimie/periodikz.py`, `chimie/composes.py`, `chimie/composes.csv`,
  `chimie/solutions.py`, `chimie/test-Tikz/` : repris dans `src/pychimie` et
  `exemples/periodikz`.

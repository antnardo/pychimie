# Changelog

Le format suit [Keep a Changelog](https://keepachangelog.com/fr/1.1.0/) et le
[versionnage sémantique](https://semver.org/lang/fr/).

## [Unreleased]

### Ajouté

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

### Corrigé (par rapport à `chimie/thermochimie.py`)

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

### Supprimé

- `chimie/thermochimie.py`.

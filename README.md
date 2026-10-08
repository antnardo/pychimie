# pychimie

Thermochimie (avancement, équilibre, rupture), données atomiques et tableau périodique en TikZ

Code d'origine importé tel quel le 8 octobre 2026 ; voir `../AUDIT.md`.

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

`Systeme(reaction, initial, T=…, P=…, V=…, V_gaz=…, inertes=…)` : `P` est la pression
imposée (bar), `V` le volume de la solution (L) quand il y a des solutés, `V_gaz` le
volume du gaz s'il est fixé à la place de P, `inertes` les constituants inactifs
(`{"Ar(g)": 1}`).

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

`src/pychimie/donnees/thermo298.csv` donne ΔfH°, S°m et C°p à 298,15 K pour les 26
espèces des exercices d'équilibre (N₂, H₂, O₂, NH₃, H₂O, CO, CO₂, C, CH₄, CH₃OH, SO₂,
SO₃, NO₂, N₂O₄, CaCO₃, CaO, Fe et ses oxydes, Cu, CuO, Zn, ZnO). Les valeurs viennent des
tables NBS (Wagman et al., 1982) et, pour FeO, des tables NIST-JANAF (Chase, 1998) ;
`src/pychimie/donnees/SOURCES.md` donne les références et la licence. Le fichier est
écrit par `scripts/ecrire_donnees_thermo.py`.

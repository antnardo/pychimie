# Données thermodynamiques à 298,15 K

`thermo298.csv` donne, pour chaque espèce, l'enthalpie standard de formation ΔfH°
(kJ/mol), l'entropie molaire standard S°m (J/K/mol) et la capacité thermique molaire
standard à pression constante C°p,m (J/K/mol), à 298,15 K sous P° = 1 bar. Le fichier
est écrit par `scripts/ecrire_donnees_thermo.py` : on corrige le script et on le
relance, jamais le CSV à la main.

## Espèces

Celles des exercices d'équilibre de CPGE : N₂, H₂, O₂, NH₃, H₂O (liquide et gaz), CO,
CO₂, C (graphite), CH₄, CH₃OH (liquide et gaz), SO₂, SO₃, NO₂, N₂O₄, CaCO₃ (calcite),
CaO, Fe, FeO, Fe₂O₃, Fe₃O₄, Cu, CuO, Zn, ZnO. La colonne `phase` vaut `g`, `l` ou `s` ;
`C(s)` est le graphite, état standard de référence du carbone.

## Sources

| Code | Référence | Espèces |
| --- | --- | --- |
| `NBS82` | D. D. Wagman, W. H. Evans, V. B. Parker, R. H. Schumm, I. Halow, S. M. Bailey, K. L. Churney, R. L. Nuttall, *The NBS tables of chemical thermodynamic properties. Selected values for inorganic and C1 and C2 organic substances in SI units*, J. Phys. Chem. Ref. Data 11, suppl. 2 (1982) | toutes sauf FeO |
| `JANAF98` | M. W. Chase Jr., *NIST-JANAF Thermochemical Tables*, 4e édition, J. Phys. Chem. Ref. Data, monographie 9 (1998) | FeO |

Les tables NBS sont celles que reprennent les manuels (et les annexes de l'Atkins) ; les
scripts de cours en tirent leurs valeurs — pour l'équilibre de Boudouard, ΔrH° =
172,5 kJ/mol et ΔrS° = 175,8 J/K/mol. Elles ne donnent pour FeO que ΔfH° : la ligne
entière vient alors de JANAF, pour ne pas mêler deux tables dans une même espèce. D'une
table à l'autre, les valeurs diffèrent couramment de quelques dixièmes de kJ/mol et de
quelques dixièmes de J/K/mol (CODATA 1989 donne par exemple S°(CaO) = 38,1 J/K/mol, la
table NBS 39,75).

## Licence

Les valeurs numériques sont des faits mesurés, que le droit d'auteur ne protège pas.
Les 26 lignes retenues, trois grandeurs chacune, ne sont pas une partie substantielle
des tables d'origine au sens du droit des bases de données (articles L.342-1 et
suivants du code de la propriété intellectuelle, directive 96/9/CE). Le fichier
`thermo298.csv` se diffuse donc librement avec `pychimie`, sous la licence du dépôt,
à condition de citer les sources ci-dessus.

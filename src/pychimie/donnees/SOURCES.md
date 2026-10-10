# Sources et licences des données

Les données de pychimie sont des faits — masses, températures, rayons, noms des
éléments —, relevés dans des sources à licence permissive ou publiés par l'IUPAC.
Aucune n'est reprise de Wikipédia, dont le texte et les tableaux sont sous
CC BY-SA : l'ancienne table `chimie/atomes.csv`, qui en semblait extraite, a été
retirée le 8 octobre 2026, et seuls les noms français des éléments en ont été
gardés.

## Fichiers

| Fichier | Livré dans le paquet | Contenu | Produit par |
| --- | --- | --- | --- |
| `src/pychimie/donnees/elements.csv` | oui | 118 éléments, 26 colonnes | `scripts/generer_elements.py` |
| `scripts/sources/noms_fr.csv` | non | Z, symbole, nom français | extraction unique de l'ancien `atomes.csv` |
| `scripts/sources/ciaaw_2024_abrege.csv` | non | masses atomiques standard abrégées | `scripts/telecharger_ciaaw.py` |
| `src/pychimie/donnees/thermo298.csv` | oui | ΔfH°, S°m, C°p à 298,15 K, 127 espèces | `scripts/ecrire_donnees_thermo.py` |
| `src/pychimie/donnees/LICENSE-mendeleev.txt` | oui | licence MIT de mendeleev | copie de sa distribution |

Pour régénérer la table :

```bash
uv run python scripts/telecharger_ciaaw.py            # réseau : table CIAAW en vigueur
uv run --group donnees python scripts/generer_elements.py
```

## Sources

### Masses atomiques — CIAAW (IUPAC)

- Source : [table des masses atomiques standard abrégées](https://www.ciaaw.org/abridged-atomic-weights.htm)
  de la Commission on Isotopic Abundances and Atomic Weights, relevée le
  8 octobre 2026 : « Abridged Standard Atomic Weights 2024 ».
- C'est la table IUPAC 2021 (Prohaska *et al.*, *Pure Appl. Chem.* 94, 573, 2022),
  avec les révisions de 2024 (Gd, Lu, Zr).
- Valeurs abrégées à cinq chiffres significatifs, avec leur incertitude : celles du
  tableau périodique de l'IUPAC. Elles donnent une valeur unique aux douze éléments
  dont la masse standard est un intervalle (H, Li, B, C, N, O, Mg, Si, S, Cl, Br, Tl).
- Licence : des valeurs de mesure, qui ne sont pas protégées par le droit d'auteur ;
  source citée.

### Propriétés numériques — mendeleev

- Source : le paquet [mendeleev](https://github.com/lmmentel/mendeleev) 1.3.0,
  licence MIT, « Copyright (c) 2015 Lukasz Mentel ». Il n'est qu'une dépendance de
  développement (groupe `donnees`) : pychimie lit la table figée, sans lui.
- On n'en prend que des valeurs numériques isolées et des identifiants, jamais de
  texte (ni description, ni usages, ni origine du nom).
- Ses propres références, d'après ses métadonnées :

| Colonne | Grandeur | Référence citée par mendeleev |
| --- | --- | --- |
| `configuration` | configuration électronique de l'état fondamental | — |
| `electronegativite` | électronégativité de Pauling | CRC Handbook (Haynes 2014) ; Kr (3,00), absent de mendeleev, repris du CRC |
| `rayon_atomique` | rayon atomique empirique, pm | Slater 1964 |
| `rayon_covalent` | rayon covalent, pm | Cordero *et al.* 2008 |
| `rayon_vdw` | rayon de van der Waals, pm | CRC Handbook (Haynes 2014) |
| `energie_ionisation` | première énergie d'ionisation, eV | NIST Atomic Spectra Database |
| `T_fusion`, `T_ebullition`, `T_sublimation` | sous 101,325 kPa, K | CRC Handbook (Haynes 2016) |
| `masse_volumique` | kg/m³ ; gaz parfait à 25 °C pour un gaz | CRC Handbook (Haynes 2014), Wikipédia EN |
| `nombres_oxydation` | nombres d'oxydation principaux | Wikipédia EN |
| `isotope_reference` | isotope de plus longue période | NUBASE2020 (Kondev 2021) |
| `cas` | numéro CAS du corps simple | — |
| `decouverte` | année de la découverte | — |

Deux colonnes renvoient en partie à Wikipédia en anglais. Ce sont des listes de
valeurs factuelles (un nombre d'oxydation, une masse volumique), que l'on trouve dans
tous les ouvrages de référence, et non une reprise de texte ou de mise en forme.

### Noms français

- `noms_fr.csv` : nom français de chaque élément, extrait par script de l'ancien
  `atomes.csv`. Un nom d'élément est un fait (la nomenclature française officielle),
  non une œuvre.

### Règles écrites dans `generer_elements.py`

Ces colonnes ne viennent d'aucune table : elles se déduisent de Z ou des colonnes
précédentes, par des règles commentées dans le script.

| Colonne | Règle |
| --- | --- |
| `periode`, `groupe`, `bloc` | tableau à 18 colonnes ; Sc, Y, Lu, Lr au groupe 3 (rapport IUPAC 2021) ; La–Yb et Ac–No au bloc f, sans groupe |
| `famille` | classement de mendeleev traduit ; La–Lu lanthanides et Ac–Lr actinides (IUPAC) |
| `etat` | à 25 °C sous 1 atm, d'après les températures de changement d'état |
| `radioactif` | aucun isotope stable : Tc, Pm et Z ≥ 83 |
| `origine` | primordial, issu de désintégrations, ou synthétique |

Pour At, Fr et les éléments de Z ≥ 100, dont on n'a jamais eu d'échantillon
macroscopique, les températures, la masse volumique et l'état restent vides : les
sources n'en donnent que des estimations.

## Licence de la table

`elements.csv` est une compilation de faits, distribuée avec pychimie. Elle reprend
des valeurs de mendeleev : sa licence MIT demande de joindre la mention de copyright
et la notice de permission à toute partie substantielle. Elles sont dans
`LICENSE-mendeleev.txt`, à côté de la table.

## Liens vers Wikipédia

`pychimie.periodikz` produit des hyperliens vers les articles de Wikipédia en
français. Un lien n'est pas une reprise de contenu : seul le titre de l'article
(le nom de l'élément, ou « Mercure (chimie) » et « Tantale (chimie) ») figure dans
le code.

## Données thermodynamiques à 298,15 K (`thermo298.csv`)

`thermo298.csv` donne, pour chaque espèce, l'enthalpie standard de formation ΔfH°
(kJ/mol), l'entropie molaire standard S°m (J/K/mol) et la capacité thermique molaire
standard à pression constante C°p,m (J/K/mol), à 298,15 K sous P° = 1 bar. Le fichier
est écrit par `scripts/ecrire_donnees_thermo.py` : on corrige le script et on le
relance, jamais le CSV à la main.

### Espèces

127 espèces, dont 13 sans C°p :

- gaz simples et composés de l'azote, du soufre et du phosphore : N₂, H₂, O₂, Ar, He,
  NH₃, NO, NO₂, N₂O, N₂O₄, HNO₃(l), NH₄Cl, NH₄NO₃, urée, S, H₂S, SO₂, SO₃, H₂SO₄(l),
  P, PCl₃, PCl₅ ;
- eau (liquide et vapeur), halogènes et hydracides : F₂, Cl₂, Br₂ (liquide et gaz),
  I₂ (solide et gaz), HF, HCl, HBr, HI ;
- carbone et composés organiques jusqu'à C₂, où s'arrête la table NBS : C (graphite),
  CO, CO₂, CH₄, C₂H₂, C₂H₄, C₂H₆, CH₃OH et C₂H₅OH (liquide et gaz), CH₃CHO(g),
  CH₃COOH(l) ;
- métaux, oxydes, sulfures et sels : Na, NaCl, NaOH, Na₂CO₃, NaHCO₃, K, KCl, Mg, MgO,
  Ca, CaO, Ca(OH)₂, CaCO₃, CaCl₂, CaSO₄, BaSO₄, Al, Al₂O₃, Si, SiO₂, Ti, TiO₂, Cr,
  Cr₂O₃, Mn, MnO, MnO₂, Fe, FeO, Fe₂O₃, Fe₃O₄, FeS₂, Ni, NiO, Cu, Cu₂O, CuO, CuSO₄,
  CuSO₄·5H₂O, Zn, ZnO, ZnS, Ag, Ag₂O, AgCl, Sn, SnO₂, Pb, PbO, PbO₂, Hg, HgO ;
- espèces en solution aqueuse : H⁺, OH⁻, Na⁺, K⁺, Mg²⁺, Ca²⁺, Ba²⁺, Al³⁺, Fe²⁺, Fe³⁺,
  Cu²⁺, Zn²⁺, Ag⁺, Pb²⁺, NH₄⁺, F⁻, Cl⁻, Br⁻, I⁻, NO₃⁻, SO₄²⁻, CO₃²⁻, HCO₃⁻, MnO₄⁻,
  CH₃COO⁻, CH₃COOH, NH₃, CO₂.

Conventions :

- une seule variété par espèce condensée, celle de l'état standard de référence ou
  la plus courante : graphite pour C (le diamant aurait la même clé `C(s)`), soufre
  rhombique, phosphore blanc, étain blanc, quartz α, rutile, corindon, anhydrite,
  blende, pyrite, PbO jaune, HgO rouge ;
- les espèces `aq` sont à l'état standard de la table NBS, la solution idéale
  hypothétique à 1 mol/kg, que le cours prend pour c° = 1 mol/L (l'écart, RT ln ρ_eau,
  vaut 7 J/mol) ; les grandeurs des ions reposent sur la convention
  ΔfH° = S° = C°p = 0 pour H⁺(aq) ;
- la table NBS ne donne pas de C°p pour la plupart des ions : la cellule est vide,
  `Cp0` vaut `None`, et `Reaction.depuis_tables(kirchhoff=True)` refuse l'espèce.

Les formules se lisent avec `pychimie.formules` et s'écrivent comme dans une équation
(`SO4^2-(aq)`, `Fe^3+(aq)`, `CuSO4·5H2O(s)`).

Écartées faute d'une valeur sûre dans ces éditions : propane, butane, benzène et
glucose (absents de la table NBS, qui s'arrête aux C₂ ; les valeurs de Pedley 1986 et
des tables TRC diffèrent de l'ordre du kJ/mol), éthanal liquide (sans C°p dans la table
NBS).

### Sources

| Code | Référence | Espèces |
| --- | --- | --- |
| `NBS82` | D. D. Wagman, W. H. Evans, V. B. Parker, R. H. Schumm, I. Halow, S. M. Bailey, K. L. Churney, R. L. Nuttall, *The NBS tables of chemical thermodynamic properties. Selected values for inorganic and C1 and C2 organic substances in SI units*, J. Phys. Chem. Ref. Data 11, suppl. 2 (1982) | toutes sauf FeO (126 espèces) |
| `JANAF98` | M. W. Chase Jr., *NIST-JANAF Thermochemical Tables*, 4e édition, J. Phys. Chem. Ref. Data, monographie 9 (1998) | FeO |

Les tables NBS sont celles que reprennent les manuels (et les annexes de l'Atkins) ; les
scripts de cours en tirent leurs valeurs — pour l'équilibre de Boudouard, ΔrH° =
172,5 kJ/mol et ΔrS° = 175,8 J/K/mol. Elles ne donnent pour FeO que ΔfH° : la ligne
entière vient alors de JANAF, pour ne pas mêler deux tables dans une même espèce. D'une
table à l'autre, les valeurs diffèrent couramment de quelques dixièmes de kJ/mol et de
quelques dixièmes de J/K/mol (les *CODATA key values* de 1989 — J. D. Cox, D. D.
Wagman, V. A. Medvedev, *CODATA Key Values for Thermodynamics*, Hemisphere, 1989 —
donnent par exemple S°(CaO) = 38,1 J/K/mol, la table NBS 39,75) ; on n'a pas mêlé les
deux.

Vérification par la loi de Hess (`tests/test_donnees.py`) : H₂ + Cl₂ → 2 HCl, −184,6
kJ/mol ; thermite 2 Al + Fe₂O₃, −851,5 kJ/mol ; dissolution de NaCl, +3,9 kJ/mol ;
neutralisation H⁺ + OH⁻, −55,8 kJ/mol ; PCl₅ → PCl₃ + Cl₂, 87,9 kJ/mol et
170,3 J/K/mol ; pKs(AgCl) = 9,75, pKe = 14,00, pKs(CaCO₃) = 8,3 à 298 K.

### Licence

Les valeurs numériques sont des faits mesurés, que le droit d'auteur ne protège pas.
Les 127 lignes retenues, trois grandeurs chacune, ne sont pas une partie substantielle
des tables d'origine au sens du droit des bases de données (articles L.342-1 et
suivants du code de la propriété intellectuelle, directive 96/9/CE). Le fichier
`thermo298.csv` se diffuse donc librement avec `pychimie`, sous la licence du dépôt,
à condition de citer les sources ci-dessus.

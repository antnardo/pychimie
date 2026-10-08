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
| `electronegativite` | électronégativité de Pauling | CRC Handbook (Haynes 2014) |
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

"""Écrit ``src/pychimie/donnees/thermo298.csv`` : ΔfH°, S°m et C°p à 298,15 K et 1 bar.

Les valeurs sont recopiées ici, une espèce par ligne, avec la table d'où elles viennent :
le CSV ne se corrige pas à la main, on corrige ce script et on le relance. Il ne garde
que les espèces des exercices d'équilibre de CPGE ; on en ajoute une en ajoutant une
ligne, avec sa source.

Sources (voir ``src/pychimie/donnees/SOURCES.md``) :

- ``NBS82`` : D. D. Wagman et al., *The NBS tables of chemical thermodynamic
  properties*, J. Phys. Chem. Ref. Data 11, suppl. 2 (1982) — la table des manuels, et
  celle d'où viennent les ΔrH° et ΔrS° des scripts de cours (Boudouard : 172,5 kJ/mol
  et 175,8 J/K/mol) ;
- ``JANAF98`` : M. W. Chase, *NIST-JANAF Thermochemical Tables*, 4e éd., J. Phys. Chem.
  Ref. Data, monographie 9 (1998), pour FeO, dont la table NBS ne donne ni S° ni C°p.

    uv run python scripts/ecrire_donnees_thermo.py
"""

import csv
from pathlib import Path

SORTIE = Path(__file__).resolve().parents[1] / "src" / "pychimie" / "donnees" / "thermo298.csv"

ENTETE = ("formule", "phase", "nom", "DfH0_kJ_mol", "Sm0_J_K_mol", "Cp0_J_K_mol", "source")

# formule, phase, nom, ΔfH° (kJ/mol), S°m (J/K/mol), C°p,m (J/K/mol), source
ESPECES = (
    ("N2", "g", "diazote", 0.0, 191.61, 29.125, "NBS82"),
    ("H2", "g", "dihydrogène", 0.0, 130.684, 28.824, "NBS82"),
    ("O2", "g", "dioxygène", 0.0, 205.138, 29.355, "NBS82"),
    ("NH3", "g", "ammoniac", -46.11, 192.45, 35.06, "NBS82"),
    ("H2O", "l", "eau liquide", -285.830, 69.91, 75.291, "NBS82"),
    ("H2O", "g", "eau vapeur", -241.818, 188.825, 33.577, "NBS82"),
    ("CO", "g", "monoxyde de carbone", -110.525, 197.674, 29.142, "NBS82"),
    ("CO2", "g", "dioxyde de carbone", -393.509, 213.74, 37.11, "NBS82"),
    ("C", "s", "carbone graphite", 0.0, 5.740, 8.527, "NBS82"),
    ("CH4", "g", "méthane", -74.81, 186.264, 35.309, "NBS82"),
    ("CH3OH", "l", "méthanol liquide", -238.66, 126.8, 81.6, "NBS82"),
    ("CH3OH", "g", "méthanol gazeux", -200.66, 239.81, 43.89, "NBS82"),
    ("SO2", "g", "dioxyde de soufre", -296.830, 248.22, 39.87, "NBS82"),
    ("SO3", "g", "trioxyde de soufre", -395.72, 256.76, 50.67, "NBS82"),
    ("NO2", "g", "dioxyde d'azote", 33.18, 240.06, 37.20, "NBS82"),
    ("N2O4", "g", "tétraoxyde de diazote", 9.16, 304.29, 77.28, "NBS82"),
    ("CaCO3", "s", "carbonate de calcium (calcite)", -1206.92, 92.9, 81.88, "NBS82"),
    ("CaO", "s", "oxyde de calcium", -635.09, 39.75, 42.80, "NBS82"),
    ("Fe", "s", "fer α", 0.0, 27.28, 25.10, "NBS82"),
    ("FeO", "s", "oxyde de fer(II)", -272.04, 60.75, 49.92, "JANAF98"),
    ("Fe2O3", "s", "oxyde de fer(III) (hématite)", -824.2, 87.40, 103.85, "NBS82"),
    ("Fe3O4", "s", "magnétite", -1118.4, 146.4, 143.43, "NBS82"),
    ("Cu", "s", "cuivre", 0.0, 33.150, 24.435, "NBS82"),
    ("CuO", "s", "oxyde de cuivre(II)", -157.3, 42.63, 42.30, "NBS82"),
    ("Zn", "s", "zinc", 0.0, 41.63, 25.40, "NBS82"),
    ("ZnO", "s", "oxyde de zinc", -348.28, 43.64, 40.25, "NBS82"),
)


def main() -> None:
    SORTIE.parent.mkdir(parents=True, exist_ok=True)
    with SORTIE.open("w", encoding="utf-8", newline="") as f:
        ecrivain = csv.writer(f, lineterminator="\n")
        ecrivain.writerow(ENTETE)
        ecrivain.writerows(ESPECES)
    print(f"{SORTIE} : {len(ESPECES)} espèces")


if __name__ == "__main__":
    main()

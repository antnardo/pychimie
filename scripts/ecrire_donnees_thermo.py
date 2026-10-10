"""Écrit ``src/pychimie/donnees/thermo298.csv`` : ΔfH°, S°m et C°p à 298,15 K et 1 bar.

Les valeurs sont recopiées ici, une espèce par ligne, avec la table d'où elles viennent :
le CSV ne se corrige pas à la main, on corrige ce script et on le relance. Les espèces
sont celles que rencontrent les sujets de CPGE : corps simples, oxydes et sels des
métaux usuels, halogènes et hydracides, petits composés du carbone (la table NBS
s'arrête aux C₂), ions en solution aqueuse. On en ajoute une en ajoutant une ligne,
avec sa source ; une valeur dont on n'est pas sûr ne s'ajoute pas.

Sources (voir ``src/pychimie/donnees/SOURCES.md``) :

- ``NBS82`` : D. D. Wagman et al., *The NBS tables of chemical thermodynamic
  properties*, J. Phys. Chem. Ref. Data 11, suppl. 2 (1982) — la table des manuels, et
  celle d'où viennent les ΔrH° et ΔrS° des scripts de cours (Boudouard : 172,5 kJ/mol
  et 175,8 J/K/mol) ;
- ``JANAF98`` : M. W. Chase, *NIST-JANAF Thermochemical Tables*, 4e éd., J. Phys. Chem.
  Ref. Data, monographie 9 (1998), pour FeO, dont la table NBS ne donne ni S° ni C°p.

Conventions :

- phase ``s``, ``l``, ``g`` : corps pur sous P° = 1 bar ; une seule variété par
  espèce (graphite pour C, soufre rhombique, phosphore blanc, étain blanc, quartz α,
  rutile, corindon, anhydrite, blende, pyrite, PbO jaune) ;
- phase ``aq`` : soluté dans l'eau, à l'état standard de la table NBS (solution idéale
  hypothétique à 1 mol/kg), que le cours de CPGE prend pour c° = 1 mol/L ; les
  grandeurs des ions reposent sur la convention ΔfH°(H⁺, aq) = S°(H⁺, aq) =
  C°p(H⁺, aq) = 0 ;
- C°p vide (``None``) : la table NBS ne le donne pas pour la plupart des ions ;
  ``Reaction.depuis_tables(kirchhoff=True)`` refuse alors l'espèce.

    uv run python scripts/ecrire_donnees_thermo.py
"""

import csv
from pathlib import Path

SORTIE = Path(__file__).resolve().parents[1] / "src" / "pychimie" / "donnees" / "thermo298.csv"

ENTETE = ("formule", "phase", "nom", "DfH0_kJ_mol", "Sm0_J_K_mol", "Cp0_J_K_mol", "source")

Ligne = tuple[str, str, str, float, float, float | None, str]

# formule, phase, nom, ΔfH° (kJ/mol), S°m (J/K/mol), C°p,m (J/K/mol) ou None, source
ESPECES: tuple[Ligne, ...] = (
    # --- gaz simples et composés de l'azote, du soufre, du phosphore ---------------------
    ("N2", "g", "diazote", 0.0, 191.61, 29.125, "NBS82"),
    ("H2", "g", "dihydrogène", 0.0, 130.684, 28.824, "NBS82"),
    ("O2", "g", "dioxygène", 0.0, 205.138, 29.355, "NBS82"),
    ("Ar", "g", "argon", 0.0, 154.843, 20.786, "NBS82"),
    ("He", "g", "hélium", 0.0, 126.150, 20.786, "NBS82"),
    ("NH3", "g", "ammoniac", -46.11, 192.45, 35.06, "NBS82"),
    ("NO", "g", "monoxyde d'azote", 90.25, 210.761, 29.844, "NBS82"),
    ("NO2", "g", "dioxyde d'azote", 33.18, 240.06, 37.20, "NBS82"),
    ("N2O", "g", "protoxyde d'azote", 82.05, 219.85, 38.45, "NBS82"),
    ("N2O4", "g", "tétraoxyde de diazote", 9.16, 304.29, 77.28, "NBS82"),
    ("HNO3", "l", "acide nitrique liquide", -174.10, 155.60, 109.87, "NBS82"),
    ("NH4Cl", "s", "chlorure d'ammonium", -314.43, 94.6, 84.1, "NBS82"),
    ("NH4NO3", "s", "nitrate d'ammonium", -365.56, 151.08, 139.3, "NBS82"),
    ("CO(NH2)2", "s", "urée", -333.51, 104.60, 93.14, "NBS82"),
    ("S", "s", "soufre rhombique", 0.0, 31.80, 22.64, "NBS82"),
    ("H2S", "g", "sulfure d'hydrogène", -20.63, 205.79, 34.23, "NBS82"),
    ("SO2", "g", "dioxyde de soufre", -296.830, 248.22, 39.87, "NBS82"),
    ("SO3", "g", "trioxyde de soufre", -395.72, 256.76, 50.67, "NBS82"),
    ("H2SO4", "l", "acide sulfurique liquide", -813.989, 156.904, 138.91, "NBS82"),
    ("P", "s", "phosphore blanc", 0.0, 41.09, 23.840, "NBS82"),
    ("PCl3", "g", "trichlorure de phosphore", -287.0, 311.78, 71.84, "NBS82"),
    ("PCl5", "g", "pentachlorure de phosphore", -374.9, 364.6, 112.8, "NBS82"),
    # --- eau, halogènes et hydracides --------------------------------------------------
    ("H2O", "l", "eau liquide", -285.830, 69.91, 75.291, "NBS82"),
    ("H2O", "g", "eau vapeur", -241.818, 188.825, 33.577, "NBS82"),
    ("F2", "g", "difluor", 0.0, 202.78, 31.30, "NBS82"),
    ("Cl2", "g", "dichlore", 0.0, 223.066, 33.907, "NBS82"),
    ("Br2", "l", "dibrome liquide", 0.0, 152.231, 75.689, "NBS82"),
    ("Br2", "g", "dibrome gazeux", 30.907, 245.463, 36.02, "NBS82"),
    ("I2", "s", "diiode solide", 0.0, 116.135, 54.438, "NBS82"),
    ("I2", "g", "diiode gazeux", 62.438, 260.69, 36.90, "NBS82"),
    ("HF", "g", "fluorure d'hydrogène", -271.1, 173.779, 29.133, "NBS82"),
    ("HCl", "g", "chlorure d'hydrogène", -92.307, 186.908, 29.12, "NBS82"),
    ("HBr", "g", "bromure d'hydrogène", -36.40, 198.695, 29.142, "NBS82"),
    ("HI", "g", "iodure d'hydrogène", 26.48, 206.594, 29.158, "NBS82"),
    # --- carbone et petits composés organiques (la table NBS s'arrête aux C₂) ---------
    ("C", "s", "carbone graphite", 0.0, 5.740, 8.527, "NBS82"),
    ("CO", "g", "monoxyde de carbone", -110.525, 197.674, 29.142, "NBS82"),
    ("CO2", "g", "dioxyde de carbone", -393.509, 213.74, 37.11, "NBS82"),
    ("CH4", "g", "méthane", -74.81, 186.264, 35.309, "NBS82"),
    ("C2H2", "g", "éthyne (acétylène)", 226.73, 200.94, 43.93, "NBS82"),
    ("C2H4", "g", "éthène (éthylène)", 52.26, 219.56, 43.56, "NBS82"),
    ("C2H6", "g", "éthane", -84.68, 229.60, 52.63, "NBS82"),
    ("CH3OH", "l", "méthanol liquide", -238.66, 126.8, 81.6, "NBS82"),
    ("CH3OH", "g", "méthanol gazeux", -200.66, 239.81, 43.89, "NBS82"),
    ("C2H5OH", "l", "éthanol liquide", -277.69, 160.7, 111.46, "NBS82"),
    ("C2H5OH", "g", "éthanol gazeux", -235.10, 282.70, 65.44, "NBS82"),
    ("CH3CHO", "g", "éthanal gazeux", -166.19, 250.3, 57.3, "NBS82"),
    ("CH3COOH", "l", "acide éthanoïque liquide", -484.5, 159.8, 124.3, "NBS82"),
    # --- métaux, oxydes, sulfures, sels ------------------------------------------------
    ("Na", "s", "sodium", 0.0, 51.21, 28.24, "NBS82"),
    ("NaCl", "s", "chlorure de sodium", -411.153, 72.13, 50.50, "NBS82"),
    ("NaOH", "s", "hydroxyde de sodium", -425.609, 64.455, 59.54, "NBS82"),
    ("Na2CO3", "s", "carbonate de sodium", -1130.68, 134.98, 112.30, "NBS82"),
    ("NaHCO3", "s", "hydrogénocarbonate de sodium", -950.81, 101.7, 87.61, "NBS82"),
    ("K", "s", "potassium", 0.0, 64.18, 29.58, "NBS82"),
    ("KCl", "s", "chlorure de potassium", -436.747, 82.59, 51.30, "NBS82"),
    ("Mg", "s", "magnésium", 0.0, 32.68, 24.89, "NBS82"),
    ("MgO", "s", "oxyde de magnésium", -601.70, 26.94, 37.15, "NBS82"),
    ("Ca", "s", "calcium", 0.0, 41.42, 25.31, "NBS82"),
    ("CaO", "s", "oxyde de calcium", -635.09, 39.75, 42.80, "NBS82"),
    ("Ca(OH)2", "s", "hydroxyde de calcium", -986.09, 83.39, 87.49, "NBS82"),
    ("CaCO3", "s", "carbonate de calcium (calcite)", -1206.92, 92.9, 81.88, "NBS82"),
    ("CaCl2", "s", "chlorure de calcium", -795.8, 104.6, 72.59, "NBS82"),
    ("CaSO4", "s", "sulfate de calcium (anhydrite)", -1434.11, 106.7, 99.66, "NBS82"),
    ("BaSO4", "s", "sulfate de baryum", -1473.2, 132.2, 101.75, "NBS82"),
    ("Al", "s", "aluminium", 0.0, 28.33, 24.35, "NBS82"),
    ("Al2O3", "s", "alumine (corindon)", -1675.7, 50.92, 79.04, "NBS82"),
    ("Si", "s", "silicium", 0.0, 18.83, 20.00, "NBS82"),
    ("SiO2", "s", "silice (quartz α)", -910.94, 41.84, 44.43, "NBS82"),
    ("Ti", "s", "titane", 0.0, 30.63, 25.02, "NBS82"),
    ("TiO2", "s", "dioxyde de titane (rutile)", -944.7, 50.33, 55.48, "NBS82"),
    ("Cr", "s", "chrome", 0.0, 23.77, 23.35, "NBS82"),
    ("Cr2O3", "s", "oxyde de chrome(III)", -1139.7, 81.2, 118.7, "NBS82"),
    ("Mn", "s", "manganèse", 0.0, 32.01, 26.32, "NBS82"),
    ("MnO", "s", "oxyde de manganèse(II)", -385.22, 59.71, 45.44, "NBS82"),
    ("MnO2", "s", "dioxyde de manganèse", -520.03, 53.05, 54.14, "NBS82"),
    ("Fe", "s", "fer α", 0.0, 27.28, 25.10, "NBS82"),
    ("FeO", "s", "oxyde de fer(II)", -272.04, 60.75, 49.92, "JANAF98"),
    ("Fe2O3", "s", "oxyde de fer(III) (hématite)", -824.2, 87.40, 103.85, "NBS82"),
    ("Fe3O4", "s", "magnétite", -1118.4, 146.4, 143.43, "NBS82"),
    ("FeS2", "s", "pyrite", -178.2, 52.93, 62.17, "NBS82"),
    ("Ni", "s", "nickel", 0.0, 29.87, 26.07, "NBS82"),
    ("NiO", "s", "oxyde de nickel(II)", -239.7, 37.99, 44.31, "NBS82"),
    ("Cu", "s", "cuivre", 0.0, 33.150, 24.435, "NBS82"),
    ("Cu2O", "s", "oxyde de cuivre(I)", -168.6, 93.14, 63.64, "NBS82"),
    ("CuO", "s", "oxyde de cuivre(II)", -157.3, 42.63, 42.30, "NBS82"),
    ("CuSO4", "s", "sulfate de cuivre(II) anhydre", -771.36, 109.0, 100.0, "NBS82"),
    ("CuSO4·5H2O", "s", "sulfate de cuivre(II) pentahydraté", -2279.65, 300.4, 280.0, "NBS82"),
    ("Zn", "s", "zinc", 0.0, 41.63, 25.40, "NBS82"),
    ("ZnO", "s", "oxyde de zinc", -348.28, 43.64, 40.25, "NBS82"),
    ("ZnS", "s", "sulfure de zinc (blende)", -205.98, 57.7, 46.0, "NBS82"),
    ("Ag", "s", "argent", 0.0, 42.55, 25.351, "NBS82"),
    ("Ag2O", "s", "oxyde d'argent", -31.05, 121.3, 65.86, "NBS82"),
    ("AgCl", "s", "chlorure d'argent", -127.068, 96.2, 50.79, "NBS82"),
    ("Sn", "s", "étain blanc", 0.0, 51.55, 26.99, "NBS82"),
    ("SnO2", "s", "dioxyde d'étain (cassitérite)", -580.7, 52.3, 52.59, "NBS82"),
    ("Pb", "s", "plomb", 0.0, 64.81, 26.44, "NBS82"),
    ("PbO", "s", "oxyde de plomb(II) jaune", -217.32, 68.70, 45.77, "NBS82"),
    ("PbO2", "s", "dioxyde de plomb", -277.4, 68.6, 64.64, "NBS82"),
    ("Hg", "l", "mercure", 0.0, 76.02, 27.983, "NBS82"),
    ("HgO", "s", "oxyde de mercure(II) rouge", -90.83, 70.29, 44.06, "NBS82"),
    # --- espèces en solution aqueuse (convention H⁺ : 0) -------------------------------
    ("H+", "aq", "ion hydrogène (convention)", 0.0, 0.0, 0.0, "NBS82"),
    ("OH-", "aq", "ion hydroxyde", -229.994, -10.75, -148.5, "NBS82"),
    ("Na+", "aq", "ion sodium", -240.12, 59.0, 46.4, "NBS82"),
    ("K+", "aq", "ion potassium", -252.38, 102.5, 21.8, "NBS82"),
    ("Mg^2+", "aq", "ion magnésium", -466.85, -138.1, None, "NBS82"),
    ("Ca^2+", "aq", "ion calcium", -542.83, -53.1, None, "NBS82"),
    ("Ba^2+", "aq", "ion baryum", -537.64, 9.6, None, "NBS82"),
    ("Al^3+", "aq", "ion aluminium", -531.0, -321.7, None, "NBS82"),
    ("Fe^2+", "aq", "ion fer(II)", -89.1, -137.7, None, "NBS82"),
    ("Fe^3+", "aq", "ion fer(III)", -48.5, -315.9, None, "NBS82"),
    ("Cu^2+", "aq", "ion cuivre(II)", 64.77, -99.6, None, "NBS82"),
    ("Zn^2+", "aq", "ion zinc", -153.89, -112.1, 46.0, "NBS82"),
    ("Ag+", "aq", "ion argent", 105.579, 72.68, 21.8, "NBS82"),
    ("Pb^2+", "aq", "ion plomb(II)", -1.7, 10.5, None, "NBS82"),
    ("NH4+", "aq", "ion ammonium", -132.51, 113.4, 79.9, "NBS82"),
    ("F-", "aq", "ion fluorure", -332.63, -13.8, -106.7, "NBS82"),
    ("Cl-", "aq", "ion chlorure", -167.159, 56.5, -136.4, "NBS82"),
    ("Br-", "aq", "ion bromure", -121.55, 82.4, -141.8, "NBS82"),
    ("I-", "aq", "ion iodure", -55.19, 111.3, -142.3, "NBS82"),
    ("NO3-", "aq", "ion nitrate", -205.0, 146.4, -86.6, "NBS82"),
    ("SO4^2-", "aq", "ion sulfate", -909.27, 20.1, -293.0, "NBS82"),
    ("CO3^2-", "aq", "ion carbonate", -677.14, -56.9, None, "NBS82"),
    ("HCO3-", "aq", "ion hydrogénocarbonate", -691.99, 91.2, None, "NBS82"),
    ("MnO4-", "aq", "ion permanganate", -541.4, 191.2, -82.0, "NBS82"),
    ("CH3COO-", "aq", "ion éthanoate", -486.01, 86.6, -6.3, "NBS82"),
    ("CH3COOH", "aq", "acide éthanoïque dissous", -485.76, 178.7, None, "NBS82"),
    ("NH3", "aq", "ammoniac dissous", -80.29, 111.3, None, "NBS82"),
    ("CO2", "aq", "dioxyde de carbone dissous", -413.80, 117.6, None, "NBS82"),
)


def main() -> None:
    SORTIE.parent.mkdir(parents=True, exist_ok=True)
    with SORTIE.open("w", encoding="utf-8", newline="") as f:
        ecrivain = csv.writer(f, lineterminator="\n")
        ecrivain.writerow(ENTETE)
        for ligne in ESPECES:
            ecrivain.writerow("" if v is None else v for v in ligne)
    sans_cp = sum(1 for ligne in ESPECES if ligne[5] is None)
    print(f"{SORTIE} : {len(ESPECES)} espèces, dont {sans_cp} sans C°p")


if __name__ == "__main__":
    main()

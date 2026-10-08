"""Construit ``src/pychimie/donnees/elements.csv``, la table figée des 118 éléments.

La table est générée une fois pour toutes et versionnée : pychimie la lit avec le
seul module ``csv`` de la bibliothèque standard, sans dépendre de mendeleev (et de
pandas, SQLAlchemy, pint…) à l'exécution. On ne relance ce script que pour suivre
une nouvelle table de la CIAAW ou une nouvelle version de mendeleev.

Trois sources, toutes à licence permissive ou de faits — le détail colonne par
colonne est dans ``src/pychimie/donnees/SOURCES.md`` :

- ``scripts/sources/noms_fr.csv`` : les noms français des éléments ;
- ``scripts/sources/ciaaw_<année>_abrege.csv`` : les masses atomiques standard
  abrégées de la CIAAW (``telecharger_ciaaw.py``) ;
- le paquet mendeleev (licence MIT), pour les autres propriétés numériques.

Le reste — place dans le tableau, bloc, famille, état à 25 °C, origine, radioactivité
— se déduit de Z ou des données précédentes par des règles écrites ci-dessous, pour
que chaque choix de classement soit lisible et discutable ici plutôt qu'enfoui dans
une table.

Usage : ``uv run --group donnees python scripts/generer_elements.py``.
"""

import csv
import math
import re
import sys
from pathlib import Path

from mendeleev.fetch import fetch_table

RACINE = Path(__file__).resolve().parent.parent
SOURCES = RACINE / "scripts" / "sources"
SORTIE = RACINE / "src" / "pychimie" / "donnees" / "elements.csv"

COLONNES = (
    "Z",
    "symbole",
    "nom",
    "periode",
    "groupe",
    "bloc",
    "famille",
    "masse_molaire",
    "incertitude_masse",
    "isotope_reference",
    "configuration",
    "electronegativite",
    "rayon_atomique",
    "rayon_covalent",
    "rayon_vdw",
    "energie_ionisation",
    "T_fusion",
    "T_ebullition",
    "T_sublimation",
    "masse_volumique",
    "etat",
    "nombres_oxydation",
    "radioactif",
    "origine",
    "decouverte",
    "cas",
)

# Fin de chaque période : Z du gaz noble qui la ferme.
FIN_PERIODE = (2, 10, 18, 36, 54, 86, 118)

# Familles de mendeleev, traduites. Le classement de mendeleev est celui des tableaux
# usuels ; il range Po parmi les métalloïdes et At parmi les halogènes.
FAMILLES = {
    "Nonmetals": "non-métal",
    "Noble gases": "gaz noble",
    "Alkali metals": "métal alcalin",
    "Alkaline earth metals": "métal alcalino-terreux",
    "Metalloids": "métalloïde",
    "Halogens": "halogène",
    "Poor metals": "métal pauvre",
    "Transition metals": "métal de transition",
    "Lanthanides": "lanthanide",
    "Actinides": "actinide",
}

# Éléments dont aucun échantillon macroscopique n'a jamais été obtenu : leurs
# températures de changement d'état et leur masse volumique ne sont que des
# estimations, et leur état à 25 °C n'est pas connu. On les laisse vides plutôt que
# de mêler prédictions et mesures dans une même colonne.
SANS_ECHANTILLON = frozenset({85, 87, *range(100, 119)})

T_AMBIANTE = 298.15  # K, 25 °C : l'état « standard » des tableaux.

# Électronégativité de Pauling absente de mendeleev mais donnée par sa source, le
# CRC Handbook (table « Electronegativity of the Elements ») : le krypton, que KrF2
# rend utile en classe. Les estimations des lanthanides et actinides isolés ne
# sont pas ajoutées.
PAULING_COMPLEMENTS = {36: 3.00}

ANNEE = 3.15576e7  # s, année julienne
DUREES = {
    "ysec": 1e-24,
    "zsec": 1e-21,
    "asec": 1e-18,
    "psec": 1e-12,
    "nsec": 1e-9,
    "usec": 1e-6,
    "msec": 1e-3,
    "sec": 1.0,
    "minute": 60.0,
    "hour": 3600.0,
    "day": 86400.0,
    "year": ANNEE,
    "kyear": 1e3 * ANNEE,
    "Myear": 1e6 * ANNEE,
    "Gyear": 1e9 * ANNEE,
    "Tyear": 1e12 * ANNEE,
    "Pyear": 1e15 * ANNEE,
    "Eyear": 1e18 * ANNEE,
    "Zyear": 1e21 * ANNEE,
    "Yyear": 1e24 * ANNEE,
}


def periode_groupe_bloc(Z: int) -> tuple[int, int | None, str]:
    """Place de Z dans un tableau à 18 colonnes, lanthanides et actinides à part.

    Le groupe 3 reçoit Sc, Y, Lu et Lr, comme le recommande le rapport IUPAC de 2021
    sur la composition du groupe 3 ; La–Yb et Ac–No forment les deux lignes du bloc f,
    sans numéro de groupe.
    """
    periode = next(p for p, fin in enumerate(FIN_PERIODE, start=1) if fin >= Z)
    debut = FIN_PERIODE[periode - 2] + 1 if periode > 1 else 1
    rang = Z - debut  # 0 pour le premier élément de la période
    longueur = FIN_PERIODE[periode - 1] - debut + 1
    if periode == 1:
        return periode, (1 if Z == 1 else 18), "s"
    if rang < 2:
        return periode, rang + 1, "s"
    if longueur == 32 and rang < 16:
        return periode, None, "f"
    groupe = 18 - (longueur - 1 - rang)
    return periode, groupe, "d" if groupe <= 12 else "p"


def origine(Z: int) -> str:
    """Présence dans la nature : depuis la formation de la Terre, par filiation, ou pas.

    Th et U sont primordiaux ; Tc, Pm, Po à Ac, Pa, Np et Pu n'existent dans la nature
    qu'à l'état de traces, produits par désintégration ou capture neutronique.
    """
    if Z in (43, 61) or 84 <= Z <= 89 or Z in (91, 93, 94):
        return "désintégration"
    if Z >= 95:
        return "synthétique"
    return "primordial"


def radioactif(Z: int) -> bool:
    """Aucun isotope stable : Tc, Pm, et tout à partir du bismuth (²⁰⁹Bi se désintègre)."""
    return Z in (43, 61) or Z >= 83


def nombre(x: object) -> float | None:
    """Une valeur de pandas, ``None`` si elle manque."""
    if x is None:
        return None
    x = float(x)  # type: ignore[arg-type]
    return None if math.isnan(x) else x


def texte(x: float | int | str | None) -> str:
    """Écriture d'une cellule : sans artefact d'arrondi binaire, vide si absente."""
    if x is None:
        return ""
    if isinstance(x, bool):
        return "1" if x else "0"
    if isinstance(x, float):
        return format(x, ".10g")
    return str(x)


def configuration(econf: str) -> str:
    """« [Kr] 4d4 5s » → « [Kr] 4d4 5s1 » : mendeleev omet l'exposant 1."""
    return " ".join(
        sous_couche if re.fullmatch(r"\[\w+\]|\d[spdf]\d+", sous_couche) else sous_couche + "1"
        for sous_couche in econf.split()
    )


def changements_etat(lignes: list[dict]) -> tuple[float | None, float | None, float | None]:
    """(T_fusion, T_ebullition, T_sublimation) à 1 atm, en K.

    Quand mendeleev donne plusieurs variétés allotropiques, on garde celle qui fond
    réellement sous 1 atm : ni un point de transition entre variétés solides, ni un
    point de sublimation (P blanc, S monoclinique, Se gris, Sn blanc). Pour C et As,
    qui se subliment sous 1 atm, la température de « fusion » de la source est celle
    du point triple, sous pression : elle n'a pas sa place dans la colonne T_fusion,
    et l'ébullition devient une sublimation.
    """
    if not lignes:
        return None, None, None
    fusions = [r for r in lignes if r["is_transition"] != 1.0 and r["is_sublimation_point"] != 1.0]
    if fusions:
        r = fusions[0]
        return nombre(r["melting_point"]), nombre(r["boiling_point"]), None
    return None, None, nombre(lignes[0]["boiling_point"])


def etat(T_fus: float | None, T_eb: float | None, T_sub: float | None) -> str | None:
    """État du corps simple à 25 °C sous 1 atm : s, l ou g, ``None`` si inconnu."""
    if T_eb is not None and T_eb < T_AMBIANTE:
        return "g"
    if T_fus is not None and T_fus < T_AMBIANTE:
        return "l"
    if T_fus is not None or T_sub is not None:
        return "s"
    return None


def lire(chemin: Path) -> dict[int, dict[str, str]]:
    with chemin.open(encoding="utf-8", newline="") as f:
        return {int(r["Z"]): r for r in csv.DictReader(f)}


def main() -> None:
    ciaaw_fichiers = sorted(SOURCES.glob("ciaaw_*_abrege.csv"))
    if not ciaaw_fichiers:
        sys.exit("lancer d'abord scripts/telecharger_ciaaw.py")
    ciaaw = lire(ciaaw_fichiers[-1])
    noms = lire(SOURCES / "noms_fr.csv")

    elements = fetch_table("elements").set_index("atomic_number")
    series = fetch_table("series").set_index("id")["name"].to_dict()
    transitions = fetch_table("phasetransitions").to_dict("records")
    oxydation = fetch_table("oxidationstates")
    ionisation = fetch_table("ionizationenergies")
    isotopes = fetch_table("isotopes")

    lignes = []
    for Z in range(1, 119):
        e = elements.loc[Z]
        if noms[Z]["symbole"] != e.symbol or ciaaw[Z]["symbole"] != e.symbol:
            sys.exit(f"symboles discordants pour Z = {Z}")
        periode, groupe, bloc = periode_groupe_bloc(Z)

        famille = FAMILLES[series[e.series_id]]
        if 57 <= Z <= 71:  # les quinze lanthanides, Lu compris (IUPAC)
            famille = "lanthanide"
        elif 89 <= Z <= 103:
            famille = "actinide"

        # Les masses gardent l'écriture de la CIAAW : « 1.0080 » a cinq chiffres
        # significatifs, que le passage par un flottant effacerait.
        masse = ciaaw[Z]["masse"] or None
        isotope = None
        if masse is None:
            # Pas de masse standard : nombre de masse de l'isotope de plus longue
            # période, celui que les tableaux impriment entre crochets.
            connus = isotopes[(isotopes.atomic_number == Z) & isotopes.half_life.notna()]
            isotope = int(
                max(
                    connus.itertuples(),
                    key=lambda i: i.half_life * DUREES[i.half_life_unit],
                ).mass_number
            )

        T_fus, T_eb, T_sub = changements_etat([r for r in transitions if r["atomic_number"] == Z])
        rho = nombre(e.density)
        if Z in SANS_ECHANTILLON:
            T_fus = T_eb = T_sub = rho = None

        etats_ox = sorted(
            int(n)
            for n in oxydation[
                (oxydation.atomic_number == Z) & (oxydation.category == "main")
            ].oxidation_state
        )
        ei = ionisation[(ionisation.atomic_number == Z) & (ionisation.ion_charge == 0)]
        annee = nombre(e.discovery_year)

        lignes.append(
            {
                "Z": Z,
                "symbole": e.symbol,
                "nom": noms[Z]["nom"],
                "periode": periode,
                "groupe": groupe,
                "bloc": bloc,
                "famille": famille,
                "masse_molaire": masse,
                "incertitude_masse": ciaaw[Z]["incertitude"] or None,
                "isotope_reference": isotope,
                "configuration": configuration(e.electronic_configuration),
                "electronegativite": nombre(e.en_pauling) or PAULING_COMPLEMENTS.get(Z),
                "rayon_atomique": nombre(e.atomic_radius),
                "rayon_covalent": nombre(e.covalent_radius_cordero),
                "rayon_vdw": nombre(e.vdw_radius),
                "energie_ionisation": nombre(ei.ionization_energy.iloc[0]) if len(ei) else None,
                "T_fusion": T_fus,
                "T_ebullition": T_eb,
                "T_sublimation": T_sub,
                "masse_volumique": None if rho is None else round(rho * 1000, 6),
                "etat": etat(T_fus, T_eb, T_sub),
                "nombres_oxydation": " ".join(str(n) for n in etats_ox),
                "radioactif": radioactif(Z),
                "origine": origine(Z),
                "decouverte": None if annee is None else int(annee),
                "cas": e.cas,
            }
        )

    with SORTIE.open("w", encoding="utf-8", newline="") as f:
        ecrivain = csv.writer(f, lineterminator="\n")
        ecrivain.writerow(COLONNES)
        for ligne in lignes:
            ecrivain.writerow(texte(ligne[c]) for c in COLONNES)
    print(f"{SORTIE.relative_to(RACINE)} : {len(lignes)} éléments, {len(COLONNES)} colonnes")


if __name__ == "__main__":
    main()

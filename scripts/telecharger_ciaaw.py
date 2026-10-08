"""Relève la table des masses atomiques standard abrégées de la CIAAW.

La Commission on Isotopic Abundances and Atomic Weights (IUPAC) publie sur
https://www.ciaaw.org la table en vigueur ; ce script la recopie dans
``scripts/sources/ciaaw_<année>_abrege.csv``, que ``generer_elements.py`` lit sans
réseau. Les valeurs abrégées (cinq chiffres significatifs, avec leur incertitude)
sont celles que l'IUPAC imprime sur son tableau périodique : ce sont elles qui servent
à calculer une masse molaire en classe. Les douze éléments dont la masse standard est
un intervalle (H, Li, B, C, N, O, Mg, Si, S, Cl, Br, Tl) y reçoivent une valeur
unique, ce qui évite d'avoir à choisir soi-même dans l'intervalle.

Une masse atomique est un fait de mesure : la table n'est pas une œuvre protégée, et
sa source est citée dans ``src/pychimie/donnees/SOURCES.md``.

Usage : ``uv run python scripts/telecharger_ciaaw.py``.
"""

import csv
import re
import sys
import urllib.request
from html.parser import HTMLParser
from pathlib import Path

URL = "https://www.ciaaw.org/abridged-atomic-weights.htm"
DOSSIER = Path(__file__).resolve().parent / "sources"


class _Tableau(HTMLParser):
    """Collecte les cellules de chaque ligne ``<tr>`` d'une page HTML."""

    def __init__(self) -> None:
        super().__init__()
        self.lignes: list[list[str]] = []
        self._ligne: list[str] | None = None
        self._cellule: list[str] | None = None

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag == "tr":
            self._ligne = []
        elif tag in ("td", "th"):
            self._cellule = []

    def handle_endtag(self, tag: str) -> None:
        if tag in ("td", "th") and self._cellule is not None and self._ligne is not None:
            self._ligne.append(" ".join("".join(self._cellule).split()))
            self._cellule = None
        elif tag == "tr" and self._ligne is not None:
            self.lignes.append(self._ligne)
            self._ligne = None

    def handle_data(self, data: str) -> None:
        if self._cellule is not None:
            self._cellule.append(data)


def _valeur(texte: str) -> tuple[str, str]:
    """« 1.0080 ± 0.0002 » → (« 1.0080 », « 0.0002 ») ; « — » → (« », « »)."""
    texte = texte.replace(" ", " ").replace(" ", "")
    if texte in ("—", "-", ""):
        return "", ""
    m = re.fullmatch(r"([0-9.]+)±([0-9.]+)", texte)
    if m is None:
        raise ValueError(f"valeur CIAAW inattendue : {texte!r}")
    return m.group(1), m.group(2)


def main() -> None:
    requete = urllib.request.Request(URL, headers={"User-Agent": "pychimie (enseignement)"})
    with urllib.request.urlopen(requete) as reponse:
        page = reponse.read().decode("utf-8")
    annee = re.search(r"Abridged Standard Atomic Weights (\d{4})", page)
    if annee is None:
        sys.exit("année de la table introuvable dans la page de la CIAAW")
    analyse = _Tableau()
    analyse.feed(page)
    lignes = [ligne for ligne in analyse.lignes if ligne and ligne[0].isdigit()]
    if [int(ligne[0]) for ligne in lignes] != list(range(1, 119)):
        sys.exit("la table de la CIAAW ne couvre pas Z = 1 à 118 dans l'ordre")
    sortie = DOSSIER / f"ciaaw_{annee.group(1)}_abrege.csv"
    with sortie.open("w", encoding="utf-8", newline="") as f:
        ecrivain = csv.writer(f, lineterminator="\n")
        ecrivain.writerow(["Z", "symbole", "masse", "incertitude", "notes"])
        for Z, symbole, _nom, valeur, notes, *_ in lignes:
            ecrivain.writerow([Z, symbole, *_valeur(valeur), notes])
    print(f"{sortie.name} : {len(lignes)} éléments")


if __name__ == "__main__":
    main()

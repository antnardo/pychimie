"""Les 118 éléments chimiques, lus dans une table figée, en français.

La table ``donnees/elements.csv`` est générée par ``scripts/generer_elements.py`` à
partir de sources à licence permissive (masses de la CIAAW, paquet mendeleev) et de
règles écrites — voir ``donnees/SOURCES.md``. Elle se lit ici avec le seul module
``csv`` : pychimie n'emporte pas mendeleev, pandas et SQLAlchemy pour connaître la
masse du fer.

Chaque grandeur a une seule unité, dite dans le nom de la colonne ou ci-dessous, et
une valeur absente vaut ``None`` : l'ancienne table mêlait g/L, g/cm³ et kg/m³ dans
une même colonne, et perdait en silence toute température écrite avec un séparateur
de milliers. La conversion de chaque cellule se déduit de l'annotation du champ de
:class:`Element` : un seul endroit dit ce que contient chaque colonne.
"""

import csv
import dataclasses
import types
import unicodedata
from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, get_args

__all__ = [
    "CARACTERES",
    "FAMILLES",
    "Element",
    "TableauPeriodique",
    "element",
    "tableau_periodique",
]

DONNEES = Path(__file__).resolve().parent / "donnees" / "elements.csv"

FAMILLES = (
    "métal alcalin",
    "métal alcalino-terreux",
    "lanthanide",
    "actinide",
    "métal de transition",
    "métal pauvre",
    "métalloïde",
    "non-métal",
    "halogène",
    "gaz noble",
)
"""Les familles, dans l'ordre où le tableau les présente de gauche à droite."""

CARACTERES = ("métal", "métalloïde", "non-métal")

_NON_METAUX = frozenset({"non-métal", "halogène", "gaz noble"})

T_ZERO_CELSIUS = 273.15  # K


@dataclass(frozen=True, slots=True)
class Element:
    """Un élément chimique et les propriétés de son atome et de son corps simple.

    Unités : masses molaires en g/mol, rayons en pm, énergie d'ionisation en eV,
    températures en K (sous 1 atm), masse volumique en kg/m³ (à 25 °C pour un gaz).
    """

    Z: int
    symbole: str
    nom: str
    periode: int
    groupe: int | None
    """Groupe 1 à 18 ; ``None`` pour La–Yb et Ac–No, rangés dans les lignes du bloc f."""
    bloc: str
    famille: str
    masse_molaire: float | None
    """Masse atomique standard abrégée de la CIAAW ; ``None`` sans masse standard."""
    incertitude_masse: float | None
    isotope_reference: int | None
    """Sans masse standard : nombre de masse de l'isotope de plus longue période."""
    configuration: str
    electronegativite: float | None
    """Électronégativité de Pauling."""
    rayon_atomique: float | None
    rayon_covalent: float | None
    rayon_vdw: float | None
    energie_ionisation: float | None
    T_fusion: float | None
    T_ebullition: float | None
    T_sublimation: float | None
    """Seulement pour C et As, qui se subliment sous 1 atm au lieu de fondre."""
    masse_volumique: float | None
    etat: str | None
    """« s », « l » ou « g » à 25 °C sous 1 atm ; ``None`` sans échantillon connu."""
    nombres_oxydation: tuple[int, ...]
    """Les nombres d'oxydation principaux, dans l'ordre croissant."""
    radioactif: bool
    origine: str
    """« primordial », « désintégration » ou « synthétique »."""
    decouverte: int | None
    """Année de la découverte ; ``None`` pour les éléments connus depuis l'Antiquité."""
    cas: str

    def __str__(self) -> str:
        return f"{self.symbole} ({self.nom.lower()}, Z = {self.Z})"

    @property
    def masse(self) -> float:
        """La masse à employer dans un calcul : standard, ou à défaut celle de l'isotope.

        Pour Tc, Pm et les éléments lourds, le nombre de masse de l'isotope de plus
        longue période tient lieu de masse molaire, comme les tableaux l'impriment
        entre crochets.
        """
        if self.masse_molaire is not None:
            return self.masse_molaire
        assert self.isotope_reference is not None
        return float(self.isotope_reference)

    @property
    def caractere(self) -> str:
        """Caractère métallique : « métal », « métalloïde » ou « non-métal »."""
        if self.famille == "métalloïde":
            return "métalloïde"
        return "non-métal" if self.famille in _NON_METAUX else "métal"

    @property
    def cle(self) -> str:
        """Le nom en minuscules sans accents, pour un nom de fichier ou une clé."""
        return _cle(self.nom)

    @property
    def theta_fusion(self) -> float | None:
        """Température de fusion en °C."""
        return _celsius(self.T_fusion)

    @property
    def theta_ebullition(self) -> float | None:
        """Température d'ébullition en °C."""
        return _celsius(self.T_ebullition)

    @property
    def theta_sublimation(self) -> float | None:
        """Température de sublimation en °C."""
        return _celsius(self.T_sublimation)


def _cle(nom: str) -> str:
    decompose = unicodedata.normalize("NFKD", nom.strip().lower())
    return "".join(c for c in decompose if not unicodedata.combining(c))


def _celsius(T: float | None) -> float | None:
    return None if T is None else round(T - T_ZERO_CELSIUS, 6)


def _convertir(texte: str, annotation: Any) -> Any:
    """Convertit une cellule selon l'annotation du champ qui la reçoit."""
    if annotation is bool:
        return texte == "1"
    if annotation in (int, float, str):
        return annotation(texte)
    if annotation == tuple[int, ...]:
        return tuple(int(n) for n in texte.split())
    options = get_args(annotation)
    if types.NoneType in options:
        if texte == "":
            return None
        (type_,) = (t for t in options if t is not types.NoneType)
        return _convertir(texte, type_)
    raise TypeError(f"annotation non prise en charge : {annotation!r}")


class TableauPeriodique:
    """Les éléments, accessibles par Z, par symbole ou par nom.

    ``tableau[26]``, ``tableau["Fe"]`` et ``tableau["fer"]`` désignent le même élément :
    le symbole se lit tel quel (Co n'est pas CO), le nom sans tenir compte de la casse
    ni des accents.
    """

    def __init__(self, elements: Iterable[Element]) -> None:
        self._elements = tuple(sorted(elements, key=lambda e: e.Z))
        self._index: dict[int | str, Element] = {}
        for e in self._elements:
            self._index[e.Z] = e
            self._index[e.symbole] = e
            self._index[e.cle] = e

    @classmethod
    def depuis_csv(cls, chemin: Path = DONNEES) -> "TableauPeriodique":
        """Lit une table au format de ``donnees/elements.csv``."""
        champs = {f.name: f.type for f in dataclasses.fields(Element)}
        with chemin.open(encoding="utf-8", newline="") as f:
            lignes = list(csv.DictReader(f))
        return cls(
            Element(**{nom: _convertir(ligne[nom], type_) for nom, type_ in champs.items()})
            for ligne in lignes
        )

    def __getitem__(self, cle: int | str) -> Element:
        demande = cle
        if isinstance(cle, str) and cle not in self._index:
            cle = _cle(cle)
        try:
            return self._index[cle]
        except KeyError:
            raise KeyError(
                f"élément inconnu : {demande!r} (un numéro atomique, un symbole avec sa casse "
                "comme « Fe », ou un nom français, sans tenir compte de la casse ni des accents)"
            ) from None

    def __contains__(self, cle: object) -> bool:
        try:
            self[cle]  # type: ignore[index]
        except KeyError:
            return False
        return True

    def __iter__(self) -> Iterator[Element]:
        return iter(self._elements)

    def __len__(self) -> int:
        return len(self._elements)

    def filtrer(
        self, predicat: Callable[[Element], bool] | None = None, **egalites: object
    ) -> list[Element]:
        """Les éléments qui vérifient ``predicat`` et dont les attributs valent ``egalites``.

        >>> [e.symbole for e in tableau_periodique().filtrer(famille="gaz noble", periode=2)]
        ['Ne']
        """
        return [
            e
            for e in self._elements
            if (predicat is None or predicat(e))
            and all(getattr(e, nom) == valeur for nom, valeur in egalites.items())
        ]

    def valeurs(self, attribut: str) -> dict[str, Any]:
        """``{symbole: valeur}`` pour un attribut, valeurs absentes comprises."""
        return {e.symbole: getattr(e, attribut) for e in self._elements}

    def normaliser(
        self, attribut: str, bornes: tuple[float, float] | None = None
    ) -> dict[int, float | None]:
        """``{Z: x}`` avec x entre 0 et 1, affine en la valeur numérique de l'attribut.

        Par défaut, 0 et 1 vont aux valeurs extrêmes de la table, pour que tout le
        dégradé serve ; ``bornes`` fixe l'échelle (par exemple ``(0, 4)`` pour
        l'électronégativité), et les valeurs hors bornes sont ramenées à 0 ou 1. Une
        valeur absente reste ``None``.
        """
        brutes = {e.Z: getattr(e, attribut) for e in self._elements}
        presentes = [v for v in brutes.values() if v is not None]
        if not presentes:
            raise ValueError(f"aucune valeur pour {attribut!r}")
        if any(isinstance(v, bool) or not isinstance(v, int | float) for v in presentes):
            raise TypeError(f"{attribut!r} n'est pas une grandeur numérique")
        bas, haut = bornes if bornes is not None else (min(presentes), max(presentes))
        if haut <= bas:
            raise ValueError(f"bornes vides pour {attribut!r} : {bas} à {haut}")
        return {
            Z: None if v is None else min(1.0, max(0.0, (v - bas) / (haut - bas)))
            for Z, v in brutes.items()
        }


@lru_cache(maxsize=1)
def tableau_periodique() -> TableauPeriodique:
    """La table livrée avec pychimie, lue une seule fois."""
    return TableauPeriodique.depuis_csv()


def element(cle: int | str) -> Element:
    """Un élément de la table livrée, par Z, symbole ou nom : ``element("Fe")``."""
    return tableau_periodique()[cle]

"""La table des ΔfH°, S°m et C°p,m à 298,15 K, lue dans ``donnees/thermo298.csv``.

La table vit dans un CSV du paquet, lisible et corrigeable sans toucher au code ; elle
est écrite par ``scripts/ecrire_donnees_thermo.py``, et ses sources et sa licence sont
dans ``donnees/SOURCES.md``. Les grandeurs sont rendues en unités SI (J/mol, J/K/mol) :
le kJ/mol des tables ne sort pas de ce module.

Une espèce se cherche par ``"formule(phase)"`` (``"H2O(g)"``, ``"C(s)"`` pour le
graphite, ``"SO4^2-(aq)"`` pour un ion). Une phase étiquetée (``l1``, ``s2``) se cherche
sous sa phase nue : les grandeurs standard sont celles du corps pur. Les espèces ``aq``
sont à l'état standard des tables (solution idéale à 1 mol/kg, prise pour
c° = 1 mol/L), avec la convention ΔfH° = S° = C°p = 0 pour H⁺(aq) ; la plupart des
ions n'ont pas de C°p dans la table : ``Cp0`` vaut alors ``None``.
"""

import csv
from collections.abc import Iterator
from dataclasses import dataclass
from functools import cache
from importlib import resources
from pathlib import Path

from pychimie.equation import Constituant

__all__ = ["DonneesThermo", "TableThermo", "table_thermo"]


@dataclass(frozen=True, slots=True)
class DonneesThermo:
    """Les grandeurs standard d'une espèce à 298,15 K, en J/mol et J/K/mol."""

    formule: str
    phase: str
    nom: str
    DfH0: float
    Sm0: float
    Cp0: float | None
    """``None`` quand la table ne le donne pas (la plupart des ions aqueux)."""
    source: str

    @property
    def cle(self) -> str:
        return f"{self.formule}({self.phase})"


class TableThermo:
    """Un ensemble d'espèces, cherchées par ``"formule(phase)"``."""

    T_REF = 298.15

    def __init__(self, especes: tuple[DonneesThermo, ...]) -> None:
        self._especes = {e.cle: e for e in especes}
        if len(self._especes) != len(especes):
            raise ValueError("la table contient deux fois la même espèce")

    @classmethod
    def lire(cls, chemin: str | Path | None = None) -> "TableThermo":
        """Lit le CSV du paquet, ou celui de ``chemin`` (mêmes colonnes)."""
        if chemin is None:
            texte = resources.files("pychimie").joinpath("donnees/thermo298.csv").read_text("utf-8")
        else:
            texte = Path(chemin).read_text(encoding="utf-8")
        lignes = csv.DictReader(texte.splitlines())
        return cls(
            tuple(
                DonneesThermo(
                    formule=ligne["formule"],
                    phase=ligne["phase"],
                    nom=ligne["nom"],
                    DfH0=float(ligne["DfH0_kJ_mol"]) * 1e3,
                    Sm0=float(ligne["Sm0_J_K_mol"]),
                    Cp0=float(ligne["Cp0_J_K_mol"]) if ligne["Cp0_J_K_mol"] != "" else None,
                    source=ligne["source"],
                )
                for ligne in lignes
            )
        )

    def __getitem__(self, cle: str | Constituant) -> DonneesThermo:
        if isinstance(cle, Constituant):
            cle = f"{cle.formule}({cle.etat})"
        try:
            return self._especes[cle.strip()]
        except KeyError:
            raise KeyError(
                f"{cle} n'est pas dans la table ({len(self)} espèces : {', '.join(self._especes)})"
            ) from None

    def cles(self) -> tuple[str, ...]:
        """Les clés ``"formule(phase)"`` des espèces, dans l'ordre de la table."""
        return tuple(self._especes)

    def __contains__(self, cle: object) -> bool:
        return cle in self._especes

    def __iter__(self) -> Iterator[DonneesThermo]:
        return iter(self._especes.values())

    def __len__(self) -> int:
        return len(self._especes)


@cache
def table_thermo() -> TableThermo:
    """La table du paquet, lue une fois."""
    return TableThermo.lire()

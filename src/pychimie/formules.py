"""Formules brutes : composition, masse molaire, notation de Hill.

L'analyseur est une descente récursive sur la grammaire des formules écrites en
classe :

.. code-block:: text

    formule  := partie (séparateur partie)*        séparateur : · • ⋅ . *
    partie   := coefficient? groupe+ charge?       « 5H2O »
    groupe   := (symbole | (groupe+) | [groupe+] | {groupe+}) nombre?
    charge   := ^? nombre? signe | exposants        « SO4^2- », « NH4+ », « Fe³⁺ »

Une expression régulière ne suffit pas : les parenthèses s'imbriquent
(« K4[Fe(CN)6] »), et c'est leur multiplicateur qui fait la composition. Les indices
et exposants Unicode (« H₂SO₄ », « SO₄²⁻ ») sont acceptés tels qu'on les copie d'un
document.

Une charge sans ``^`` ni exposant ne peut être qu'un signe seul : « Fe3+ » se lirait
Fe₃⁺, et « SO42- » S O₄₂⁻. On écrit donc « Fe^3+ » ou « Fe³⁺ » ; un signe seul
(« NH4+ », « Cl- ») vaut une charge ±1. La masse d'un ion ignore celle des électrons
(0,000 55 g/mol par charge), très en deçà de l'incertitude des masses atomiques.
"""

from collections import Counter
from dataclasses import dataclass

from pychimie.elements import TableauPeriodique, tableau_periodique

__all__ = [
    "AnalyseurFormule",
    "Formule",
    "FormuleInvalide",
    "masse_molaire",
]

SEPARATEURS = frozenset("·•⋅.*")
OUVRANTES = {"(": ")", "[": "]", "{": "}"}
SIGNES = {"+": 1, "-": -1, "−": -1}
INDICES = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")
EXPOSANTS = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻", "0123456789+-")
CARACTERES_EXPOSANTS = frozenset("⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻")
# Pas str.isdigit() : il accepte aussi « ² », qui doit rester un exposant de charge.
CHIFFRES = frozenset("0123456789")
# Ce qui termine une suite de groupes : séparateur, début de charge, espace.
FIN_DE_GROUPE = SEPARATEURS | SIGNES.keys() | CARACTERES_EXPOSANTS | {"^", " "}


class FormuleInvalide(ValueError):  # noqa: N818 — nom français, comme toute l'API
    """La formule ne se lit pas : le message dit où et pourquoi."""


@dataclass(frozen=True, slots=True)
class Formule:
    """Une formule lue : nombre d'atomes de chaque élément et charge.

    ``composition`` garde l'ordre de première apparition des symboles :
    ``Formule.lire("CH3COOH").composition == (("C", 2), ("H", 4), ("O", 2))``.
    """

    texte: str
    composition: tuple[tuple[str, int], ...]
    charge: int = 0

    @classmethod
    def lire(cls, texte: str, tableau: TableauPeriodique | None = None) -> "Formule":
        """Analyse ``texte`` ; lève :class:`FormuleInvalide` si la formule est fautive."""
        return AnalyseurFormule(tableau).analyser(texte)

    def __str__(self) -> str:
        return self.texte

    def nombre(self, symbole: str) -> int:
        """Nombre d'atomes de l'élément ``symbole`` (0 s'il est absent)."""
        return dict(self.composition).get(symbole, 0)

    @property
    def atomes(self) -> int:
        """Nombre total d'atomes."""
        return sum(n for _, n in self.composition)

    @property
    def hill(self) -> str:
        """Formule brute en notation de Hill : C, puis H, puis l'ordre alphabétique.

        Sans carbone, tout est alphabétique, H compris : « H2O4S » pour l'acide
        sulfurique. C'est l'ordre des index de formules (CAS, catalogues).
        """
        compte = dict(self.composition)
        if "C" in compte:
            ordre = ["C", *(["H"] if "H" in compte else [])]
            ordre += sorted(s for s in compte if s not in ("C", "H"))
        else:
            ordre = sorted(compte)
        brute = "".join(s + (str(compte[s]) if compte[s] > 1 else "") for s in ordre)
        return brute + _charge_texte(self.charge)

    def masse_molaire(self, tableau: TableauPeriodique | None = None) -> float:
        """Masse molaire en g/mol, avec les masses atomiques standard abrégées."""
        t = tableau or tableau_periodique()
        return sum(n * t[s].masse for s, n in self.composition)

    def fractions_massiques(self, tableau: TableauPeriodique | None = None) -> dict[str, float]:
        """Part de la masse molaire apportée par chaque élément (somme égale à 1)."""
        t = tableau or tableau_periodique()
        M = self.masse_molaire(t)
        return {s: n * t[s].masse / M for s, n in self.composition}


def _charge_texte(charge: int) -> str:
    if charge == 0:
        return ""
    signe = "+" if charge > 0 else "-"
    return "^" + (str(abs(charge)) if abs(charge) > 1 else "") + signe


class AnalyseurFormule:
    """Lit une formule brute par descente récursive.

    Les symboles sont vérifiés contre la table, casse comprise : « Xy » est refusé,
    et « Co » reste le cobalt — l'oxyde de carbone s'écrit « CO ».
    """

    def __init__(self, tableau: TableauPeriodique | None = None) -> None:
        self._symboles = frozenset(e.symbole for e in (tableau or tableau_periodique()))
        self._texte = ""
        self._pos = 0

    def analyser(self, texte: str) -> Formule:
        """Analyse une formule complète, hydrates et charge compris."""
        self._texte = texte.strip().translate(INDICES)
        self._pos = 0
        if not self._texte:
            raise FormuleInvalide("formule vide")
        total: Counter[str] = Counter()
        charge = 0
        while True:
            coefficient, composition, charge_partie = self._partie()
            for symbole, n in composition.items():
                total[symbole] += coefficient * n
            charge += coefficient * charge_partie
            self._espaces()
            if self._fin():
                break
            if self._courant() not in SEPARATEURS:
                self._erreur("séparateur ou fin de formule attendu")
            self._pos += 1
            self._espaces()
        return Formule(texte=texte.strip(), composition=tuple(total.items()), charge=charge)

    # Grammaire -----------------------------------------------------------------

    def _partie(self) -> tuple[int, Counter[str], int]:
        coefficient = self._entier(defaut=1)
        if coefficient == 0:
            self._erreur("coefficient nul")
        self._espaces()
        composition = self._suite(fermante=None)
        return coefficient, composition, self._charge()

    def _suite(self, fermante: str | None) -> Counter[str]:
        composition: Counter[str] = Counter()
        while not self._fin():
            c = self._courant()
            if c in OUVRANTES:
                self._pos += 1
                interne = self._suite(fermante=OUVRANTES[c])
                if self._fin() or self._courant() != OUVRANTES[c]:
                    self._erreur(f"« {OUVRANTES[c]} » attendue")
                self._pos += 1
                multiplicateur = self._multiplicateur()
                for symbole, n in interne.items():
                    composition[symbole] += multiplicateur * n
            elif c.isupper():
                symbole = self._symbole()
                composition[symbole] += self._multiplicateur()
            elif c in FIN_DE_GROUPE or c == fermante:
                break
            else:
                self._erreur(f"caractère inattendu « {c} »")
        if not composition:
            self._erreur("groupe vide")
        return composition

    def _symbole(self) -> str:
        debut = self._pos
        self._pos += 1
        if not self._fin() and self._courant().islower():
            self._pos += 1
        symbole = self._texte[debut : self._pos]
        if symbole not in self._symboles:
            self._pos = debut
            self._erreur(f"symbole inconnu « {symbole} »")
        return symbole

    def _multiplicateur(self) -> int:
        n = self._entier(defaut=1)
        if n == 0:
            self._erreur("indice nul")
        return n

    def _charge(self) -> int:
        """« ^2- », « ^+ », « + », « ²⁻ », ou rien."""
        if self._fin():
            return 0
        if self._courant() in CARACTERES_EXPOSANTS:
            debut = self._pos
            while not self._fin() and self._courant() in CARACTERES_EXPOSANTS:
                self._pos += 1
            return self._lire_charge(self._texte[debut : self._pos].translate(EXPOSANTS), debut)
        if self._courant() == "^":
            self._pos += 1
            accolade = not self._fin() and self._courant() == "{"
            self._pos += accolade
            debut = self._pos
            while not self._fin() and (self._courant() in CHIFFRES or self._courant() in SIGNES):
                self._pos += 1
            charge = self._lire_charge(self._texte[debut : self._pos], debut)
            if accolade:
                if self._fin() or self._courant() != "}":
                    self._erreur("« } » attendue")
                self._pos += 1
            return charge
        if self._courant() in SIGNES:
            debut = self._pos
            while not self._fin() and self._courant() in SIGNES:
                self._pos += 1
            return self._lire_charge(self._texte[debut : self._pos], debut)
        return 0

    def _lire_charge(self, texte: str, debut: int) -> int:
        """« 2- », « -2 », « + », « ++ » → entier signé."""
        chiffres = "".join(c for c in texte if c in CHIFFRES)
        signes = [SIGNES[c] for c in texte if c in SIGNES]
        valide = signes and len(set(signes)) == 1 and (not chiffres or len(signes) == 1)
        if not valide or texte.strip("+-−").strip("0123456789") != "":
            self._pos = debut
            self._erreur("charge illisible")
        return signes[0] * (int(chiffres) if chiffres else len(signes))

    # Lecture ------------------------------------------------------------------

    def _entier(self, defaut: int) -> int:
        debut = self._pos
        while not self._fin() and self._courant() in CHIFFRES:
            self._pos += 1
        return int(self._texte[debut : self._pos]) if self._pos > debut else defaut

    def _espaces(self) -> None:
        while not self._fin() and self._courant() == " ":
            self._pos += 1

    def _courant(self) -> str:
        return self._texte[self._pos]

    def _fin(self) -> bool:
        return self._pos >= len(self._texte)

    def _erreur(self, message: str) -> None:
        raise FormuleInvalide(f"{message} en position {self._pos + 1} de « {self._texte} »")


def masse_molaire(formule: str | Formule) -> float:
    """Masse molaire en g/mol d'une formule : ``masse_molaire("CuSO4·5H2O")``."""
    if isinstance(formule, str):
        formule = Formule.lire(formule)
    return formule.masse_molaire()

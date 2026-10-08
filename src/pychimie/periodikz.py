r"""Tableau périodique en TikZ : cases, symboles et informations, élément par élément.

Le code produit ne dessine rien de lui-même : il appelle des commandes et des styles
que définit un gabarit (voir ``exemples/periodikz/``). Chaque case est un
``\dessinecase[options]{coordonnée}``, chaque symbole un ``node[symbole, …]``, et la
mise en forme — couleurs, polices, décalages — reste entièrement dans le gabarit.
Ce partage permet de changer l'allure d'un tableau sans relancer Python, et de
produire plusieurs tableaux (familles, électronégativité, rayons) avec un même code.

Coordonnées attendues du gabarit : ``<groupe>-<période>`` pour le tableau principal
(``1-1`` pour H, ``18-7`` pour Og) et ``<rang>-L`` ou ``<rang>-A`` pour les deux
lignes du bloc f, de La (``1-L``) à Yb (``14-L``) et de Ac à No.

Deux façons de colorer une case :

- :meth:`TableauTikz.cases` passe la valeur d'attributs qualitatifs à une clé de style,
  ``famille=halogène`` ou ``etat=g`` ; le gabarit définit ``famille/halogène/.style``
  avec ``/.is choice`` ;
- :meth:`TableauTikz.ombrage` passe une propriété numérique ramenée entre 0 et 100,
  ``chi=73``, qu'un style à argument utilise : ``chi/.style={fill=blue!#1!red}``.

Une valeur absente ne produit aucune clé : la case garde le style par défaut, au lieu
d'une valeur ``no`` que chaque gabarit devait prévoir.
"""

import math
import urllib.parse
from collections.abc import Callable, Iterable, Mapping
from pathlib import Path
from typing import Any

from pychimie.elements import Element, TableauPeriodique, tableau_periodique

__all__ = [
    "ARTICLES_WIKIPEDIA",
    "TableauTikz",
    "coordonnees",
    "url_wikipedia",
]

ARTICLES_WIKIPEDIA = {
    "Mercure": "Mercure (chimie)",
    "Tantale": "Tantale (chimie)",
}
"""Éléments dont le nom désigne, sur Wikipédia FR, une page d'homonymie."""


def coordonnees(e: Element) -> str:
    """Nom de la coordonnée TikZ de la case de ``e`` dans le gabarit."""
    if e.groupe is None:
        rang = e.Z - (56 if e.periode == 6 else 88)
        return f"{rang}-{'L' if e.periode == 6 else 'A'}"
    return f"{e.groupe}-{e.periode}"


def url_wikipedia(e: Element) -> str:
    """Adresse de l'article de Wikipédia en français consacré à l'élément."""
    titre = ARTICLES_WIKIPEDIA.get(e.nom, e.nom).replace(" ", "_")
    return "https://fr.wikipedia.org/wiki/" + urllib.parse.quote(titre)


def _valeur_de_style(valeur: Any) -> str | None:
    """Une valeur d'attribut écrite comme clé pgfkeys ; ``None`` si elle manque."""
    if valeur is None or valeur in ("", ()):
        return None
    if isinstance(valeur, bool):
        return "oui" if valeur else "non"
    return str(valeur)


def _options(*options: str) -> str:
    presentes = [o for o in options if o]
    return f"[{', '.join(presentes)}]" if presentes else ""


Format = str | Callable[[Any], str]


class TableauTikz:
    r"""Accumule les commandes TikZ d'un tableau périodique, couche par couche.

    Chaque méthode ajoute une couche pour les 118 éléments, dans l'ordre où on les
    appelle : les cases d'abord, puis les symboles et les informations par-dessus.

    >>> t = TableauTikz(liens=False).cases("famille").symboles("etat")
    >>> t.texte.splitlines()[0]
    '\\dessinecase[famille=non-métal]{1-1}'
    """

    def __init__(self, tableau: TableauPeriodique | None = None, *, liens: bool = True) -> None:
        self._tableau = tableau or tableau_periodique()
        self.liens = liens
        """Symboles en hyperliens vers Wikipédia (paquet hyperref requis)."""
        self._lignes: list[str] = []

    @property
    def texte(self) -> str:
        """Le code TikZ produit jusqu'ici."""
        return "".join(ligne + "\n" for ligne in self._lignes)

    def ecrire(self, chemin: str | Path) -> Path:
        """Écrit le code dans ``chemin`` (UTF-8) et rend ce chemin."""
        chemin = Path(chemin)
        chemin.write_text(self.texte, encoding="utf-8")
        return chemin

    def cases(self, *attributs: str, style: str = "") -> "TableauTikz":
        r"""Une case par élément, stylée par la valeur de chaque attribut.

        ``cases("famille")`` → ``\dessinecase[famille=gaz noble]{18-1}`` ; un booléen
        s'écrit ``oui`` ou ``non`` (``radioactif=oui``). ``style`` s'ajoute à toutes
        les cases.
        """
        for e in self._tableau:
            options = _options(style, *self._styles(e, attributs))
            self._lignes.append(rf"\dessinecase{options}{{{coordonnees(e)}}}")
        return self

    def ombrage(
        self,
        attribut: str,
        style: str,
        *,
        bornes: tuple[float, float] | None = None,
        defaut: str = "",
    ) -> "TableauTikz":
        r"""Une case par élément, ombrée selon une propriété numérique.

        La valeur est ramenée à un entier de 0 à 100 (voir
        :meth:`TableauPeriodique.normaliser` pour ``bornes``) et passée au style :
        ``ombrage("electronegativite", "chi")`` → ``\dessinecase[chi=100]{17-2}``
        pour le fluor. Un élément sans valeur reçoit le style ``defaut``.
        """
        x = self._tableau.normaliser(attribut, bornes)
        for e in self._tableau:
            v = x[e.Z]
            options = _options(defaut) if v is None else f"[{style}={round(100 * v)}]"
            self._lignes.append(rf"\dessinecase{options}{{{coordonnees(e)}}}")
        return self

    def symboles(
        self,
        *attributs: str,
        style: str = "",
        selection: Iterable[str] | None = None,
    ) -> "TableauTikz":
        r"""Le symbole de chaque élément, au style ``symbole`` du gabarit.

        Les ``attributs`` stylent le symbole comme :meth:`cases` style la case
        (``symboles("etat")`` colore les gaz et les liquides) ; ``style`` s'ajoute
        aux seuls symboles de ``selection`` si elle est donnée, à tous sinon. Chaque
        symbole devient aussi une coordonnée TikZ, ``(Fe)``, pour y pointer une
        flèche.
        """
        choisis = None if selection is None else frozenset(selection)
        for e in self._tableau:
            commun = style if choisis is None or e.symbole in choisis else ""
            options = _options("symbole", commun, *self._styles(e, attributs))
            texte = self._hyperlien(e) if self.liens else e.symbole
            self._lignes.append(
                rf"\draw ({coordonnees(e)}) coordinate({e.symbole}) node{options} {{{texte}}};"
            )
        return self

    def infos(self, attribut: str, style: str, format: Format = "{}") -> "TableauTikz":
        r"""Une information par case : ``node[style]{valeur}``.

        ``format`` est un modèle (``"\\num{{{}}}"``) ou une fonction de la valeur
        (``lambda r: f"\\SI{{{r / 100:g}}}{{\\angstrom}}"`` pour un rayon donné en
        pm). Les éléments sans valeur n'ont pas de nœud.
        """
        for e in self._tableau:
            valeur = getattr(e, attribut)
            if valeur is None or valeur == ():
                continue
            texte = format(valeur) if callable(format) else format.format(valeur)
            self._lignes.append(rf"\draw ({coordonnees(e)}) node[{style}] {{{texte}}};")
        return self

    def legende(
        self,
        attribut: str,
        *,
        colonne: int = 0,
        ligne: int = 0,
        rangs: int = 1,
        horizontale: bool = False,
        valeurs: Iterable[Any] | None = None,
        exclues: Iterable[Any] = (),
        libelles: Mapping[str, str] | None = None,
    ) -> int:
        r"""Une entrée ``\legend`` (ou ``\legendhoriz``) par valeur de l'attribut.

        Les valeurs viennent dans l'ordre du tableau, ou dans celui de ``valeurs``,
        moins ``exclues`` ; elles se répartissent sur ``rangs`` colonnes (rangs
        horizontaux si ``horizontale``), à partir de la position (``colonne``,
        ``ligne``) du gabarit. Le libellé est la valeur avec une majuscule, ou celui
        de ``libelles``, indexé par la valeur telle qu'elle s'écrit dans le style.

        Rend le nombre de positions occupées en largeur, pour placer une légende
        suivante à sa droite.
        """
        if rangs < 1:
            raise ValueError("il faut au moins un rang de légende")
        source = (getattr(e, attribut) for e in self._tableau) if valeurs is None else valeurs
        retirees = {_valeur_de_style(v) for v in exclues}
        entrees = [
            v
            for v in dict.fromkeys(_valeur_de_style(v) for v in source)
            if v is not None and v not in retirees
        ]
        if not entrees:
            return 0
        par_rang = math.ceil(len(entrees) / rangs)
        commande = "legendhoriz" if horizontale else "legend"
        for i, v in enumerate(entrees):
            rang, place = divmod(i, par_rang)
            x, y = (place, rang) if horizontale else (rang, place)
            libelle = (libelles or {}).get(v, v[:1].upper() + v[1:])
            self._lignes.append(
                rf"\{commande}{{{attribut}={v}}}{{{x + colonne}}}{{{y + ligne}}}{{{libelle}}}"
            )
        return par_rang if horizontale else math.ceil(len(entrees) / par_rang)

    # Aides -------------------------------------------------------------------

    @staticmethod
    def _styles(e: Element, attributs: Iterable[str]) -> list[str]:
        styles = []
        for attribut in attributs:
            valeur = _valeur_de_style(getattr(e, attribut))
            if valeur is not None:
                styles.append(f"{attribut}={valeur}")
        return styles

    @staticmethod
    def _hyperlien(e: Element) -> str:
        # Dans un fichier lu par \input, « % » ouvrirait un commentaire : hyperref
        # accepte \% dans l'adresse d'un \href et l'écrit « % » dans le PDF.
        url = url_wikipedia(e).replace("%", r"\%")
        return rf"\href{{{url}}}{{{e.symbole}}}"

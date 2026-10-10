"""Une réaction et sa constante d'équilibre K°(T).

Une réaction est son équation et trois grandeurs standard à la température de
référence : ΔrH°, ΔrS° et ΔrC°p. On les donne (valeurs d'un énoncé), on les tire des
tables par la loi de Hess (:meth:`Reaction.depuis_tables`), ou on part de K° à une
température (:meth:`Reaction.depuis_K0`).

**Approximation d'Ellingham et loi de Kirchhoff.** Avec ΔrC°p constant,

    ΔrH°(T) = ΔrH°(T₀) + ΔrC°p (T − T₀),   ΔrS°(T) = ΔrS°(T₀) + ΔrC°p ln(T/T₀),

et l'approximation d'Ellingham n'est que le cas ΔrC°p = 0 : un seul calcul, pas deux
modèles. K°(T) = exp(−ΔrG°(T)/RT) avec ΔrG° = ΔrH° − TΔrS°.

**Pas de valeur par défaut.** L'ancien ``Reaction`` prenait ΔrH° = 1 J/mol et
ΔrS° = 1 J/K/mol quand on oubliait de les donner, et rendait un K° faux sans rien dire.
Ici ΔrH° et ΔrS° sont des arguments obligatoires ; une réaction dont on ne connaît que
K° se construit par :meth:`Reaction.depuis_K0`, qui dit ce qu'elle suppose.

Les températures sont en kelvins et peuvent être des tableaux : ``r.K0(T)`` se trace.
"""

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike
from scipy.constants import R

from pychimie._calcul import Flottants, sortie, temperatures
from pychimie.donnees import TableThermo, table_thermo
from pychimie.equation import Equation

__all__ = ["Reaction"]


@dataclass(frozen=True, slots=True, init=False)
class Reaction:
    """Une équation, et ΔrH°, ΔrS°, ΔrC°p à ``T_ref`` (J/mol, J/K/mol)."""

    equation: Equation
    DrH0_ref: float
    DrS0_ref: float
    DrCp0: float
    T_ref: float

    def __init__(
        self,
        equation: Equation | str,
        *,
        DrH0: float,
        DrS0: float,
        DrCp0: float = 0.0,
        T_ref: float = 298.15,
    ) -> None:
        """``Reaction("N2(g) + 3 H2(g) = 2 NH3(g)", DrH0=-92.2e3, DrS0=-198.8)``.

        ``DrCp0`` vaut 0 par défaut : c'est l'approximation d'Ellingham, pas une valeur
        inventée.
        """
        if isinstance(equation, str):
            equation = Equation.lire(equation)
        if T_ref <= 0:
            raise ValueError("la température de référence est en kelvins")
        object.__setattr__(self, "equation", equation)
        object.__setattr__(self, "DrH0_ref", float(DrH0))
        object.__setattr__(self, "DrS0_ref", float(DrS0))
        object.__setattr__(self, "DrCp0", float(DrCp0))
        object.__setattr__(self, "T_ref", float(T_ref))

    @classmethod
    def depuis_K0(
        cls, equation: Equation | str, K0: float, *, T: float = 298.15, DrH0: float = 0.0
    ) -> "Reaction":
        """La réaction de constante ``K0`` à ``T``.

        Avec ``DrH0``, K° varie avec T selon la loi de van 't Hoff (ΔrH° constant). Sans,
        ΔrH° = 0 : K° ne dépend pas de la température, et ``K0()`` se passe de T.
        """
        if K0 <= 0:
            raise ValueError("K° est strictement positive")
        DrS0 = R * np.log(K0) + DrH0 / T
        return cls(equation, DrH0=DrH0, DrS0=float(DrS0), T_ref=T)

    @classmethod
    def depuis_tables(
        cls,
        equation: Equation | str,
        *,
        kirchhoff: bool = False,
        table: TableThermo | None = None,
    ) -> "Reaction":
        """ΔrH°, ΔrS° (et ΔrC°p si ``kirchhoff``) par la loi de Hess, à 298,15 K.

        ΔrX° = Σ νᵢ X°ᵢ, avec ΔfH° pour H, S°m pour S et C°p,m pour Cp. Sans
        ``kirchhoff``, ΔrC°p est mis à zéro : approximation d'Ellingham. Avec, chaque
        espèce doit avoir un C°p dans la table (la plupart des ions aqueux n'en ont pas).
        """
        if isinstance(equation, str):
            equation = Equation.lire(equation)
        table = table_thermo() if table is None else table
        donnees = [table[c] for c in equation.constituants]
        nu = equation.nu
        DrCp0 = 0.0
        if kirchhoff:
            sans_cp = [d.cle for d in donnees if d.Cp0 is None]
            if sans_cp:
                raise ValueError(
                    f"pas de C°p dans la table pour {', '.join(sans_cp)} : Kirchhoff impossible, "
                    "rester en approximation d'Ellingham (kirchhoff=False)"
                )
            DrCp0 = _somme(nu, [d.Cp0 for d in donnees if d.Cp0 is not None])
        return cls(
            equation,
            DrH0=_somme(nu, [d.DfH0 for d in donnees]),
            DrS0=_somme(nu, [d.Sm0 for d in donnees]),
            DrCp0=DrCp0,
            T_ref=TableThermo.T_REF,
        )

    @property
    def kirchhoff(self) -> bool:
        """Vrai si ΔrC°p ≠ 0 : ΔrH° et ΔrS° varient avec T."""
        return self.DrCp0 != 0

    def ellingham(self) -> "Reaction":
        """La même réaction en approximation d'Ellingham (ΔrC°p = 0)."""
        return Reaction(self.equation, DrH0=self.DrH0_ref, DrS0=self.DrS0_ref, T_ref=self.T_ref)

    def DrH0(self, T: ArrayLike | None = None) -> Flottants:
        """ΔrH°(T) en J/mol ; à ``T_ref`` sans argument."""
        if T is None:
            return self.DrH0_ref
        return sortie(self.DrH0_ref + self.DrCp0 * (temperatures(T) - self.T_ref))

    def DrS0(self, T: ArrayLike | None = None) -> Flottants:
        """ΔrS°(T) en J/K/mol ; à ``T_ref`` sans argument."""
        if T is None:
            return self.DrS0_ref
        return sortie(self.DrS0_ref + self.DrCp0 * np.log(temperatures(T) / self.T_ref))

    def DrG0(self, T: ArrayLike) -> Flottants:
        """ΔrG°(T) = ΔrH°(T) − TΔrS°(T), en J/mol."""
        T = temperatures(T)
        return sortie(np.asarray(self.DrH0(T)) - T * np.asarray(self.DrS0(T)))

    def ln_K0(self, T: ArrayLike | None = None) -> Flottants:
        """ln K°(T) : on y travaille, K° sort des flottants dès 10³⁰⁸.

        Sans T, seulement si K° ne dépend pas de la température (ΔrH° = 0, ΔrC°p = 0).
        """
        if T is None:
            if self.DrH0_ref == 0 and self.DrCp0 == 0:
                return self.DrS0_ref / R
            raise ValueError("K° dépend de la température (ΔrH° ≠ 0) : donner T")
        T = temperatures(T)
        return sortie(-np.asarray(self.DrG0(T)) / (R * T))

    def K0(self, T: ArrayLike | None = None) -> Flottants:
        """K°(T) = exp(−ΔrG°(T)/RT)."""
        return sortie(np.exp(np.asarray(self.ln_K0(T))))

    def __str__(self) -> str:
        texte = (
            f"{self.equation} : ΔrH° = {self.DrH0_ref / 1e3:.2f} kJ/mol, "
            f"ΔrS° = {self.DrS0_ref:.2f} J/K/mol"
        )
        if self.kirchhoff:
            texte += f", ΔrC°p = {self.DrCp0:.2f} J/K/mol"
        return f"{texte} (à {self.T_ref:g} K)"


def _somme(nu: Sequence[float], valeurs: Sequence[float]) -> float:
    return float(sum(n * v for n, v in zip(nu, valeurs, strict=True)))

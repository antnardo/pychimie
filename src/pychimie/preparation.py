"""Préparer une solution : pesée, dilution, solution du commerce.

Les calculs de la paillasse — quelle masse peser, quel volume prélever, quelle
concentration a l'acide sulfurique à 96 % du flacon — reposent tous sur la masse
molaire du soluté, que :mod:`pychimie.formules` calcule depuis sa formule. Les
équilibres en solution (pH, titrages) relèvent du paquet ``titrages``, dont la classe
``Solution`` décrit une solution par ses couples acide-base : ici une
:class:`Preparation` est une solution à préparer, décrite par son soluté et sa
concentration.

Unités de la chimie des solutions, plutôt que du SI : concentrations en mol/L,
concentrations massiques en g/L, volumes en L, masses en g. La densité est celle
d'un liquide par rapport à l'eau, ρ = d × 1000 g/L.
"""

from dataclasses import dataclass

from pychimie.formules import Formule

__all__ = [
    "RHO_EAU",
    "Preparation",
    "masse_a_peser",
    "volume_a_prelever",
]

RHO_EAU = 1000.0  # g/L : masse volumique de référence de la densité


def _formule(solute: str | Formule) -> Formule:
    return Formule.lire(solute) if isinstance(solute, str) else solute


def masse_a_peser(solute: str | Formule, concentration: float, volume: float) -> float:
    """Masse (g) de solide à dissoudre pour ``volume`` L à ``concentration`` mol/L.

    La formule est celle du solide pesé : un hydrate apporte son eau, et
    ``masse_a_peser("CuSO4·5H2O", 0.1, 0.1)`` donne 2,50 g, contre 1,60 g pour le
    sulfate anhydre.
    """
    return concentration * volume * _formule(solute).masse_molaire()


def volume_a_prelever(
    concentration_mere: float, concentration_fille: float, volume: float
) -> float:
    """Volume (L) de solution mère à diluer pour ``volume`` L de solution fille.

    La quantité de soluté se conserve : c_mère × V_prélevé = c_fille × V_fiole.
    """
    if concentration_fille > concentration_mere:
        raise ValueError("une dilution ne peut pas augmenter la concentration")
    return concentration_fille * volume / concentration_mere


@dataclass(frozen=True, slots=True)
class Preparation:
    """Une solution à préparer : un soluté et sa concentration en quantité de matière.

    ``densite`` n'est connue que pour une solution du commerce : c'est elle qui relie
    le titre massique à la concentration.
    """

    solute: Formule
    concentration: float
    densite: float | None = None

    @classmethod
    def molaire(cls, solute: str | Formule, concentration: float) -> "Preparation":
        """Solution de ``concentration`` mol/L."""
        return cls(_formule(solute), concentration)

    @classmethod
    def massique(cls, solute: str | Formule, concentration_massique: float) -> "Preparation":
        """Solution de concentration massique donnée, en g/L."""
        formule = _formule(solute)
        return cls(formule, concentration_massique / formule.masse_molaire())

    @classmethod
    def commerciale(
        cls, solute: str | Formule, titre_massique: float, densite: float
    ) -> "Preparation":
        """Solution du commerce, d'après l'étiquette : titre massique w et densité d.

        c = w × d × ρ_eau / M : 1 L de solution pèse d × 1000 g, dont la fraction w de
        soluté. ``Preparation.commerciale("H2SO4", 0.96, 1.84)`` : 18,0 mol/L.
        """
        if not 0 < titre_massique <= 1:
            raise ValueError("le titre massique est une fraction, entre 0 et 1")
        formule = _formule(solute)
        c = titre_massique * densite * RHO_EAU / formule.masse_molaire()
        return cls(formule, c, densite)

    @property
    def masse_molaire(self) -> float:
        """Masse molaire du soluté, en g/mol."""
        return self.solute.masse_molaire()

    @property
    def concentration_massique(self) -> float:
        """Concentration massique, en g/L."""
        return self.concentration * self.masse_molaire

    @property
    def titre_massique(self) -> float | None:
        """Fraction massique du soluté ; ``None`` si la densité est inconnue."""
        if self.densite is None:
            return None
        return self.concentration_massique / (self.densite * RHO_EAU)

    def quantite(self, volume: float) -> float:
        """Quantité de soluté (mol) dans ``volume`` L."""
        return self.concentration * volume

    def volume_pour(self, quantite: float) -> float:
        """Volume (L) qui contient ``quantite`` mol de soluté."""
        return quantite / self.concentration

    def diluee(self, facteur: float) -> "Preparation":
        """La solution diluée ``facteur`` fois ; la densité n'est plus connue."""
        if facteur < 1:
            raise ValueError("le facteur de dilution vaut au moins 1")
        return Preparation(self.solute, self.concentration / facteur)

    def masse_a_peser(self, volume: float) -> float:
        """Masse (g) du soluté solide à peser pour ``volume`` L de cette solution."""
        return masse_a_peser(self.solute, self.concentration, volume)

    def volume_a_prelever(self, concentration: float, volume: float) -> float:
        """Volume (L) à prélever de cette solution pour ``volume`` L à ``concentration`` mol/L."""
        return volume_a_prelever(self.concentration, concentration, volume)

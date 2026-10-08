import doctest
import re

import pytest

import pychimie.pgfperiodictable
from pychimie.pgfperiodictable import BLANC, GRIS, rgb, schema_categories, schema_degrade


def triplets(commande: str) -> list[tuple[float, ...]]:
    m = re.fullmatch(r"\\pgfPTnewColorScheme\{(\w+)\}\{(.*)\}\n", commande)
    assert m is not None
    return [tuple(float(c) for c in t.split("/")) for t in m.group(2).split(",")]


def test_doctests_passent():
    assert doctest.testmod(pychimie.pgfperiodictable).failed == 0


class TestRgb:
    @pytest.mark.parametrize(
        ("couleur", "attendu"),
        [("#ff0000", (1.0, 0.0, 0.0)), ("00FF00", (0.0, 1.0, 0.0)), ((0, 0, 1), (0.0, 0.0, 1.0))],
    )
    def test_conversions(self, couleur, attendu):
        assert rgb(couleur) == attendu

    @pytest.mark.parametrize("couleur", ["rouge", "#ff00", (255, 0, 0), (1, 0)])
    def test_couleurs_invalides(self, couleur):
        with pytest.raises(ValueError):
            rgb(couleur)


class TestCategories:
    def test_une_couleur_par_element(self):
        assert len(triplets(schema_categories("etats", "etat", {"g": "#ff0000"}))) == 118

    def test_couleur_par_valeur_et_defaut(self):
        t = triplets(schema_categories("etats", "etat", {"g": "#ff0000", "l": (0, 0, 1)}))
        assert t[0] == (1.0, 0.0, 0.0)  # H
        assert t[34] == (0.0, 0.0, 1.0)  # Br
        assert t[25] == BLANC  # Fe

    def test_booleens(self):
        t = triplets(schema_categories("radio", "radioactif", {True: "#ffff00"}, defaut=GRIS))
        assert t[42] == (1.0, 1.0, 0.0)  # Tc
        assert t[41] == GRIS  # Mo


class TestDegrade:
    def test_extremes_et_manquants(self):
        t = triplets(schema_degrade("chi", "electronegativite", "#ff0000", "#0000ff"))
        assert t[8] == (0.0, 0.0, 1.0)  # F, le plus électronégatif
        assert t[86] == (1.0, 0.0, 0.0)  # Fr
        assert t[1] == GRIS  # He

    def test_bornes(self):
        t = triplets(
            schema_degrade("chi", "electronegativite", (0, 0, 0), (1, 1, 1), bornes=(0, 4))
        )
        assert t[8] == pytest.approx((0.995,) * 3)


def test_nom_en_lettres_seules():
    with pytest.raises(ValueError, match="lettres"):
        schema_categories("mon schéma 2", "etat", {})

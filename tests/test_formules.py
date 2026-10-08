import re

import pytest

from pychimie.elements import TableauPeriodique, element
from pychimie.formules import AnalyseurFormule, Formule, FormuleInvalide, masse_molaire


class TestComposition:
    @pytest.mark.parametrize(
        ("texte", "composition"),
        [
            ("H2O", {"H": 2, "O": 1}),
            ("H2SO4", {"H": 2, "S": 1, "O": 4}),
            ("Ca(OH)2", {"Ca": 1, "O": 2, "H": 2}),
            ("(NH4)2SO4", {"N": 2, "H": 8, "S": 1, "O": 4}),
            ("K4[Fe(CN)6]", {"K": 4, "Fe": 1, "C": 6, "N": 6}),
            ("CuSO4·5H2O", {"Cu": 1, "S": 1, "O": 9, "H": 10}),
            ("CuSO4.5H2O", {"Cu": 1, "S": 1, "O": 9, "H": 10}),
            ("CuSO4 * 5 H2O", {"Cu": 1, "S": 1, "O": 9, "H": 10}),
            ("Al2(SO4)3·18H2O", {"Al": 2, "S": 3, "O": 30, "H": 36}),
            ("CH3COOH", {"C": 2, "H": 4, "O": 2}),
            ("H₂SO₄", {"H": 2, "S": 1, "O": 4}),
            ("2H2O", {"H": 4, "O": 2}),
            ("Co", {"Co": 1}),
            ("CO", {"C": 1, "O": 1}),
        ],
    )
    def test_composition(self, texte, composition):
        assert dict(Formule.lire(texte).composition) == composition

    def test_ordre_de_premiere_apparition(self):
        assert [s for s, _ in Formule.lire("CH3COOH").composition] == ["C", "H", "O"]

    def test_nombre_et_atomes(self):
        f = Formule.lire("C6H12O6")
        assert f.nombre("C") == 6
        assert f.nombre("N") == 0
        assert f.atomes == 24


class TestCharge:
    @pytest.mark.parametrize(
        ("texte", "charge"),
        [
            ("SO4^2-", -2),
            ("SO4^{2-}", -2),
            ("SO₄²⁻", -2),
            ("Fe³⁺", 3),
            ("Fe^3+", 3),
            ("NH4+", 1),
            ("Cl-", -1),
            ("Cl−", -1),
            ("O^-2", -2),
            ("H2O", 0),
            ("2NH4+", 2),
        ],
    )
    def test_charge(self, texte, charge):
        assert Formule.lire(texte).charge == charge

    def test_charge_sans_effet_sur_la_composition(self):
        assert dict(Formule.lire("SO4^2-").composition) == {"S": 1, "O": 4}


class TestMasseMolaire:
    @pytest.mark.parametrize(
        ("texte", "masse"),
        [
            ("H2O", 18.015),
            ("H2SO4", 98.072),
            ("CuSO4·5H2O", 249.677),
            ("Ca(OH)2", 74.092),
            ("NaCl", 58.44),
            ("C6H12O6", 180.156),
            ("KMnO4", 158.032),
        ],
    )
    def test_masses(self, texte, masse):
        assert masse_molaire(texte) == pytest.approx(masse, abs=1e-3)

    def test_formule_deja_lue(self):
        f = Formule.lire("H2O")
        assert masse_molaire(f) == f.masse_molaire()

    def test_element_sans_masse_standard(self):
        assert masse_molaire("TcO4-") == pytest.approx(97 + 4 * 15.999)

    def test_fractions_massiques(self):
        w = Formule.lire("H2O").fractions_massiques()
        assert sum(w.values()) == pytest.approx(1)
        assert w["O"] == pytest.approx(15.999 / 18.015, abs=1e-5)

    def test_tableau_explicite(self):
        t = TableauPeriodique([element("H"), element("O")])
        assert Formule.lire("H2O", t).masse_molaire(t) == pytest.approx(18.015, abs=1e-3)


class TestHill:
    @pytest.mark.parametrize(
        ("texte", "hill"),
        [
            ("CH3COOH", "C2H4O2"),
            ("H2SO4", "H2O4S"),
            ("CCl4", "CCl4"),
            ("NaCl", "ClNa"),
            ("SO4^2-", "O4S^2-"),
            ("NH4+", "H4N^+"),
        ],
    )
    def test_notation_de_hill(self, texte, hill):
        assert Formule.lire(texte).hill == hill


class TestErreurs:
    @pytest.mark.parametrize(
        ("texte", "message"),
        [
            ("", "vide"),
            ("Xy", "symbole inconnu « Xy »"),
            ("Or", "symbole inconnu « Or »"),
            ("h2o", "caractère inattendu"),
            ("(H2O", "« ) » attendue"),
            ("H2O)", "caractère inattendu « ) »"),
            ("H0", "indice nul"),
            ("0H2O", "coefficient nul"),
            ("Fe 2+", "séparateur ou fin"),
            ("Fe^3", "charge illisible"),
            ("SO4^+-", "charge illisible"),
            ("H2O··H2O", "groupe vide"),
            ("()", "groupe vide"),
        ],
    )
    def test_formule_fautive(self, texte, message):
        with pytest.raises(FormuleInvalide, match=re.escape(message)):
            Formule.lire(texte)

    def test_erreur_est_une_valueerror(self):
        assert issubclass(FormuleInvalide, ValueError)

    def test_analyseur_reutilisable(self):
        a = AnalyseurFormule()
        assert a.analyser("H2O").atomes == 3
        assert a.analyser("CO2").atomes == 3

    def test_str_rend_le_texte(self):
        assert str(Formule.lire(" CuSO4·5H2O ")) == "CuSO4·5H2O"

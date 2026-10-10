import pytest

from pychimie.formules import Formule
from pychimie.preparation import Preparation, masse_a_peser, volume_a_prelever


class TestFonctions:
    def test_masse_d_un_hydrate(self):
        assert masse_a_peser("CuSO4·5H2O", 0.1, 0.1) == pytest.approx(2.4968, abs=1e-4)

    def test_masse_du_sel_anhydre(self):
        assert masse_a_peser("CuSO4", 0.1, 0.1) == pytest.approx(1.5961, abs=1e-4)

    def test_volume_a_prelever(self):
        assert volume_a_prelever(1.0, 0.1, 0.1) == pytest.approx(0.01)

    def test_dilution_impossible(self):
        with pytest.raises(ValueError, match="dilution"):
            volume_a_prelever(0.1, 1.0, 0.1)

    def test_volume_en_millilitres_refuse(self):
        with pytest.raises(ValueError, match=r"litres.*100 mL = 0.1 L"):
            masse_a_peser("CuSO4·5H2O", 0.1, 100)

    def test_volume_en_millilitres_refuse_pour_la_dilution(self):
        with pytest.raises(ValueError, match="litres"):
            volume_a_prelever(1.0, 0.1, 250)

    def test_volume_nul_refuse(self):
        with pytest.raises(ValueError, match="positif"):
            masse_a_peser("NaCl", 0.1, 0)


class TestPreparation:
    def test_acide_sulfurique_du_commerce(self):
        s = Preparation.commerciale("H2SO4", 0.96, 1.84)
        assert s.concentration == pytest.approx(18.01, abs=0.01)

    def test_titre_massique_retrouve(self):
        assert Preparation.commerciale("HCl", 0.37, 1.18).titre_massique == pytest.approx(0.37)

    def test_titre_inconnu_sans_densite(self):
        assert Preparation.molaire("NaCl", 0.1).titre_massique is None

    def test_titre_hors_bornes(self):
        with pytest.raises(ValueError, match="fraction"):
            Preparation.commerciale("H2SO4", 96, 1.84)

    def test_concentration_massique(self):
        s = Preparation.massique("NaCl", 5.844)
        assert s.concentration == pytest.approx(0.1, abs=1e-5)
        assert s.concentration_massique == pytest.approx(5.844)

    def test_formule_deja_lue(self):
        assert Preparation.molaire(Formule.lire("NaCl"), 0.1).masse_molaire == pytest.approx(58.44)

    def test_quantite_et_volume(self):
        s = Preparation.molaire("NaCl", 0.2)
        assert s.quantite(0.05) == pytest.approx(0.01)
        assert s.volume_pour(0.01) == pytest.approx(0.05)

    def test_diluee_perd_la_densite(self):
        s = Preparation.commerciale("H2SO4", 0.96, 1.84).diluee(10)
        assert s.concentration == pytest.approx(1.801, abs=1e-3)
        assert s.densite is None

    def test_facteur_de_dilution_invalide(self):
        with pytest.raises(ValueError, match="facteur"):
            Preparation.molaire("NaCl", 0.1).diluee(0.5)

    def test_volume_a_prelever_depuis_la_solution(self):
        s = Preparation.molaire("NaCl", 1.0)
        assert s.volume_a_prelever(0.1, 0.1) == pytest.approx(0.01)

    def test_masse_a_peser_depuis_la_solution(self):
        s = Preparation.molaire("CuSO4·5H2O", 0.1)
        assert s.masse_a_peser(0.1) == pytest.approx(masse_a_peser("CuSO4·5H2O", 0.1, 0.1))

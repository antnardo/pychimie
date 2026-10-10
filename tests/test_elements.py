import doctest
import re
from collections import Counter

import pytest

import pychimie.elements
from pychimie.elements import FAMILLES, Element, TableauPeriodique, element, tableau_periodique


@pytest.fixture(scope="module")
def tableau() -> TableauPeriodique:
    return tableau_periodique()


def test_doctests_passent():
    assert doctest.testmod(pychimie.elements).failed == 0


class TestTable:
    def test_contient_les_118_elements_dans_l_ordre(self, tableau):
        assert [e.Z for e in tableau] == list(range(1, 119))

    def test_symboles_uniques(self, tableau):
        assert len({e.symbole for e in tableau}) == 118

    @pytest.mark.parametrize(
        ("symbole", "masse"),
        [
            ("H", 1.0080),
            ("Fe", 55.845),
            ("Zn", 65.38),
            ("Ge", 72.630),
            ("Se", 78.971),
            ("Mo", 95.95),
            ("Cd", 112.41),
            ("Dy", 162.50),
            ("Yb", 173.05),
            ("Pb", 207.2),
            ("U", 238.03),
        ],
    )
    def test_masses_standard_abregees_de_la_ciaaw(self, symbole, masse):
        assert element(symbole).masse_molaire == masse

    def test_84_elements_ont_une_masse_standard(self, tableau):
        assert sum(e.masse_molaire is not None for e in tableau) == 84

    def test_sans_masse_standard_un_isotope_de_reference(self, tableau):
        assert all((e.masse_molaire is None) == (e.isotope_reference is not None) for e in tableau)

    @pytest.mark.parametrize(("symbole", "A"), [("Tc", 97), ("Pm", 145), ("Po", 209), ("Pu", 244)])
    def test_masse_d_un_element_radioactif_est_son_nombre_de_masse(self, symbole, A):
        assert element(symbole).masse == A

    @pytest.mark.parametrize(("symbole", "cas"), [("H", "1333-74-0"), ("He", "7440-59-7")])
    def test_numeros_cas(self, symbole, cas):
        assert element(symbole).cas == cas

    def test_numeros_cas_valides_et_uniques(self, tableau):
        # Le dernier chiffre d'un numéro CAS est une somme de contrôle.
        for e in tableau:
            chiffres = e.cas.replace("-", "")
            assert re.fullmatch(r"\d{2,7}-\d{2}-\d", e.cas), e.symbole
            somme = sum(i * int(c) for i, c in enumerate(reversed(chiffres[:-1]), start=1))
            assert somme % 10 == int(chiffres[-1]), e.symbole
        assert len({e.cas for e in tableau}) == 118

    def test_configurations_comptent_z_electrons(self, tableau):
        coeurs = {"He": 2, "Ne": 10, "Ar": 18, "Kr": 36, "Xe": 54, "Rn": 86}
        for e in tableau:
            total = 0
            for terme in e.configuration.split():
                if terme.startswith("["):
                    total += coeurs[terme.strip("[]")]
                else:
                    total += int(re.fullmatch(r"\d[spdf](\d+)", terme).group(1))
            assert total == e.Z, e.symbole

    def test_configuration_du_niobium(self):
        assert element("Nb").configuration == "[Kr] 4d4 5s1"


class TestProprietesPhysiques:
    def test_temperatures_au_dela_de_1000_celsius_conservees(self):
        assert element("Fe").theta_fusion == pytest.approx(1538, abs=1)
        assert element("W").T_ebullition > 5000

    def test_ebullition_jamais_avant_fusion(self, tableau):
        for e in tableau:
            if e.T_fusion is not None and e.T_ebullition is not None:
                assert e.T_fusion < e.T_ebullition, e.symbole

    def test_arsenic_et_carbone_se_subliment(self):
        assert element("As").T_fusion is None
        assert element("As").theta_sublimation == pytest.approx(616, abs=1)
        assert element("C").T_sublimation is not None

    def test_masses_volumiques_en_kg_par_m3(self, tableau):
        for e in tableau:
            if e.masse_volumique is None:
                continue
            if e.etat == "g":
                assert e.masse_volumique < 10, e.symbole
            else:
                assert e.masse_volumique > 500, e.symbole

    @pytest.mark.parametrize(
        ("symbole", "etat"), [("H", "g"), ("Br", "l"), ("Hg", "l"), ("Ga", "s"), ("Rn", "g")]
    )
    def test_etat_a_25_celsius(self, symbole, etat):
        assert element(symbole).etat == etat

    def test_pas_d_etat_sans_echantillon(self):
        assert element("Og").etat is None
        assert element("At").masse_volumique is None


class TestClassement:
    def test_hydrogene_est_un_non_metal(self):
        assert element("H").famille == "non-métal"

    def test_metaux_alcalins(self, tableau):
        alcalins = [e.symbole for e in tableau.filtrer(famille="métal alcalin")]
        assert alcalins == ["Li", "Na", "K", "Rb", "Cs", "Fr"]

    def test_familles_connues(self, tableau):
        assert {e.famille for e in tableau} == set(FAMILLES)

    def test_quinze_lanthanides(self, tableau):
        assert len(tableau.filtrer(famille="lanthanide")) == 15

    def test_blocs(self, tableau):
        assert Counter(e.bloc for e in tableau) == {"s": 14, "p": 36, "d": 40, "f": 28}

    def test_groupe_3_contient_lu_et_lr(self, tableau):
        assert [e.symbole for e in tableau.filtrer(groupe=3)] == ["Sc", "Y", "Lu", "Lr"]

    def test_la_hors_groupe(self):
        assert element("La").groupe is None

    @pytest.mark.parametrize(
        ("symbole", "caractere"), [("Fe", "métal"), ("Si", "métalloïde"), ("Cl", "non-métal")]
    )
    def test_caractere(self, symbole, caractere):
        assert element(symbole).caractere == caractere

    def test_origines(self, tableau):
        assert Counter(e.origine for e in tableau) == {
            "primordial": 83,
            "désintégration": 11,
            "synthétique": 24,
        }

    def test_radioactifs(self, tableau):
        assert sum(e.radioactif for e in tableau) == 38
        assert element("Bi").radioactif


class TestAcces:
    @pytest.mark.parametrize("cle", [26, "Fe", "Fer", "fer", " FER "])
    def test_par_z_symbole_ou_nom(self, cle):
        assert element(cle).Z == 26

    def test_nom_sans_accents(self):
        assert element("helium").symbole == "He"

    def test_symbole_sensible_a_la_casse(self, tableau):
        assert "CO" not in tableau
        assert tableau["Co"].nom == "Cobalt"

    def test_element_inconnu(self, tableau):
        with pytest.raises(KeyError, match="inconnu"):
            tableau["Xy"]

    def test_element_inconnu_montre_l_entree_d_origine(self, tableau):
        with pytest.raises(KeyError, match=r"'FE'.*casse"):
            tableau["FE"]

    def test_cle_sans_accents(self):
        assert element("Hélium").cle == "helium"

    def test_str(self):
        assert str(element("Fe")) == "Fe (fer, Z = 26)"

    def test_valeurs_par_symbole(self, tableau):
        valeurs = tableau.valeurs("masse_molaire")
        assert valeurs["O"] == 15.999
        assert valeurs["Og"] is None

    def test_filtrer_par_predicat(self, tableau):
        tres_electronegatifs = tableau.filtrer(
            lambda e: e.electronegativite is not None and e.electronegativite > 3.4
        )
        assert [e.symbole for e in tres_electronegatifs] == ["O", "F"]

    def test_element_est_gele(self):
        with pytest.raises(AttributeError):
            element("Fe").nom = "Acier"

    def test_tableau_construit_depuis_des_elements(self):
        t = TableauPeriodique([element("O"), element("H")])
        assert [e.symbole for e in t] == ["H", "O"]
        assert len(t) == 2

    def test_lecture_d_un_autre_csv(self, tmp_path):
        source = pychimie.elements.DONNEES.read_text(encoding="utf-8").splitlines()
        (tmp_path / "t.csv").write_text("\n".join(source[:3]) + "\n", encoding="utf-8")
        t = TableauPeriodique.depuis_csv(tmp_path / "t.csv")
        assert [e.symbole for e in t] == ["H", "He"]
        assert isinstance(t[2], Element)


class TestNormaliser:
    def test_extremes_de_l_electronegativite(self, tableau):
        x = tableau.normaliser("electronegativite")
        assert x[9] == 1.0
        assert x[87] == 0.0

    def test_valeur_absente_reste_absente(self, tableau):
        assert tableau.normaliser("electronegativite")[2] is None

    def test_bornes_imposees_et_rabattues(self, tableau):
        x = tableau.normaliser("electronegativite", bornes=(1.0, 3.0))
        assert x[9] == 1.0
        assert x[6] == pytest.approx(0.775)

    def test_attribut_non_numerique(self, tableau):
        with pytest.raises(TypeError):
            tableau.normaliser("famille")

    def test_booleen_refuse(self, tableau):
        with pytest.raises(TypeError):
            tableau.normaliser("radioactif")

    def test_bornes_vides(self, tableau):
        with pytest.raises(ValueError, match="bornes"):
            tableau.normaliser("electronegativite", bornes=(2.0, 2.0))

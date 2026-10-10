import numpy as np
import pytest
from scipy.constants import R
from scipy.optimize import brentq

from pychimie import Reaction, RuptureEquilibre, Systeme

KA_ACETIQUE = 10**-4.76


@pytest.fixture
def nh3_scripts(nh3) -> Systeme:
    """Optimisation_NH3_1_P.py : 1 mol de chaque, 100 bar, 400 K."""
    return Systeme(nh3, {"N2": 1, "H2": 1, "NH3": 1}, T=400, P=100)


@pytest.fixture
def esterification() -> Systeme:
    r = Reaction.depuis_K0("CH3COOH(l1) + C2H5OH(l1) = CH3COOC2H5(l1) + H2O(l1)", 4.0)
    return Systeme(r, {"CH3COOH": 1, "C2H5OH": 1})


@pytest.fixture
def acide_acetique() -> Systeme:
    r = Reaction.depuis_K0("CH3COOH(aq) + H2O(l) = CH3COO-(aq) + H3O+(aq)", KA_ACETIQUE)
    return Systeme(r, {"CH3COOH": 0.1, "H2O": 55}, V=1.0)


def q_nh3(xi, n_N2, n_H2, n_NH3, P):
    """Q écrit à la main, comme dans Optimisation_NH3_v2.py."""
    n_tot = n_N2 + n_H2 + n_NH3 - 2 * xi
    return (n_NH3 + 2 * xi) ** 2 * n_tot**2 / ((n_N2 - xi) * (n_H2 - 3 * xi) ** 3 * P**2)


class TestBornes:
    def test_bornes_du_reactif_et_du_produit_limitants(self, nh3_scripts):
        assert (nh3_scripts.xi_min, nh3_scripts.xi_max) == (-0.5, 1 / 3)

    def test_bornes_nulles_sans_produit(self, nh3):
        s = Systeme(nh3, {"N2": 1, "H2": 3})
        assert (s.xi_min, s.xi_max) == (0.0, 1.0)

    def test_avancements_comme_liste_avancement(self, nh3):
        """Un pas retiré à chaque bout, comme l'ancien ``liste_avancement``."""
        s = Systeme(nh3, {"N2": 1, "H2": 0.01, "NH3": 1})
        xi = s.avancements(debut=0)
        assert xi[[0, -1]] == pytest.approx([0.01 / 3 / 2000, 0.01 / 3 * 1999 / 2000])


class TestQuotient:
    def test_Q_du_melange_gazeux(self, nh3_scripts):
        xi = np.linspace(-0.49, 0.33, 50)
        assert nh3_scripts.Q(xi) == pytest.approx(q_nh3(xi, 1, 1, 1, 100), rel=1e-12)

    def test_Q_garde_la_forme_du_tableau(self, nh3_scripts):
        assert nh3_scripts.Q(np.zeros((2, 3))).shape == (2, 3)

    def test_Q_d_un_scalaire_est_un_flottant(self, nh3_scripts):
        assert isinstance(nh3_scripts.Q(0.1), float)

    def test_Q_sans_temperature(self, nh3):
        assert Systeme(nh3, {"N2": 1, "H2": 1, "NH3": 1}).Q(0.0) == pytest.approx(9.0)

    def test_Q_vaut_zero_et_infini_aux_bornes(self, nh3_scripts):
        assert nh3_scripts.Q([-0.5, 1 / 3]) == pytest.approx([0.0, np.inf])

    def test_solide_pur_activite_un(self, boudouard):
        s = Systeme(boudouard, {"C": 1, "CO2": 1, "CO": 1}, P=2)
        assert s.Q([0.0, 0.5]) == pytest.approx([2 * 0.5, 2 * 2**2 / (0.5 * 2.5)])

    def test_solide_pur_garde_activite_un_quand_il_disparait(self, boudouard):
        s = Systeme(boudouard, {"C": 0.5, "CO2": 1}, P=1)
        assert s.activite("C", 0.5) == 1.0

    def test_soluté_activite_concentration(self, acide_acetique):
        assert acide_acetique.Q(0.01) == pytest.approx(0.01**2 / 0.09)

    def test_melange_liquide_fractions_molaires(self, esterification):
        assert esterification.Q(0.5) == pytest.approx(1.0)

    def test_gaz_a_volume_fixe(self, calcaire):
        s = Systeme(calcaire, {"CaCO3": 1}, T=1100, V_gaz=10)
        assert s.Q(0.1) == pytest.approx(0.1 * R * 1100 / 10e-3 / 1e5)

    def test_inerte_dilue_le_gaz(self, nh3):
        s = Systeme(nh3, {"N2": 1, "H2": 1, "NH3": 1}, inertes={"Ar(g)": 3}, P=1)
        assert s.Q(0.0) == pytest.approx(36.0)


class TestEquilibre:
    def test_Q_egal_K0(self, nh3_scripts):
        assert nh3_scripts.Q(nh3_scripts.xi_eq()) == pytest.approx(nh3_scripts.K0(), rel=1e-13)

    def test_resolution_independante_nh3(self, nh3_scripts):
        attendu = brentq(
            lambda x: np.log(q_nh3(x, 1, 1, 1, 100) / nh3_scripts.K0()),
            -0.5 + 1e-12,
            1 / 3 - 1e-12,
            xtol=1e-15,
        )
        assert nh3_scripts.xi_eq() == pytest.approx(attendu, rel=1e-12)

    @pytest.mark.parametrize("P", [0.1, 1.0, 10.0])
    @pytest.mark.parametrize("T", [250.0, 298.15, 350.0])
    def test_dimerisation_forme_close(self, dimerisation, T, P):
        """N₂O₄ = 2 NO₂ depuis 1 mol de N₂O₄ : 4ξ²P = K°(1 − ξ²)."""
        K = dimerisation.K0(T)
        s = Systeme(dimerisation, {"N2O4": 1}, T=T, P=P)
        assert s.xi_eq() == pytest.approx(np.sqrt(K / (K + 4 * P)), rel=1e-13)

    def test_esterification_deux_tiers(self, esterification):
        assert esterification.xi_eq() == pytest.approx(2 / 3, rel=1e-14)

    def test_acide_faible_trinome(self, acide_acetique):
        """La concentration de l'eau n'entre pas : solvant, activité 1."""
        attendu = (-KA_ACETIQUE + np.sqrt(KA_ACETIQUE**2 + 0.4 * KA_ACETIQUE)) / 2
        assert acide_acetique.xi_eq() == pytest.approx(attendu, rel=1e-12)

    def test_carbonate_a_volume_fixe_forme_close(self, calcaire):
        s = Systeme(calcaire, {"CaCO3": 1}, T=1100, V_gaz=10)
        attendu = calcaire.K0(1100) * 1e5 * 10e-3 / (R * 1100)
        assert s.xi_eq() == pytest.approx(attendu, rel=1e-13)

    def test_K0_enorme_quantite_residuelle_exacte(self):
        """2 SO₂ + O₂ = 2 SO₃ à 298 K, K° ≈ 10²⁴ : il reste 2u mol de SO₂ avec
        u = (1 − u)·√((3 + u)/(K°P(1 + u))), que ξ = 1 − u ne saurait représenter."""
        r = Reaction.depuis_tables("2 SO2(g) + O2(g) = 2 SO3(g)")
        s = Systeme(r, {"SO2": 2, "O2": 2}, T=298.15, P=1)
        K = r.K0(298.15)
        u = 0.0
        for _ in range(5):
            u = (1 - u) * np.sqrt((3 + u) / (K * (1 + u)))
        assert s.etat_final().n("SO2") == pytest.approx(2 * u, rel=1e-10)

    def test_K0_minuscule_avancement_exact(self, dimerisation):
        """À 100 K, K° ≈ 10⁻²⁰ : ξ ≈ √K°/2, loin sous la précision de ξ_max = 1."""
        s = Systeme(dimerisation, {"N2O4": 1}, T=100, P=1)
        K = dimerisation.K0(100)
        assert s.xi_eq() == pytest.approx(np.sqrt(K / (K + 4)), rel=1e-12)

    def test_DrG_nul_a_l_equilibre(self, nh3_scripts):
        assert nh3_scripts.DrG(nh3_scripts.xi_eq()) == pytest.approx(0.0, abs=1e-8)

    def test_DrG_negatif_avant_positif_apres(self, nh3_scripts):
        xi = nh3_scripts.xi_eq()
        assert np.sign(nh3_scripts.DrG([xi - 0.01, xi + 0.001])).tolist() == [-1.0, 1.0]

    def test_etat_final_quantites_et_fractions(self, dimerisation):
        etat = Systeme(dimerisation, {"N2O4": 1}, T=298.15, P=1).etat_final()
        xi = etat.xi
        assert (etat.n("NO2"), etat.fraction("NO2")) == pytest.approx((2 * xi, 2 * xi / (1 + xi)))

    def test_etat_final_taux_et_sens(self, nh3_scripts):
        etat = nh3_scripts.etat_final()
        assert (etat.taux_avancement, etat.sens) == (pytest.approx(etat.xi * 3), 1)

    def test_sens_inverse(self, nh3):
        etat = Systeme(nh3, {"N2": 1, "H2": 3, "NH3": 2}, T=900, P=1).etat_final()
        assert etat.sens == -1

    def test_fraction_d_un_solute_refusee(self, acide_acetique):
        with pytest.raises(ValueError, match="soluté"):
            acide_acetique.etat_final().fraction("H3O+")


class TestEnthalpieLibre:
    @pytest.mark.parametrize(
        "systeme",
        ["nh3_scripts", "acide_acetique_T", "esterification_T", "carbonate_isochore", "inerte"],
    )
    def test_derivee_de_G_est_DrG(self, request, systeme):
        s = request.getfixturevalue(systeme)
        xi = s.xi_min + (s.xi_max - s.xi_min) * np.array([0.2, 0.5, 0.8])
        h = 1e-6 * (s.xi_max - s.xi_min)
        derivee = (s.G(xi + h) - s.G(xi - h)) / (2 * h)
        assert derivee == pytest.approx(s.DrG(xi), rel=1e-6, abs=1e-3)

    def test_G_nul_a_l_etat_initial(self, nh3_scripts):
        assert nh3_scripts.G(0.0) == 0.0

    def test_minimum_de_G_est_l_equilibre(self, nh3_scripts):
        assert nh3_scripts.minimum_G()[0] == pytest.approx(nh3_scripts.xi_eq(), rel=1e-7)

    def test_minimum_de_G_sur_la_borne_en_cas_de_rupture(self, boudouard):
        s = Systeme(boudouard, {"C": 0.1, "CO2": 1}, T=1200, P=2)
        assert s.minimum_G()[0] == pytest.approx(0.1, rel=1e-7)

    def test_G_sans_temperature_refuse(self, nh3):
        with pytest.raises(ValueError, match="température"):
            Systeme(nh3, {"N2": 1, "H2": 3}).G(0.1)


@pytest.fixture
def acide_acetique_T() -> Systeme:
    r = Reaction.depuis_K0("CH3COOH(aq) + H2O(l) = CH3COO-(aq) + H3O+(aq)", KA_ACETIQUE)
    return Systeme(r, {"CH3COOH": 0.1, "H2O": 27, "CH3COO-": 0.05, "H3O+": 0.01}, V=0.5, T=298.15)


@pytest.fixture
def esterification_T() -> Systeme:
    r = Reaction.depuis_K0("CH3COOH(l1) + C2H5OH(l1) = CH3COOC2H5(l1) + H2O(l1)", 4.0)
    return Systeme(r, {"CH3COOH": 1, "C2H5OH": 2, "H2O": 0.5}, T=298.15)


@pytest.fixture
def carbonate_isochore(calcaire) -> Systeme:
    return Systeme(calcaire, {"CaCO3": 1, "CaO": 0.2, "CO2": 0.01}, T=1100, V_gaz=10)


@pytest.fixture
def inerte(nh3) -> Systeme:
    return Systeme(nh3, {"N2": 1, "H2": 2, "NH3": 0.5}, inertes={"Ar(g)": 1}, T=700, P=50)


class TestRupture:
    def test_solide_limitant_disparait(self, boudouard):
        etat = Systeme(boudouard, {"C": 0.1, "CO2": 1}, T=1200, P=2).etat_final()
        assert (etat.xi, etat.equilibre, etat.disparus) == (0.1, False, ("C(s)",))

    def test_rupture_Q_reste_sous_K0(self, boudouard):
        etat = Systeme(boudouard, {"C": 0.1, "CO2": 1}, T=1200, P=2).etat_final()
        assert etat.Q == pytest.approx(0.2**2 / (0.9 * 1.1) * 2)

    def test_quantites_nulles_exactement(self, boudouard):
        etat = Systeme(boudouard, {"C": 0.1, "CO2": 1}, T=1200, P=2).etat_final()
        assert etat.n("C") == 0.0

    def test_xi_eq_leve_rupture(self, boudouard):
        with pytest.raises(RuptureEquilibre, match="C\\(s\\) disparaît"):
            Systeme(boudouard, {"C": 0.1, "CO2": 1}, T=1200, P=2).xi_eq()

    def test_rupture_est_une_valeur_erreur(self):
        assert issubclass(RuptureEquilibre, ValueError)

    def test_rupture_vers_xi_min(self, calcaire):
        """CaO limitant dans le sens inverse : il disparaît, CaCO₃ se forme."""
        etat = Systeme(calcaire, {"CaO": 0.1, "CO2": 1}, T=1000, P=1).etat_final()
        assert (etat.xi, etat.disparus) == (-0.1, ("CaO(s)",))

    def test_carbonate_decompose_totalement_a_pression_fixee(self, calcaire):
        """Seul gaz, CO₂ a l'activité P/P° quelle que soit sa quantité : K°(1200 K) > 1."""
        etat = Systeme(calcaire, {"CaCO3": 1}, T=1200, P=1).etat_final()
        assert (etat.xi, etat.disparus) == (1.0, ("CaCO3(s)",))

    def test_carbonate_stable_sous_K0_inferieur_a_P(self, calcaire):
        etat = Systeme(calcaire, {"CaCO3": 1}, T=1000, P=1).etat_final()
        assert (etat.xi, etat.equilibre) == (0.0, False)

    def test_carbonate_a_volume_fixe_peu_de_solide(self, calcaire):
        etat = Systeme(calcaire, {"CaCO3": 0.05}, T=1100, V_gaz=10).etat_final()
        assert (etat.xi, etat.disparus) == (0.05, ("CaCO3(s)",))

    def test_gaz_qui_disparait_ne_rompt_pas(self, nh3):
        etat = Systeme(nh3, {"N2": 0.1, "H2": 3}, T=300, P=100).etat_final()
        assert etat.equilibre and 0 < etat.n("N2") < 1e-6

    def test_aucune_evolution_possible(self, calcaire):
        with pytest.raises(RuptureEquilibre, match="aucune évolution"):
            Systeme(calcaire, {"CaO": 1}, T=1000).xi_eq()


class TestConstruction:
    def test_constituant_absent_vaut_zero(self, nh3):
        assert Systeme(nh3, {"N2": 1}).initial == (1.0, 0.0, 0.0)

    def test_constituant_inconnu(self, nh3):
        with pytest.raises(KeyError):
            Systeme(nh3, {"O2": 1})

    def test_quantite_negative(self, nh3):
        with pytest.raises(ValueError, match="négative"):
            Systeme(nh3, {"N2": -1})

    def test_donne_deux_fois(self, nh3):
        with pytest.raises(ValueError, match="deux fois"):
            Systeme(nh3, {"N2": 1, "N2(g)": 2})

    def test_inerte_qui_participe(self, nh3):
        with pytest.raises(ValueError, match="participe"):
            Systeme(nh3, {"H2": 1}, inertes={"N2(g)": 1})

    def test_inerte_sans_phase(self, nh3):
        with pytest.raises(ValueError, match="phase"):
            Systeme(nh3, {"H2": 1}, inertes={"Ar": 1})

    def test_solute_sans_volume(self):
        r = Reaction.depuis_K0("CH3COOH(aq) + H2O(l) = CH3COO-(aq) + H3O+(aq)", KA_ACETIQUE)
        with pytest.raises(ValueError, match="volume V"):
            Systeme(r, {"CH3COOH": 0.1})

    def test_volume_de_gaz_sans_temperature(self, nh3):
        with pytest.raises(ValueError, match="demande T"):
            Systeme(nh3, {"N2": 1}, V_gaz=1)

    @pytest.mark.parametrize("cle", ["T", "P", "V", "V_gaz"])
    def test_condition_non_positive(self, nh3, cle):
        with pytest.raises(ValueError, match="strictement positif"):
            Systeme(nh3, {"N2": 1}, **{"T": 300, cle: 0.0})

    def test_avec_remplace_la_composition(self, nh3):
        s = Systeme(nh3, {"N2": 1, "H2": 3, "NH3": 1}, T=500).avec(initial={"H2": 2}, P=10)
        assert (s.initial, s.T, s.P) == ((0.0, 2.0, 0.0), 500.0, 10.0)

    def test_avec_parametre_inconnu(self, nh3):
        with pytest.raises(TypeError):
            Systeme(nh3, {"N2": 1}).avec(Temperature=300)

    def test_pression_a_volume_fixe(self, calcaire):
        s = Systeme(calcaire, {"CaCO3": 1}, T=1100, V_gaz=10)
        assert s.pression(0.2) == pytest.approx(0.2 * R * 1100 / 10e-3 / 1e5)


class TestGardesUnites:
    def test_pression_en_pascals_refusee(self, nh3):
        with pytest.raises(ValueError, match=r"bar.*2e\+07 Pa = 200 bar"):
            Systeme(nh3, {"N2": 1}, T=700, P=200e5)

    def test_pression_limite_acceptee(self, nh3):
        assert Systeme(nh3, {"N2": 1}, T=700, P=1e4).P == 1e4

    def test_temperature_en_celsius_refusee(self, nh3):
        with pytest.raises(ValueError, match=r"kelvins.*25 °C = 298.15 K"):
            Systeme(nh3, {"N2": 1}, T=25)

    def test_temperature_limite_acceptee(self, nh3):
        assert Systeme(nh3, {"N2": 1}, T=100).T == 100.0

    def test_avec_passe_par_les_memes_gardes(self, nh3_scripts):
        with pytest.raises(ValueError, match="bar"):
            nh3_scripts.avec(P=1e5)


class TestK0DrG0:
    """``Systeme.K0(T=None)`` et ``DrG0(T=None)`` : la signature de ``Reaction``."""

    def test_K0_a_la_temperature_du_systeme(self, nh3_scripts, nh3):
        assert nh3_scripts.K0() == pytest.approx(nh3.K0(400), rel=1e-15)

    def test_K0_a_une_autre_temperature(self, nh3_scripts, nh3):
        assert nh3_scripts.K0(800) == pytest.approx(nh3.K0(800), rel=1e-15)

    def test_K0_accepte_un_tableau(self, nh3_scripts, nh3):
        T = np.array([400.0, 600.0, 800.0])
        assert nh3_scripts.K0(T) == pytest.approx(nh3.K0(T), rel=1e-15)

    def test_DrG0_a_la_temperature_du_systeme(self, nh3_scripts, nh3):
        assert nh3_scripts.DrG0() == pytest.approx(nh3.DrG0(400), rel=1e-15)

    def test_DrG0_a_une_autre_temperature(self, nh3_scripts, nh3):
        assert nh3_scripts.DrG0(800) == pytest.approx(nh3.DrG0(800), rel=1e-15)

    def test_sans_temperature_K0_demande_T(self, nh3):
        with pytest.raises(ValueError, match="donner T"):
            Systeme(nh3, {"N2": 1}).K0()

    def test_sans_temperature_DrG0_demande_T(self, nh3):
        with pytest.raises(ValueError, match="demande la température"):
            Systeme(nh3, {"N2": 1}).DrG0()

    def test_K0_independante_de_T_se_passe_de_T(self, esterification):
        assert esterification.K0() == pytest.approx(4.0)

    def test_DrG_coherent_avec_DrG0(self, nh3_scripts):
        xi = 0.1
        attendu = nh3_scripts.DrG0() + R * 400 * np.log(nh3_scripts.Q(xi))
        assert nh3_scripts.DrG(xi) == pytest.approx(attendu, rel=1e-12)

import numpy as np
import pytest
from scipy.constants import R
from scipy.integrate import quad

from pychimie import Equation, Reaction, table_thermo


class TestHess:
    def test_synthese_ammoniac(self, nh3_tables):
        # 2 × (−46,11) kJ/mol ; 2 × 192,45 − 191,61 − 3 × 130,684 J/K/mol
        assert (nh3_tables.DrH0(), nh3_tables.DrS0()) == pytest.approx((-92.22e3, -198.762))

    def test_DrCp_avec_kirchhoff(self):
        r = Reaction.depuis_tables("N2(g) + 3 H2(g) = 2 NH3(g)", kirchhoff=True)
        assert r.DrCp0 == pytest.approx(2 * 35.06 - 29.125 - 3 * 28.824)

    def test_sans_kirchhoff_DrCp_nul(self, nh3_tables):
        assert nh3_tables.DrCp0 == 0.0

    def test_boudouard_redonne_les_valeurs_du_cours(self):
        r = Reaction.depuis_tables("C(s) + CO2(g) = 2 CO(g)")
        assert (r.DrH0() / 1e3, r.DrS0()) == pytest.approx((172.5, 175.8), abs=0.1)

    def test_combinaison_lineaire_des_equations(self):
        """CO + ½O₂ = CO₂ est (C + O₂ = CO₂) − (C + ½O₂ = CO) : Hess est linéaire."""
        combustion_co = Reaction.depuis_tables("CO(g) + 1/2 O2(g) = CO2(g)", kirchhoff=True)
        c_co2 = Reaction.depuis_tables("C(s) + O2(g) = CO2(g)", kirchhoff=True)
        c_co = Reaction.depuis_tables("C(s) + 1/2 O2(g) = CO(g)", kirchhoff=True)
        for grandeur in ("DrH0_ref", "DrS0_ref", "DrCp0"):
            attendu = getattr(c_co2, grandeur) - getattr(c_co, grandeur)
            assert getattr(combustion_co, grandeur) == pytest.approx(attendu)

    def test_phase_etiquetee_prend_le_corps_pur(self):
        pur = Reaction.depuis_tables("CO(g) + 2 H2(g) = CH3OH(l)")
        etiquete = Reaction.depuis_tables("CO(g) + 2 H2(g) = CH3OH(l1)")
        assert etiquete.DrH0() == pur.DrH0()

    def test_espece_absente_des_tables(self):
        with pytest.raises(KeyError):
            Reaction.depuis_tables("NH4Cl(s) = NH3(g) + HCl(g)")

    def test_autre_table(self):
        table = table_thermo()
        assert Reaction.depuis_tables("N2O4(g) = 2 NO2(g)", table=table).DrH0() == pytest.approx(
            57.2e3
        )


class TestEllingham:
    def test_K0_exp_moins_DrG0_sur_RT(self, boudouard):
        attendu = np.exp(-(172.5e3 - 1200 * 175.8) / (R * 1200))
        assert boudouard.K0(1200) == pytest.approx(attendu, rel=1e-14)

    def test_K0_vaut_1_a_la_temperature_d_inversion(self, boudouard):
        assert boudouard.K0(172.5e3 / 175.8) == pytest.approx(1.0, rel=1e-12)

    def test_DrH0_independant_de_T(self, boudouard):
        assert boudouard.DrH0(1500) == boudouard.DrH0()

    def test_K0_se_calcule_sur_un_tableau(self, boudouard):
        T = np.linspace(500, 1500, 11)
        assert boudouard.K0(T) == pytest.approx([boudouard.K0(t) for t in T], rel=1e-14)

    def test_K0_scalaire_rend_un_flottant(self, boudouard):
        assert isinstance(boudouard.K0(1000.0), float)

    def test_ln_K0_reste_fini_quand_K0_deborde(self):
        r = Reaction("A(g) = B(g)", DrH0=-3e6, DrS0=0.0)
        assert np.isfinite(r.ln_K0(298.15)) and r.ln_K0(298.15) > 1000


class TestKirchhoff:
    @pytest.fixture
    def nh3_cp(self):
        return Reaction.depuis_tables("N2(g) + 3 H2(g) = 2 NH3(g)", kirchhoff=True)

    @pytest.mark.parametrize("T", [400.0, 700.0, 1000.0])
    def test_ln_K0_suit_van_t_hoff_integree(self, nh3_cp, T):
        """d ln K°/dT = ΔrH°(T)/RT², intégrée numériquement depuis 298,15 K."""
        integrale, _ = quad(lambda t: nh3_cp.DrH0(t) / (R * t**2), 298.15, T)
        assert nh3_cp.ln_K0(T) == pytest.approx(nh3_cp.ln_K0(298.15) + integrale, rel=1e-10)

    def test_DrH0_affine_en_T(self, nh3_cp):
        assert nh3_cp.DrH0(398.15) - nh3_cp.DrH0(298.15) == pytest.approx(100 * nh3_cp.DrCp0)

    def test_DrS0_en_logarithme(self, nh3_cp):
        assert nh3_cp.DrS0(596.3) - nh3_cp.DrS0() == pytest.approx(nh3_cp.DrCp0 * np.log(2))

    def test_egal_a_ellingham_a_la_reference(self, nh3_cp):
        assert nh3_cp.K0(298.15) == pytest.approx(nh3_cp.ellingham().K0(298.15), rel=1e-14)

    def test_ellingham_annule_DrCp(self, nh3_cp):
        assert (nh3_cp.kirchhoff, nh3_cp.ellingham().kirchhoff) == (True, False)

    def test_ecart_a_ellingham(self, nh3_cp):
        """ln(K°/K°_Ellingham) = ΔrC°p/R · (ln(T/T₀) + T₀/T − 1) : à 800 K, K° est divisé
        par 7 pour NH₃ — l'approximation d'Ellingham n'est pas anodine."""
        T, T0 = 800.0, 298.15
        ecart = nh3_cp.ln_K0(T) - nh3_cp.ellingham().ln_K0(T)
        assert ecart == pytest.approx(nh3_cp.DrCp0 / R * (np.log(T / T0) + T0 / T - 1), rel=1e-12)
        assert np.exp(ecart) == pytest.approx(1 / 7, rel=0.03)


class TestDepuisK0:
    def test_sans_DrH0_K0_independant_de_T(self):
        r = Reaction.depuis_K0("A(l1) + B(l1) = C(l1) + D(l1)", 4.0)
        assert (r.K0(), r.K0(350.0), r.K0(800.0)) == pytest.approx((4.0, 4.0, 4.0))

    def test_avec_DrH0_loi_de_van_t_hoff(self):
        r = Reaction.depuis_K0("N2O4(g) = 2 NO2(g)", 0.15, T=298.15, DrH0=57.2e3)
        attendu = 0.15 * np.exp(-57.2e3 / R * (1 / 350 - 1 / 298.15))
        assert r.K0(350) == pytest.approx(attendu, rel=1e-12)

    def test_K0_donnee_retrouvee_a_sa_temperature(self):
        assert Reaction.depuis_K0("A(g) = B(g)", 2.5e-3, T=500, DrH0=-1e4).K0(500) == pytest.approx(
            2.5e-3, rel=1e-14
        )

    def test_K0_negative_refusee(self):
        with pytest.raises(ValueError):
            Reaction.depuis_K0("A(g) = B(g)", -1.0)


class TestErreurs:
    def test_temperature_en_celsius_refusee(self, nh3):
        with pytest.raises(ValueError, match="kelvins"):
            nh3.K0(25)

    def test_tableau_de_temperatures_avec_une_valeur_basse(self, nh3):
        with pytest.raises(ValueError, match="kelvins"):
            nh3.K0([300, 50, 400])

    def test_pas_de_DrH0_DrS0_par_defaut(self):
        """L'ancien Reaction prenait ΔrH° = ΔrS° = 1 en silence."""
        with pytest.raises(TypeError, match="DrS0"):
            Reaction("N2(g) + 3 H2(g) = 2 NH3(g)", DrH0=-92e3)

    def test_K0_sans_T_si_DrH0_non_nul(self, nh3):
        with pytest.raises(ValueError, match="donner T"):
            nh3.K0()

    @pytest.mark.parametrize("T", [0.0, -10.0])
    def test_temperature_non_positive(self, nh3, T):
        with pytest.raises(ValueError, match="kelvins"):
            nh3.K0(T)


def test_equation_lue_ou_donnee(nh3):
    assert Reaction(Equation.lire("N2(g) + 3 H2(g) = 2 NH3(g)"), DrH0=-92e3, DrS0=-175.6) == nh3


def test_str_donne_les_grandeurs(nh3_tables):
    assert str(nh3_tables).startswith("N2(g) + 3 H2(g) = 2 NH3(g) : ΔrH° = -92.22 kJ/mol")

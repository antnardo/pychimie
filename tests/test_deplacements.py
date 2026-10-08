import numpy as np
import pytest
from scipy.optimize import minimize_scalar

from pychimie import Reaction, Systeme, balayage, perturber


@pytest.fixture
def nh3_equilibre(nh3_tables) -> Systeme:
    return Systeme(nh3_tables, {"N2": 1, "H2": 3}, T=700, P=200)


class TestLoisDeModeration:
    def test_hausse_de_T_exothermique_sens_inverse(self, nh3_equilibre):
        p = perturber(nh3_equilibre, T=800)
        assert (p.sens, p.deplacement < 0) == (-1, True)

    def test_hausse_de_T_endothermique_sens_direct(self, boudouard):
        p = perturber(Systeme(boudouard, {"C": 2, "CO2": 1}, T=900, P=1), T=1000)
        assert (p.sens, p.deplacement > 0) == (1, True)

    def test_hausse_de_P_vers_moins_de_gaz(self, nh3_equilibre):
        assert perturber(nh3_equilibre, P=300).sens == 1

    def test_hausse_de_P_dimerisation(self, dimerisation):
        assert perturber(Systeme(dimerisation, {"N2O4": 1}, T=298.15), P=10).sens == -1

    def test_P_sans_effet_si_delta_nu_gaz_nul(self):
        r = Reaction.depuis_tables("CO(g) + H2O(g) = CO2(g) + H2(g)")
        p = perturber(Systeme(r, {"CO": 1, "H2O": 1}, T=1000, P=1), P=50)
        assert p.deplacement == pytest.approx(0.0, abs=1e-14)

    def test_inerte_a_pression_fixee_comme_une_detente(self, nh3_equilibre):
        assert perturber(nh3_equilibre, ajout={"Ar(g)": 1}).sens == -1

    def test_inerte_a_volume_fixe_sans_effet(self, nh3_tables):
        s = Systeme(nh3_tables, {"N2": 1, "H2": 3}, T=700, V_gaz=10)
        p = perturber(s, ajout={"Ar(g)": 1})
        assert (p.sens, p.deplacement) == (0, 0.0)

    def test_inerte_deja_present_s_ajoute(self, nh3_tables):
        s = Systeme(nh3_tables, {"N2": 1, "H2": 3}, T=700, P=200, inertes={"Ar(g)": 1})
        assert perturber(s, ajout={"Ar(g)": 1}).systeme.inertes[0][1] == 2.0

    def test_reactif_ajoute_en_defaut_sens_direct(self, nh3_tables):
        s = Systeme(nh3_tables, {"N2": 1, "H2": 3}, T=700, P=200)
        assert perturber(s, ajout={"N2": 0.1}).sens == 1

    def test_N2_ajoute_quand_x_N2_depasse_un_demi_sens_inverse(self, nh3_tables):
        s = Systeme(nh3_tables, {"N2": 3, "H2": 1}, T=700, P=200)
        assert s.etat_final().fraction("N2") > 0.5
        assert perturber(s, ajout={"N2": 0.1}).sens == -1

    @pytest.mark.parametrize(
        "changement",
        [
            {"T": 650},
            {"T": 750},
            {"P": 100},
            {"P": 400},
            {"ajout": {"H2": 0.5}},
            {"ajout": {"NH3": 0.2}},
        ],
    )
    def test_sens_est_le_signe_du_deplacement(self, nh3_equilibre, changement):
        p = perturber(nh3_equilibre, **changement)
        assert p.sens == np.sign(p.deplacement)


class TestPerturbation:
    def test_part_de_la_composition_d_equilibre(self, nh3_equilibre):
        p = perturber(nh3_equilibre, T=800)
        assert p.systeme.initial == pytest.approx(p.avant.quantites)

    def test_Q_juste_apres_est_l_ancien_K0_si_seul_T_change(self, nh3_equilibre):
        p = perturber(nh3_equilibre, T=800)
        assert (p.Q, p.K0) == pytest.approx((nh3_equilibre.K0, nh3_equilibre.avec(T=800).K0))

    def test_nouvel_equilibre_independant_du_chemin(self, nh3_equilibre):
        apres = perturber(nh3_equilibre, T=800).apres
        direct = nh3_equilibre.avec(T=800).etat_final()
        assert apres.quantites == pytest.approx(direct.quantites, rel=1e-12)

    def test_ajout_d_un_nom_inconnu_sans_phase(self, nh3_equilibre):
        with pytest.raises(KeyError):
            perturber(nh3_equilibre, ajout={"Ar": 1})

    def test_apres_une_rupture(self, boudouard):
        """Le carbone a disparu : refroidir fait reformer du carbone (sens inverse)."""
        s = Systeme(boudouard, {"C": 0.1, "CO2": 1}, T=1200, P=1)
        p = perturber(s, T=800)
        assert (p.avant.equilibre, p.sens, p.apres.equilibre) == (False, -1, True)


class TestBalayage:
    def test_proportions_stoechiometriques_maximisent_x_NH3(self, nh3):
        """Le résultat de Optimisation_NH3_3_r.py : x(NH₃) est maximal pour N₂:H₂ = 1:3."""
        s = Systeme(nh3, {"N2": 1, "H2": 1}, T=450, P=50)

        def moins_x_nh3(x):
            return -s.avec(initial={"N2": x, "H2": 1 - x}).etat_final().fraction("NH3")

        res = minimize_scalar(
            moins_x_nh3, bounds=(0.05, 0.95), method="bounded", options={"xatol": 1e-9}
        )
        assert res.x == pytest.approx(0.25, abs=1e-6)

    def test_taux_d_avancement_decroit_avec_T(self, nh3_equilibre):
        b = balayage(lambda T: nh3_equilibre.avec(T=T), np.linspace(500, 900, 9))
        assert np.all(np.diff(b.taux_avancement) < 0)

    def test_taux_d_avancement_croit_avec_P(self, nh3_equilibre):
        b = balayage(lambda P: nh3_equilibre.avec(P=P), [1, 10, 100, 1000])
        assert np.all(np.diff(b.xi) > 0)

    def test_marque_les_ruptures(self, boudouard):
        s = Systeme(boudouard, {"C": 0.1, "CO2": 1}, P=1)
        b = balayage(lambda T: s.avec(T=T), [700, 1300])
        assert b.equilibre.tolist() == [True, False]

    def test_tableaux_alignes_sur_les_valeurs(self, nh3_equilibre):
        b = balayage(lambda n: nh3_equilibre.avec(inertes={"Ar(g)": n}), [0, 1, 2])
        assert [len(b.valeurs), len(b.xi), len(b.n("Ar(g)")), len(b.fraction("NH3"))] == [3] * 4

    def test_inerte_croissant_diminue_xi(self, nh3_equilibre):
        b = balayage(lambda n: nh3_equilibre.avec(inertes={"Ar(g)": n}), [0, 1, 5, 20])
        assert np.all(np.diff(b.xi) < 0)

    def test_xi_max(self, nh3_equilibre):
        b = balayage(lambda x: nh3_equilibre.avec(initial={"N2": x, "H2": 1 - x}), [0.1, 0.5])
        assert b.xi_max == pytest.approx([0.1, 1 / 6])

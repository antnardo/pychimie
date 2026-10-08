"""L'ancien ``chimie/thermochimie.py`` est remplacé par ``pychimie`` : rien ne s'est perdu.

Les valeurs attendues ont été calculées le 8 octobre 2026 avec l'ancien module, avant
sa suppression, en remplaçant le ``R_GP`` en ``Decimal`` de ``constantes`` (qui le
faisait planter dès ``Ko``) par le même R en flottant. Chaque test refait le calcul
avec ``pychimie`` et nomme l'API qui remplace l'ancienne :

- ``Reaction(dic, DrHo, DrSo)`` et ``set_CI`` : ``Reaction(équation, DrH0=, DrS0=)``
  et ``Systeme(reaction, initial, T=, P=)`` ;
- ``xmin()``, ``xmax()`` : ``Systeme.xi_min``, ``Systeme.xi_max`` ;
- ``Q(x, P)``, ``DrG(x, P, T)``, ``Ko(T)`` : ``Systeme.Q``, ``Systeme.DrG``,
  ``Reaction.K0`` (ou ``Systeme.K0``) ;
- ``xeq(P, T)`` : ``Systeme.xi_eq()`` ; ``liste_avancement`` : ``Systeme.avancements`` ;
- ``compute_G`` (cassé : ``quad`` rend un tuple) : ``Systeme.G`` ;
- ``is_rupture`` (ne testait que le réactif limitant) : ``Systeme.etat_final()``.

L'ancien ``xeq`` s'arrêtait à la tolérance par défaut de ``bisect`` (2·10⁻¹²) : l'accord
est à cette précision-là.
"""

import numpy as np
import pytest

from pychimie import Reaction, RuptureEquilibre, Systeme


@pytest.fixture
def boudouard_script(boudouard) -> Systeme:
    return Systeme(boudouard, {"C": 2, "CO2": 1, "CO": 1}, T=1200, P=2)


@pytest.fixture
def nh3_script(nh3) -> Systeme:
    return Systeme(nh3, {"N2": 1, "H2": 1, "NH3": 1}, T=400, P=100)


def test_bornes_reproduisent_xmin_xmax(boudouard_script, nh3_script):
    bornes = [
        boudouard_script.xi_min,
        boudouard_script.xi_max,
        nh3_script.xi_min,
        nh3_script.xi_max,
    ]
    assert bornes == [-0.5, 1.0, -0.5, 0.3333333333333333]


def test_K0_reproduit_Ko(boudouard, nh3):
    assert (boudouard.K0(1200), nh3.K0(400)) == pytest.approx(
        (47.21581873428936, 694.2388838533769), rel=1e-14
    )


def test_xi_eq_reproduit_xeq(boudouard_script, nh3_script):
    attendu = (0.887085600769502, 0.3284046982824331)
    assert (boudouard_script.xi_eq(), nh3_script.xi_eq()) == pytest.approx(attendu, abs=2e-12)


def test_Q_reproduit_Q(boudouard_script, nh3_script):
    attendu_b = [0.035714285714285705, 1.0, 6.400000000000001, 54.06896551724139]
    attendu_n = [
        3.874637758935277e-06,
        0.0009000000000000003,
        0.02587812500000001,
        2252.0243820895485,
    ]
    assert boudouard_script.Q([-0.4, 0.0, 0.5, 0.9]) == pytest.approx(attendu_b, rel=1e-13)
    assert nh3_script.Q([-0.4, 0.0, 0.2, 0.33]) == pytest.approx(attendu_n, rel=1e-13)


def test_DrG_reproduit_DrG(boudouard_script):
    attendu = [-71706.58780307203, -38460.0, -19939.05570114242, 1352.244828893854]
    assert boudouard_script.DrG([-0.4, 0.0, 0.5, 0.9]) == pytest.approx(attendu, rel=1e-12)


@pytest.mark.parametrize(
    ("x", "attendu"),
    [
        (0.1, 0.09999852369048387),
        (0.25, 0.23698655015995396),
        (0.5, 0.1626281466223693),
        (0.9, 0.03209093750145562),
    ],
)
def test_proportions_reproduisent_Optimisation_NH3_3_r(nh3, x, attendu):
    s = Systeme(nh3, {"N2": x, "H2": 1 - x}, T=450, P=50)
    assert s.xi_eq() == pytest.approx(attendu, abs=2e-12)


def test_avancements_et_Q_reproduisent_Optimisation_NH3_2_ni(nh3):
    s = Systeme(nh3, {"N2": 1, "H2": 0.01, "NH3": 1})
    xi = s.avancements(debut=0)
    assert xi[[0, 1, -1]] == pytest.approx(
        [1.6666666666666669e-06, 3.332499583124896e-06, 0.003331666666666667], rel=1e-14
    )
    assert s.Q([1e-4, 1e-3]) == pytest.approx([4428000.008319519, 11814170.976387467], rel=1e-13)


def test_rupture_que_l_ancien_xeq_ne_savait_pas_traiter(boudouard):
    """L'ancien ``xeq`` levait ``ValueError`` (« bisect ne fonctionne pas ») ; l'état
    final est la disparition du carbone, et ``xi_eq`` le dit."""
    s = Systeme(boudouard, {"C": 0.1, "CO2": 1, "CO": 0}, T=1200, P=2)
    with pytest.raises(RuptureEquilibre, match="C\\(s\\)"):
        s.xi_eq()


def test_G_remplace_compute_G(nh3_script):
    """``compute_G`` plantait (``quad`` rend un tuple) ; G se dérive en ΔrG."""
    xi = np.linspace(-0.4, 0.3, 8)
    h = 1e-7
    assert (nh3_script.G(xi + h) - nh3_script.G(xi - h)) / (2 * h) == pytest.approx(
        nh3_script.DrG(xi), rel=1e-6
    )


def test_DrHo_DrSo_ne_valent_plus_1_par_defaut():
    with pytest.raises(TypeError):
        Reaction("C(s) + CO2(g) = 2 CO(g)")

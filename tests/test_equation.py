import pytest

from pychimie import Constituant, Equation


class TestLire:
    def test_coefficients_entiers_signes_par_membre(self):
        eq = Equation.lire("N2(g) + 3 H2(g) = 2 NH3(g)")
        assert eq.nu == (-1.0, -3.0, 2.0)

    def test_noms_portent_la_phase(self):
        eq = Equation.lire("C(s) + CO2(g) = 2 CO(g)")
        assert eq.noms == ("C(s)", "CO2(g)", "CO(g)")

    @pytest.mark.parametrize(
        ("texte", "nu"),
        [
            ("CO(g) + 1/2 O2(g) = CO2(g)", (-1.0, -0.5, 1.0)),
            ("CO(g) + 0.5 O2(g) = CO2(g)", (-1.0, -0.5, 1.0)),
            ("CO(g) + 0,5 O2(g) = CO2(g)", (-1.0, -0.5, 1.0)),
            ("2SO2(g) + O2(g) = 2SO3(g)", (-2.0, -1.0, 2.0)),
        ],
    )
    def test_coefficients_fractionnaires_ou_colles(self, texte, nu):
        assert Equation.lire(texte).nu == nu

    @pytest.mark.parametrize("fleche", ["=", "⇌", "->", "→"])
    def test_separateurs_des_membres(self, fleche):
        assert Equation.lire(f"N2O4(g) {fleche} 2 NO2(g)").nu == (-1.0, 2.0)

    def test_parentheses_dans_la_formule(self):
        eq = Equation.lire("Ca(OH)2(s) = Ca2+(aq) + 2 HO-(aq)")
        assert [(c.formule, c.phase) for c in eq.constituants] == [
            ("Ca(OH)2", "s"),
            ("Ca2+", "aq"),
            ("HO-", "aq"),
        ]

    def test_charges_restent_dans_les_formules(self):
        eq = Equation.lire("CH3COOH(aq) + H2O(l) = CH3COO-(aq) + H3O+(aq)")
        assert eq.noms[-1] == "H3O+(aq)"

    def test_str_relit_la_meme_equation(self):
        eq = Equation.lire("CO(g) + 1/2 O2(g) = CO2(g)")
        assert Equation.lire(str(eq)) == eq

    def test_str_ecrit_les_fractions(self):
        assert str(Equation.lire("CO(g) + 0.5 O2(g) = CO2(g)")) == "CO(g) + 1/2 O2(g) = CO2(g)"

    @pytest.mark.parametrize(
        "texte",
        [
            "N2(g) + 3 H2(g)",
            "N2 + 3 H2(g) = 2 NH3(g)",
            "N2(x) + 3 H2(g) = 2 NH3(g)",
            "N2(g) + N2(g) = 2 N(g)",
            "N2(g) = N2(g) + H2(g)",
        ],
    )
    def test_equation_mal_ecrite_refusee(self, texte):
        with pytest.raises(ValueError):
            Equation.lire(texte)


class TestConstituant:
    @pytest.mark.parametrize(
        ("phase", "etat"), [("g", "g"), ("aq", "aq"), ("aq1", "aq"), ("s", "s"), ("l2", "l")]
    )
    def test_etat_ignore_l_etiquette(self, phase, etat):
        assert Constituant("X", phase).etat == etat

    @pytest.mark.parametrize(
        ("phase", "melange"),
        [("g", "g"), ("s", "Fe(s)"), ("l", "Fe(l)"), ("l1", "l1"), ("aq", None), ("aq2", None)],
    )
    def test_melange_seul_sans_etiquette(self, phase, melange):
        assert Constituant("Fe", phase).melange == melange

    def test_lire_un_inactif(self):
        assert Constituant.lire("Ar(g)") == Constituant("Ar", "g", 0.0)

    @pytest.mark.parametrize("texte", ["Ar", "2 Ar(g)"])
    def test_lire_refuse_sans_phase_ou_avec_coefficient(self, texte):
        with pytest.raises(ValueError):
            Constituant.lire(texte)


class TestIndice:
    @pytest.fixture
    def vaporisation(self):
        return Equation.lire("CH4(g) + 2 O2(g) = CO2(g) + 2 H2O(g)")

    @pytest.mark.parametrize(("nom", "i"), [("CH4", 0), ("O2(g)", 1), (" H2O ", 3)])
    def test_par_formule_ou_par_nom(self, vaporisation, nom, i):
        assert vaporisation.indice(nom) == i

    def test_formule_ambigue_demande_la_phase(self):
        eq = Equation.lire("H2O(l) = H2O(g)")
        with pytest.raises(KeyError, match="ambigu"):
            eq.indice("H2O")

    def test_formule_ambigue_trouvee_avec_la_phase(self):
        assert Equation.lire("H2O(l) = H2O(g)").indice("H2O(g)") == 1

    def test_constituant_absent(self, vaporisation):
        with pytest.raises(KeyError, match="n'est pas un constituant"):
            vaporisation.indice("NH3")


@pytest.mark.parametrize(
    ("texte", "delta"),
    [
        ("N2(g) + 3 H2(g) = 2 NH3(g)", -2.0),
        ("C(s) + CO2(g) = 2 CO(g)", 1.0),
        ("CO(g) + H2O(g) = CO2(g) + H2(g)", 0.0),
    ],
)
def test_delta_nu_gaz_ignore_les_phases_condensees(texte, delta):
    assert Equation.lire(texte).delta_nu_gaz == delta

import importlib.util
import math
from pathlib import Path

import pytest

from pychimie import Constituant, Formule, Reaction, TableThermo, table_thermo

RACINE = Path(__file__).resolve().parents[1]
# Les espèces que l'audit du 10 octobre 2026 relevait absentes : dix équations
# classiques des sujets, trois seulement étaient couvertes.
ESPECES_DES_SUJETS = [
    "NaCl(s)", "Na+(aq)", "Cl-(aq)", "HCl(g)", "NH4Cl(s)", "Al(s)", "Al2O3(s)",
    "PCl5(g)", "PCl3(g)", "NO(g)", "I2(s)", "HI(g)", "C2H4(g)", "C2H5OH(l)",
    "C2H5OH(g)", "H+(aq)", "OH-(aq)", "SO4^2-(aq)", "Fe^2+(aq)", "Fe^3+(aq)",
]  # fmt: skip


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location(
        "ecrire_donnees_thermo", RACINE / "scripts" / "ecrire_donnees_thermo.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestTable:
    def test_contient_les_especes_des_sujets(self):
        assert all(cle in table_thermo() for cle in ESPECES_DES_SUJETS)

    def test_au_moins_quatre_vingts_especes(self):
        assert len(table_thermo()) >= 80

    def test_cles_dans_l_ordre_du_script(self, script):
        assert table_thermo().cles() == tuple(f"{f}({p})" for f, p, *_ in script.ESPECES)

    def test_toutes_les_formules_se_lisent(self):
        assert all(Formule.lire(e.formule).atomes > 0 for e in table_thermo())

    def test_csv_nomme_les_especes_du_script(self, script):
        lu = [(e.formule, e.phase, e.nom, e.source) for e in table_thermo()]
        assert lu == [(f, p, nom, source) for f, p, nom, *_, source in script.ESPECES]

    def test_csv_donne_les_valeurs_du_script(self, script):
        lu = [(e.DfH0 / 1e3, e.Sm0) for e in table_thermo()]
        assert sum(lu, ()) == pytest.approx(
            sum((tuple(ligne[3:5]) for ligne in script.ESPECES), ())
        )
        assert [e.Cp0 for e in table_thermo()] == [ligne[5] for ligne in script.ESPECES]

    def test_enthalpie_en_joules_par_mole(self):
        assert table_thermo()["NH3(g)"].DfH0 == pytest.approx(-46110.0)

    @pytest.mark.parametrize(
        "cle", ["N2(g)", "H2(g)", "O2(g)", "C(s)", "Fe(s)", "Cu(s)", "Zn(s)", "Al(s)", "Na(s)"]
    )
    def test_corps_simples_de_reference_ont_DfH0_nul(self, cle):
        assert table_thermo()[cle].DfH0 == 0.0

    def test_convention_de_l_ion_hydrogene(self):
        h = table_thermo()["H+(aq)"]
        assert (h.DfH0, h.Sm0, h.Cp0) == (0.0, 0.0, 0.0)

    def test_ion_sans_cp(self):
        assert table_thermo()["Fe^2+(aq)"].Cp0 is None

    def test_phase_etiquetee_cherchee_sous_sa_phase_nue(self):
        assert table_thermo()[Constituant("CH3OH", "l1")].cle == "CH3OH(l)"

    def test_espece_absente_nomme_les_especes_disponibles(self):
        with pytest.raises(KeyError, match="NH3\\(g\\)"):
            table_thermo()["XeF4(g)"]

    def test_chaque_source_est_documentee(self):
        sources = (RACINE / "src" / "pychimie" / "donnees" / "SOURCES.md").read_text("utf-8")
        assert all(f"`{e.source}`" in sources for e in table_thermo())

    def test_lire_un_autre_csv(self, tmp_path):
        chemin = tmp_path / "t.csv"
        chemin.write_text(
            "formule,phase,nom,DfH0_kJ_mol,Sm0_J_K_mol,Cp0_J_K_mol,source\n"
            "Ar,g,argon,0,154.8,20.8,X\nXe+,aq,ion,1,2,,X\n",
            encoding="utf-8",
        )
        t = TableThermo.lire(chemin)
        assert (t["Ar(g)"].Sm0, t["Xe+(aq)"].Cp0) == (154.8, None)

    def test_doublon_refuse(self):
        e = table_thermo()["N2(g)"]
        with pytest.raises(ValueError):
            TableThermo((e, e))


class TestHessValeursDeCours:
    """La loi de Hess sur la table retrouve les valeurs des cours et des manuels."""

    @pytest.mark.parametrize(
        ("equation", "DrH0_kJ", "DrS0"),
        [
            ("H2(g) + Cl2(g) = 2 HCl(g)", -184.6, 20.1),
            ("2 Al(s) + Fe2O3(s) = Al2O3(s) + 2 Fe(s)", -851.5, -38.6),
            ("NaCl(s) = Na+(aq) + Cl-(aq)", 3.9, 43.4),
            ("H+(aq) + OH-(aq) = H2O(l)", -55.8, 80.7),
            ("PCl5(g) = PCl3(g) + Cl2(g)", 87.9, 170.3),
            ("N2(g) + 3 H2(g) = 2 NH3(g)", -92.2, -198.8),
            ("2 SO2(g) + O2(g) = 2 SO3(g)", -197.8, -188.1),
            ("CaCO3(s) = CaO(s) + CO2(g)", 178.3, 160.6),
            ("2 ZnS(s) + 3 O2(g) = 2 ZnO(s) + 2 SO2(g)", -878.3, -147.1),
        ],
    )
    def test_enthalpie_et_entropie(self, equation, DrH0_kJ, DrS0):
        r = Reaction.depuis_tables(equation)
        assert (r.DrH0_ref / 1e3, r.DrS0_ref) == pytest.approx((DrH0_kJ, DrS0), abs=0.06)

    @pytest.mark.parametrize(
        ("equation", "pK"),
        [
            ("AgCl(s) = Ag+(aq) + Cl-(aq)", 9.75),
            ("H2O(l) = H+(aq) + OH-(aq)", 14.0),
            ("CaCO3(s) = Ca^2+(aq) + CO3^2-(aq)", 8.3),
        ],
    )
    def test_constantes_a_298K(self, equation, pK):
        r = Reaction.depuis_tables(equation)
        assert -math.log10(float(r.K0(298.15))) == pytest.approx(pK, abs=0.02)

    def test_kirchhoff_refuse_sans_cp(self):
        with pytest.raises(ValueError, match=r"C°p.*Fe\^2\+\(aq\)"):
            Reaction.depuis_tables("Fe^2+(aq) + Ag+(aq) = Fe^3+(aq) + Ag(s)", kirchhoff=True)

    def test_kirchhoff_passe_avec_cp(self):
        r = Reaction.depuis_tables("NaCl(s) = Na+(aq) + Cl-(aq)", kirchhoff=True)
        assert r.DrCp0 == pytest.approx(46.4 - 136.4 - 50.50)

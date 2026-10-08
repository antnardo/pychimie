import importlib.util
from pathlib import Path

import pytest

from pychimie import Constituant, TableThermo, table_thermo

RACINE = Path(__file__).resolve().parents[1]
ESPECES_CPGE = [
    "N2(g)",
    "H2(g)",
    "O2(g)",
    "NH3(g)",
    "H2O(l)",
    "H2O(g)",
    "CO(g)",
    "CO2(g)",
    "C(s)",
    "CH4(g)",
    "CH3OH(l)",
    "CH3OH(g)",
    "SO2(g)",
    "SO3(g)",
    "NO2(g)",
    "N2O4(g)",
    "CaCO3(s)",
    "CaO(s)",
    "Fe(s)",
    "FeO(s)",
    "Fe2O3(s)",
    "Fe3O4(s)",
    "Cu(s)",
    "CuO(s)",
    "Zn(s)",
    "ZnO(s)",
]


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location(
        "ecrire_donnees_thermo", RACINE / "scripts" / "ecrire_donnees_thermo.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_table_contient_les_especes_des_exercices():
    assert sorted(e.cle for e in table_thermo()) == sorted(ESPECES_CPGE)


def test_csv_nomme_les_especes_du_script(script):
    lu = [(e.formule, e.phase, e.nom, e.source) for e in table_thermo()]
    assert lu == [(f, p, nom, source) for f, p, nom, *_, source in script.ESPECES]


def test_csv_donne_les_valeurs_du_script(script):
    lu = [(e.DfH0 / 1e3, e.Sm0, e.Cp0) for e in table_thermo()]
    assert sum(lu, ()) == pytest.approx(sum((tuple(ligne[3:6]) for ligne in script.ESPECES), ()))


def test_enthalpie_en_joules_par_mole():
    assert table_thermo()["NH3(g)"].DfH0 == pytest.approx(-46110.0)


@pytest.mark.parametrize("cle", ["N2(g)", "H2(g)", "O2(g)", "C(s)", "Fe(s)", "Cu(s)", "Zn(s)"])
def test_corps_simples_de_reference_ont_DfH0_nul(cle):
    assert table_thermo()[cle].DfH0 == 0.0


def test_phase_etiquetee_cherchee_sous_sa_phase_nue():
    assert table_thermo()[Constituant("CH3OH", "l1")].cle == "CH3OH(l)"


def test_espece_absente_nomme_les_especes_disponibles():
    with pytest.raises(KeyError, match="NH3\\(g\\)"):
        table_thermo()["NH4Cl(s)"]


def test_chaque_source_est_documentee():
    sources = (RACINE / "src" / "pychimie" / "donnees" / "SOURCES.md").read_text("utf-8")
    assert all(f"`{e.source}`" in sources for e in table_thermo())


def test_lire_un_autre_csv(tmp_path):
    chemin = tmp_path / "t.csv"
    chemin.write_text(
        "formule,phase,nom,DfH0_kJ_mol,Sm0_J_K_mol,Cp0_J_K_mol,source\nAr,g,argon,0,154.8,20.8,X\n",
        encoding="utf-8",
    )
    assert TableThermo.lire(chemin)["Ar(g)"].Sm0 == 154.8


def test_doublon_refuse():
    e = table_thermo()["N2(g)"]
    with pytest.raises(ValueError):
        TableThermo((e, e))

import pytest

from pychimie.composes import Compose
from pychimie.formules import Formule


@pytest.mark.parametrize(
    ("texte", "formule", "phase"),
    [
        ("H2O(l)", "H2O", "l"),
        ("CO2(g)", "CO2", "g"),
        ("NaCl(s)", "NaCl", "s"),
        ("Cu^2+(aq)", "Cu^2+", "aq"),
        ("H2O (g)", "H2O", "g"),
        ("Ca(OH)2", "Ca(OH)2", None),
        ("Ca(OH)2(s)", "Ca(OH)2", "s"),
    ],
)
def test_lire_separe_la_phase(texte, formule, phase):
    c = Compose.lire(texte)
    assert (c.formule.texte, c.phase) == (formule, phase)


def test_masse_molaire():
    assert Compose.lire("H2O(l)").masse_molaire == pytest.approx(18.015, abs=1e-3)


def test_str_reconstruit_la_notation():
    assert str(Compose.lire("Cu^2+(aq)", nom="ion cuivre(II)")) == "Cu^2+(aq)"


def test_nom():
    assert Compose.lire("H2O(l)", nom="eau").nom == "eau"


def test_phase_inconnue():
    with pytest.raises(ValueError, match="phase inconnue"):
        Compose(Formule.lire("H2O"), phase="vapeur")

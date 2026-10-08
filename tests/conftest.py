import pytest

from pychimie import Reaction


@pytest.fixture
def nh3() -> Reaction:
    """Les valeurs des scripts élèves (Optimisation_NH3_*.py)."""
    return Reaction("N2(g) + 3 H2(g) = 2 NH3(g)", DrH0=-92e3, DrS0=-175.6)


@pytest.fixture
def nh3_tables() -> Reaction:
    return Reaction.depuis_tables("N2(g) + 3 H2(g) = 2 NH3(g)")


@pytest.fixture
def boudouard() -> Reaction:
    """Les valeurs de Optimisation_Boudouard.py."""
    return Reaction("C(s) + CO2(g) = 2 CO(g)", DrH0=172.5e3, DrS0=175.8)


@pytest.fixture
def calcaire() -> Reaction:
    return Reaction.depuis_tables("CaCO3(s) = CaO(s) + CO2(g)")


@pytest.fixture
def dimerisation() -> Reaction:
    return Reaction.depuis_tables("N2O4(g) = 2 NO2(g)")

import doctest

import pychimie


def test_docstring_du_paquet():
    resultat = doctest.testmod(pychimie)
    assert (resultat.failed, resultat.attempted > 0) == (0, True)

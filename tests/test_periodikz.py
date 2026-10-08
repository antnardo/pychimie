import doctest

import pytest

import pychimie.periodikz
from pychimie.elements import TableauPeriodique, element
from pychimie.periodikz import TableauTikz, coordonnees, url_wikipedia


def lignes(t: TableauTikz) -> list[str]:
    return t.texte.splitlines()


def ligne_de(t: TableauTikz, coordonnee: str) -> str:
    (ligne,) = [x for x in lignes(t) if "{" + coordonnee + "}" in x or f"({coordonnee})" in x]
    return ligne


def test_doctests_passent():
    assert doctest.testmod(pychimie.periodikz).failed == 0


@pytest.mark.parametrize(
    ("symbole", "coordonnee"),
    [
        ("H", "1-1"),
        ("He", "18-1"),
        ("Fe", "8-4"),
        ("La", "1-L"),
        ("Yb", "14-L"),
        ("Lu", "3-6"),
        ("Ac", "1-A"),
        ("No", "14-A"),
        ("Lr", "3-7"),
        ("Og", "18-7"),
    ],
)
def test_coordonnees(symbole, coordonnee):
    assert coordonnees(element(symbole)) == coordonnee


def test_coordonnees_toutes_distinctes():
    t = TableauTikz(liens=False).cases()
    assert len(set(lignes(t))) == 118


class TestUrl:
    def test_accents_encodes(self):
        assert url_wikipedia(element("He")) == "https://fr.wikipedia.org/wiki/H%C3%A9lium"

    @pytest.mark.parametrize("symbole", ["Hg", "Ta"])
    def test_homonymes_vers_l_article_de_chimie(self, symbole):
        assert url_wikipedia(element(symbole)).endswith("_%28chimie%29")


class TestCases:
    def test_style_selon_un_attribut(self):
        t = TableauTikz().cases("famille")
        assert ligne_de(t, "17-3") == r"\dessinecase[famille=halogène]{17-3}"

    def test_plusieurs_attributs_et_style_commun(self):
        t = TableauTikz().cases("bloc", "etat", style="thick")
        assert ligne_de(t, "1-1") == r"\dessinecase[thick, bloc=s, etat=g]{1-1}"

    def test_valeur_absente_sans_cle(self):
        t = TableauTikz().cases("etat")
        assert ligne_de(t, "18-7") == r"\dessinecase{18-7}"

    def test_booleen_en_oui_non(self):
        t = TableauTikz().cases("radioactif")
        assert ligne_de(t, "6-4") == r"\dessinecase[radioactif=non]{6-4}"
        assert ligne_de(t, "7-5") == r"\dessinecase[radioactif=oui]{7-5}"

    def test_propriete_calculee(self):
        t = TableauTikz().cases("caractere")
        assert ligne_de(t, "14-3") == r"\dessinecase[caractere=métalloïde]{14-3}"


class TestOmbrage:
    def test_extremes(self):
        t = TableauTikz().ombrage("electronegativite", "chi")
        assert ligne_de(t, "17-2") == r"\dessinecase[chi=100]{17-2}"
        assert ligne_de(t, "1-7") == r"\dessinecase[chi=0]{1-7}"

    def test_valeur_absente_style_par_defaut(self):
        t = TableauTikz().ombrage("electronegativite", "chi", defaut="no chi")
        assert ligne_de(t, "18-1") == r"\dessinecase[no chi]{18-1}"

    def test_valeur_absente_sans_defaut(self):
        t = TableauTikz().ombrage("electronegativite", "chi")
        assert ligne_de(t, "18-1") == r"\dessinecase{18-1}"

    def test_ombrage_est_un_entier(self):
        t = TableauTikz().ombrage("rayon_covalent", "r", bornes=(0, 300))
        assert ligne_de(t, "1-1") == r"\dessinecase[r=10]{1-1}"


class TestSymboles:
    def test_sans_lien(self):
        t = TableauTikz(liens=False).symboles("etat")
        assert ligne_de(t, "1-1") == r"\draw (1-1) coordinate(H) node[symbole, etat=g] {H};"

    def test_lien_wikipedia_echappe(self):
        t = TableauTikz().symboles()
        assert ligne_de(t, "18-1") == (
            r"\draw (18-1) coordinate(He) node[symbole] "
            r"{\href{https://fr.wikipedia.org/wiki/H\%C3\%A9lium}{He}};"
        )

    def test_style_sur_une_selection(self):
        t = TableauTikz(liens=False).symboles(style="emph", selection=["H", "O"])
        assert "node[symbole, emph]" in ligne_de(t, "16-2")
        assert "node[symbole]" in ligne_de(t, "8-4")

    def test_style_sur_tous_sans_selection(self):
        t = TableauTikz(liens=False).symboles(style="emph")
        assert all("emph" in x for x in lignes(t))


class TestInfos:
    def test_modele(self):
        t = TableauTikz().infos("masse_molaire", "M info", r"\num{{{}}}")
        assert ligne_de(t, "8-4") == r"\draw (8-4) node[M info] {\num{55.845}};"

    def test_fonction(self):
        t = TableauTikz().infos("rayon_covalent", "info", lambda r: f"{r / 100:g} Å")
        assert ligne_de(t, "1-1") == r"\draw (1-1) node[info] {0.31 Å};"

    def test_valeurs_absentes_ignorees(self):
        t = TableauTikz().infos("masse_molaire", "M")
        assert len(lignes(t)) == 84

    def test_tuple_vide_ignore(self):
        t = TableauTikz().infos("nombres_oxydation", "ox")
        assert len(lignes(t)) == 100


class TestLegende:
    def test_horizontale_sur_un_rang(self):
        t = TableauTikz()
        n = t.legende("bloc", colonne=1, ligne=2, horizontale=True)
        assert n == 4
        assert lignes(t) == [
            r"\legendhoriz{bloc=s}{1}{2}{S}",
            r"\legendhoriz{bloc=p}{2}{2}{P}",
            r"\legendhoriz{bloc=d}{3}{2}{D}",
            r"\legendhoriz{bloc=f}{4}{2}{F}",
        ]

    def test_verticale_sur_deux_colonnes_arrondi_superieur(self):
        t = TableauTikz()
        n = t.legende("etat", rangs=2)
        assert n == 2
        assert lignes(t) == [
            r"\legend{etat=g}{0}{0}{G}",
            r"\legend{etat=s}{0}{1}{S}",
            r"\legend{etat=l}{1}{0}{L}",
        ]

    def test_plus_de_rangs_que_de_valeurs(self):
        t = TableauTikz()
        assert t.legende("etat", rangs=10) == 3

    def test_exclusion_d_une_valeur_absente_sans_erreur(self):
        t = TableauTikz()
        t.legende("etat", exclues=["l", "plasma"])
        assert [x.split("}{")[0] for x in lignes(t)] == [r"\legend{etat=g", r"\legend{etat=s"]

    def test_ordre_et_libelles_imposes(self):
        t = TableauTikz()
        t.legende("radioactif", valeurs=[True, False], libelles={"oui": "Radioactif"})
        assert lignes(t) == [
            r"\legend{radioactif=oui}{0}{0}{Radioactif}",
            r"\legend{radioactif=non}{0}{1}{Non}",
        ]

    def test_aucune_valeur(self):
        assert TableauTikz().legende("etat", valeurs=[]) == 0

    def test_rangs_invalides(self):
        with pytest.raises(ValueError, match="rang"):
            TableauTikz().legende("etat", rangs=0)


def test_ecrire_en_utf8(tmp_path):
    chemin = TableauTikz().cases("famille").ecrire(tmp_path / "t.tex")
    assert "métal alcalin" in chemin.read_text(encoding="utf-8")


def test_tableau_partiel():
    t = TableauTikz(TableauPeriodique([element("H"), element("He")]), liens=False).cases()
    assert lignes(t) == [r"\dessinecase{1-1}", r"\dessinecase{18-1}"]


def test_couches_dans_l_ordre_des_appels():
    t = TableauTikz(liens=False).cases().symboles()
    assert lignes(t)[0].startswith(r"\dessinecase")
    assert lignes(t)[-1].startswith(r"\draw")

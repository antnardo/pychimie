"""Masses molaires et préparation de solutions, comme à la paillasse.

python preparation.py
"""

from pychimie import Formule, Solution, element, masse_a_peser, masse_molaire

for texte in ("H2O", "H2SO4", "Ca(OH)2", "CuSO4·5H2O", "K4[Fe(CN)6]", "SO₄²⁻"):
    f = Formule.lire(texte)
    print(f"M({texte}) = {f.masse_molaire():.2f} g/mol   (Hill : {f.hill})")

fer = element("fer")
print(
    f"\n{fer} : M = {fer.masse_molaire} ± {fer.incertitude_masse} g/mol, "
    f"fusion à {fer.theta_fusion:.0f} °C"
)

# 100 mL de sulfate de cuivre à 0,10 mol/L, à partir du sel pentahydraté
m = masse_a_peser("CuSO4·5H2O", 0.10, 0.100)
print(f"\nPeser {m:.3f} g de CuSO4·5H2O pour 100 mL à 0,10 mol/L")

# 250 mL d'acide sulfurique à 0,50 mol/L, à partir de l'acide du commerce
commerce = Solution.commerciale("H2SO4", titre_massique=0.96, densite=1.84)
V = commerce.volume_a_prelever(0.50, 0.250)
print(
    f"H2SO4 à 96 % (d = 1,84) : {commerce.concentration:.1f} mol/L ; "
    f"prélever {1000 * V:.1f} mL pour 250 mL à 0,50 mol/L"
)

print(
    f"\nPart massique du cuivre dans CuSO4·5H2O : "
    f"{Formule.lire('CuSO4·5H2O').fractions_massiques()['Cu']:.1%}"
)
print(f"M(NaCl) = {masse_molaire('NaCl'):.2f} g/mol")

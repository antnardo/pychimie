"""K°(T) de la synthèse de l'ammoniac : Ellingham, Kirchhoff et polynômes NASA.

Trois modèles de la même constante, de 300 à 1000 K :

- Ellingham : ΔrH° et ΔrS° constants, pris dans la table à 298,15 K ;
- Kirchhoff : ΔrC°p constant, pris dans la même table ;
- polynômes NASA à sept coefficients, ceux de GRI-Mech 3.0 (``thermo30.dat``), qui
  donnent C°p(T), H°(T) et S°(T) de chaque gaz sur deux plages, 200–1000 K et au-delà.

Les coefficients sont recopiés ci-dessous, sans dépendance à Cantera. Source :
G. P. Smith, D. M. Golden, M. Frenklach, N. W. Moriarty, B. Eiteneer, M. Goldenberg,
C. T. Bowman, R. K. Hanson, S. Song, W. C. Gardiner Jr., V. V. Lissianski, Z. Qin,
*GRI-Mech 3.0* (1999), fichier ``thermo30.dat`` ; origine des données par espèce :
N₂ « 121286 », H₂ « TPIS78 », NH₃ « J 6/77 » (tables JANAF de juin 1977).

À 700 K, les trois modèles donnent 3,2·10⁻⁴, 6,9·10⁻⁵ et 8,8·10⁻⁵ : l'approximation
d'Ellingham surestime K° d'un facteur 4, Kirchhoff à ΔrC°p constant ramène l'écart à
20 % ; au-delà, c'est la variation de ΔrC°p avec T que seuls les polynômes suivent.

    python nh3_ellingham_kirchhoff_nasa.py
"""

import matplotlib.pyplot as plt
import numpy as np

from pychimie import Reaction

# a1 … a7 : C°p/R = a1 + a2 T + a3 T² + a4 T³ + a5 T⁴,
# H°/RT = a1 + a2 T/2 + a3 T²/3 + a4 T³/4 + a5 T⁴/5 + a6/T,
# S°/R = a1 ln T + a2 T + a3 T²/2 + a4 T³/3 + a5 T⁴/4 + a7.
# Première ligne : plage basse (jusqu'à 1000 K) ; seconde : plage haute.
NASA = {
    "N2": (
        (3.298677, 1.4082404e-3, -3.963222e-6, 5.641515e-9, -2.444854e-12, -1020.8999, 3.950372),
        (2.92664, 1.4879768e-3, -5.684760e-7, 1.0097038e-10, -6.753351e-15, -922.7977, 5.980528),
    ),
    "H2": (
        (
            2.34433112,
            7.98052075e-3,
            -1.94781510e-5,
            2.01572094e-8,
            -7.37611761e-12,
            -917.935173,
            0.683010238,
        ),
        (
            3.33727920,
            -4.94024731e-5,
            4.99456778e-7,
            -1.79566394e-10,
            2.00255376e-14,
            -950.158922,
            -3.20502331,
        ),
    ),
    "NH3": (
        (
            4.2860274,
            -4.6605230e-3,
            2.1718513e-5,
            -2.2808887e-8,
            8.2638046e-12,
            -6741.7285,
            -0.62537277,
        ),
        (
            2.6344521,
            5.6662560e-3,
            -1.7278676e-6,
            2.3867161e-10,
            -1.2578786e-14,
            -6544.6958,
            6.5662928,
        ),
    ),
}
T_MILIEU = 1000.0  # K : limite entre les deux plages de GRI-Mech 3.0
NU = {"N2": -1, "H2": -3, "NH3": 2}


def coefficients(espece: str, T: np.ndarray) -> np.ndarray:
    """Les sept coefficients de chaque température, selon sa plage : tableau (7, n)."""
    bas, haut = (np.array(c) for c in NASA[espece])
    return np.where(T < T_MILIEU, bas[:, None], haut[:, None])


def g_sur_RT(espece: str, T: np.ndarray) -> np.ndarray:
    """G°/RT = H°/RT − S°/R d'un gaz, par les polynômes NASA."""
    a1, a2, a3, a4, a5, a6, a7 = coefficients(espece, T)
    h = a1 + a2 * T / 2 + a3 * T**2 / 3 + a4 * T**3 / 4 + a5 * T**4 / 5 + a6 / T
    s = a1 * np.log(T) + a2 * T + a3 * T**2 / 2 + a4 * T**3 / 3 + a5 * T**4 / 4 + a7
    return h - s


def ln_K0_nasa(T: np.ndarray) -> np.ndarray:
    """ln K°(T) = −ΔrG°/RT = −Σ νᵢ G°ᵢ/RT."""
    return -sum(nu * g_sur_RT(espece, T) for espece, nu in NU.items())


ellingham = Reaction.depuis_tables("N2(g) + 3 H2(g) = 2 NH3(g)")
kirchhoff = Reaction.depuis_tables("N2(g) + 3 H2(g) = 2 NH3(g)", kirchhoff=True)
print(ellingham)
print(kirchhoff)

T = np.linspace(300, 1000, 141)
modeles = {
    "Ellingham (ΔrC°p = 0)": ellingham.ln_K0(T),
    "Kirchhoff (ΔrC°p constant)": kirchhoff.ln_K0(T),
    "polynômes NASA (GRI-Mech 3.0)": ln_K0_nasa(T),
}

print(f"\n{'T (K)':>6} {'Ellingham':>12} {'Kirchhoff':>12} {'NASA':>12}")
for t in (298.15, 400, 500, 700, 900, 1000):
    valeurs = [np.exp(float(ellingham.ln_K0(t))), np.exp(float(kirchhoff.ln_K0(t)))]
    valeurs.append(float(np.exp(ln_K0_nasa(np.array([t]))[0])))
    print(f"{t:6.0f} " + " ".join(f"{k:12.3g}" for k in valeurs))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
for nom, lnK in modeles.items():
    ax1.plot(T, lnK / np.log(10), label=nom)
ax1.set(xlabel="$T$ (K)", ylabel=r"$\log_{10} K^\circ$")
ax1.legend()
ax1.grid(True)
# L'écart à NASA, en log10 : ce que coûte chaque approximation
for nom, lnK in list(modeles.items())[:2]:
    ax2.plot(T, (lnK - modeles["polynômes NASA (GRI-Mech 3.0)"]) / np.log(10), label=nom)
ax2.axhline(0, color="k", lw=0.5)
ax2.set(xlabel="$T$ (K)", ylabel=r"$\log_{10}(K^\circ / K^\circ_\mathrm{NASA})$")
ax2.legend()
ax2.grid(True)
fig.tight_layout()
plt.show()

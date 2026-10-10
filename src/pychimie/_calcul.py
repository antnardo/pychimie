"""Deux outils de calcul partagés : un flottant ou un tableau, et une température."""

import numpy as np
from numpy.typing import ArrayLike, NDArray

__all__ = ["Flottants", "sortie", "temperatures"]

Flottants = float | NDArray[np.float64]


def sortie(x: NDArray[np.float64]) -> Flottants:
    """Un tableau, ou un flottant si le calcul portait sur un scalaire."""
    return x[()] if x.ndim == 0 else x


T_MIN = 100.0  # K : au-dessous, la température a toutes les chances d'être en °C.


def temperatures(T: ArrayLike) -> NDArray[np.float64]:
    """``T`` en tableau de kelvins, refusé au-dessous de 100 K (des degrés Celsius ?)."""
    T = np.asarray(T, dtype=float)
    if np.any(T <= 0):
        raise ValueError("la température est en kelvins, strictement positive")
    if np.any(T < T_MIN):
        bas = float(np.min(T))
        raise ValueError(
            f"T est en kelvins : {bas:g} K ? Une température en degrés Celsius s'écrit "
            f"θ + 273,15 ({bas:g} °C = {bas + 273.15:g} K)"
        )
    return T

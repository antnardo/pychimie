"""Deux outils de calcul partagés : un flottant ou un tableau, et une température."""

import numpy as np
from numpy.typing import ArrayLike, NDArray

__all__ = ["Flottants", "sortie", "temperatures"]

Flottants = float | NDArray[np.float64]


def sortie(x: NDArray[np.float64]) -> Flottants:
    """Un tableau, ou un flottant si le calcul portait sur un scalaire."""
    return x[()] if x.ndim == 0 else x


def temperatures(T: ArrayLike) -> NDArray[np.float64]:
    """``T`` en tableau de kelvins, refusé s'il n'est pas strictement positif."""
    T = np.asarray(T, dtype=float)
    if np.any(T <= 0):
        raise ValueError("la température est en kelvins, strictement positive")
    return T

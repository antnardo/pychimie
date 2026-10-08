import numpy as np


class Solution:
    def __init__(self, M: list, c_massique=None, c_molaire=None, densite=None, x_massique=None, x_molaire=None):
        self.M = M  # liste de couples coef stoechio/m [(1, 12.02), (2, 16)]
        if c_molaire is not None:
            self.c_molaire = c_molaire

    @property
    def c_massique(self):
        pass

    @property
    def x_massique(self):
        pass

    @property
    def x_molaire(self):
        pass

    @property
    def densite(self):
        pass
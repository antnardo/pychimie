from constantes import R_GP
from scipy.optimize import bisect
from scipy.integrate import quad

import numpy as np


class Reaction:
    ''' Classe d'abstraction d'une réaction chimique
    dic = {'nom constituant': (coef steochio, 'phase')}
    Phases : s, l, aq, g
    s, aq, et l peuvent avoir des suffixes pour disambiguer

    Mettre les CI après coup :
    Reaction({'C': (-1, 's'), 'CO2': (-1, 'g'), 'CO': (2, 'g')})
    R.set_CI({'C': 1, 'CO2': 1, 'CO': 1})

    DrHo et DrSo peuvent être définis si la connaissance de la cte d'eq est importante
    Sinon Ko peut suffire si T n'est pas pertinent

    Méthodes :
        noms et paramètres explicites

    Unités : T en K, P en bar, V en L

    TODO :
    Faire des tracés auto classiques : Q, DrG, Ko(T), G...
    Faire une base de données de DfHo et Smo pour calculer automatiquement DrHo et DrSo
    Gérer les ruptures d'équilibre ???

    '''
    N = 2000
    EPS = 1e-8  # précision numérique

    def __init__(self, dic, DrHo=1, DrSo=1, Ko=None):
        self.constituants = [k for k, v in dic.items()]
        self.coefs = [v[0] for k, v in dic.items()]
        self.phases = [v[1] for k, v in dic.items()]
        self.CI = [0 for _ in self.constituants]
        self.DrHo = DrHo
        self.DrSo = DrSo
        self._Ko = Ko

    def index(self, C):
        return self.constituants.index(C)

    def set_CI(self, dic):
        for C, n in dic.items():
            self.CI[self.index(C)] = n

    def parametres(self, C):
        i = self.index(C)
        return self.phases[i], self.coefs[i], self.CI[i]

    def xmax(self):
        return min([-n0/nu for n0, nu in zip(self.CI, self.coefs) if nu < 0])

    def xmin(self):
        return max([-n0/nu for n0, nu in zip(self.CI, self.coefs) if nu > 0])

    def is_rupture_max(self):
        L = sorted([(-n0/nu, phase) for n0, nu, phase in zip(self.CI, self.coefs, self.phases) if nu < 0],
                   key=lambda v: v[0])
        return L[0][1][0] == 's'

    def is_rupture_min(self):
        L = sorted([(-n0/nu, phase) for n0, nu, phase in zip(self.CI, self.coefs, self.phases) if nu > 0],
                   key=lambda v: v[0])
        return L[-1][1][0] == 's'

    def is_rupture(self):
        if self.is_rupture_max():
            print("réactif limitant solide ou solvant")
            return True
        elif self.is_rupture_max():
            print("produit limitant solide ou solvant")
            return True
        return False

    def n(self, C, x):
        phase, nu, n0 = self.parametres(C)
        return n0 + nu*x

    def ntot(self, phase, x):
        constituants_phase = [C for C, ph in zip(self.constituants, self.phases) if ph == phase]
        return sum([self.n(C, x) for C in constituants_phase])

    def a(self, C, x, P=1, V=1):
        phase, nu, n0 = self.parametres(C)
        if phase[0] == 's':  # inclus le solvant...
            return 1
        elif phase[0] == 'l':
            return self.n(C, x)/self.ntot(phase, x)
        elif phase == 'aq':
            return self.n(C, x)/V
        elif phase == 'g':
            return self.n(C, x)/self.ntot(phase, x)*P

    def f(self, C, x):
        return self.a(C, x, P=1, V=1)

    def Q(self, x, P=1, V=1):
        @np.vectorize
        def f(x):
            q = 1
            for nu, C in zip(self.coefs, self.constituants):
                q *= self.a(C, x, P, V)**nu
            return q
        return f(x)

    def DrGo(self, T):
        return self.DrHo - T*self.DrSo

    def DrG(self, x, P, T, V=1):
        return R_GP*T*np.log(self.Q(x, P, V)) + self.DrGo(T)

    def liste_avancement(self, xmin=None):
        xmax = self.xmax()
        if xmin is None:
            xmin = self.xmin()
        delta = (xmax-xmin)/self.N
        return np.linspace(xmin+delta, xmax-delta, self.N)  # delta pour éviter les asymptotes

    def compute_G(self, P, T, V=1):
        xi = self.liste_avancement()
        x0 = xi[0]

        def func(x):
            return self.DrG(x, P, T, V)
        G = np.zeros_like(xi)
        for i in range(1, len(xi)):
            x = xi[i]
            G[i] = quad(func, x0, x)
        return (xi, G)

    def Ko(self, T=None):
        if self._Ko is not None:
            return self._Ko
        assert T is not None
        return np.exp(-self.DrGo(T)/(R_GP*T))

    def xeq(self, P=1, T=None):
        def func(x):
            return self.Q(x, P) - self.Ko(T)
        if self.is_rupture():
            print("Attention, une solution n'existe pas forcément")
        xmax = self.xmax()
        xmin = self.xmin()
        d = max(xmax, -xmin)*1e-8  # numerical precision
        # bisect method to avoid going out of the interval
        try:
            return bisect(func, xmin+d, xmax-d)
        except ValueError:
            print('bisect ne fonctionne pas')
            print(P, T, xmax, xmin, self.CI)
            raise ValueError

import re
from pathlib import Path
from decimal import Decimal, InvalidOperation

this_file = Path(__file__)
DATA_ATOMES = this_file.parent/"atomes.csv"

REPHASES = re.compile('((?P<a>[0-9]+,[0-9]+) à )?(?P<pm>([0-9]+,[0-9]+) ± )?(?P<final>~?[0-9]+,?[0-9]*) (?P<unit>[a-z\-0-9·Å]+)( \((?P<phase>.*?)\))?')
RECAS = re.compile('(?P<cas>[0-9]+-[0-9][0-9]-[0-9])( \((?P<phase>.*?)\))?')


class ValeurPhase:
    def __init__(self, value, unit=None, pm=0, phase=None):
        self.value = value
        self.unit = unit if unit is not None else ''
        self.pm = pm
        self.phase = phase

    def __str__(self):
        pm = f' ± {self.pm}' if self.pm > 0 else ''
        phase = f' (self.phase)' if self.phase is not None else ''
        return f'{self.value}{pm} {self.unit}{phase}'

    def __format__(self, format_spec):
        return format_spec.format(self.value)

    def __repr__(self):
        return f'ValeurPhase({self.value},unit="{self.unit}", pm={self.pm}, phase="{self.phase}"'

    def __eq__(self, x):
        if type(x) in [float, int]:
            return x == self.value
        elif isinstance(x, ValeurPhase):
            return x.value == self.value
        else:
            return NotImplementedError

    def __lt__(self, x):
        if type(x) in [float, int]:
            return self.value < x
        elif isinstance(x, ValeurPhase):
            return self.value < x.value
        else:
            return NotImplementedError

    def __gt__(self, x):
        return not (self < x and self == x)

    def __truediv__(self, x):
        if type(x) in [float, int]:
            return self.value / x
        elif isinstance(x, ValeurPhase):
            return self.value / x.value
        else:
            return NotImplementedError

    def __rtruediv__(self, x):
        return 1 / (self/x)


    @classmethod
    def split(cls, s):
        L = []
        if s is None:
            return L
        for m in REPHASES.finditer(s):
            unit = m.group('unit')
            value = float(m.group('final').replace(',', '.').replace('~', ''))
            if m.group('pm') is not None:
                value = float(m.group('pm')[:-3].replace(',', '.'))
                pm = float(m.group('final').replace(',', '.').replace('~', ''))
            elif m.group('a') is not None:
                a = float(m.group('a').replace(',', '.'))
                value = (value+a)/2
                pm = (value-a)/2
            else:
                pm = 0
            phase = m.group('phase')
            L.append(cls(value, unit, pm, phase))
        return L


def fmt_atomes(key, val):
    if key in ['ligne', 'colonne', 'Z']:
        return int(val)
    elif key == 'colonne_2':
        return val
    elif key in ['masse_volumique', 'ionization_energy']:
        return ValeurPhase.split(val)
    elif key == 'CAS':
        try:
            return int(val.split(' ')[0])
        except ValueError:
            return None
    elif key == 'M':
        try:
            return Decimal(val.strip('()').replace(',', '.'))
        except (ValueError, InvalidOperation):
            return None
    elif key[:5] == 'rayon':
        return ValeurPhase.split(val)
    elif key == 'OA':
        return val
    elif key == "oxydation_states":
        L = val.split(', ')
        M = []
        for state in L:
            if state == '':
                return []
            if state[0] == '±':
                x = int(state[1:])
                M.append(x)
                M.append(-x)
            else:
                M.append(int(state))
        return M
    elif key == 'chi':
        try:
            return Decimal(val.replace(',', '.'))
        except (ValueError, InvalidOperation):
            return None
    elif key in ['Tfus', 'Teb']:
        try:
            return Decimal(val[:-3].replace(',', '.'))
        except (ValueError, InvalidOperation):
            return None
    elif key == 'radioactif':
        return val == '1'
    return val


class Atome:
    def __init__(self, **kwargs):
        self.list_attributes = []
        for key, val in kwargs.items():
            val = fmt_atomes(key, val)
            setattr(self, key, val)
            if type(val) is list:
                self.list_attributes.append(key)

    def __repr__(self):
        return f'Atome({self.symbole}, '+', '.join([f'{key}={value}' for key, value in self.__dict__.items()])+')'

    def __hash__(self):
        return hash(self.Z)


class Tableau:
    def __init__(self, cls=Atome):
        self.A = []
        with open(DATA_ATOMES, 'r') as f:
            header = f.readline().strip().split(';')
            for line in f.readlines():
                L = line.strip().split(";")
                self.A.append(cls(**{header[i]: val.strip() for i, val in enumerate(L)}))
        self.A.sort(key=lambda A: A.Z)
        assert self.A[0].Z == 1 and self.A[-1].Z == len(self.A)
        self.list_attributes = self.A[0].list_attributes.copy()

    def __getitem__(self, i):
        return self.A[i-1]

    def __iter__(self):
        yield from self.A

    def match(self, **kwargs):
        assert len(kwargs) == 1
        key = list(kwargs.keys())[0]
        return [A for A in self if getattr(A, key) == kwargs[key]]

    def match_leq(self, **kwargs):
        assert len(kwargs) == 1
        key = list(kwargs.keys())[0]
        return [A for A in self if getattr(A, key) <= kwargs[key]]

    def match_geq(self, **kwargs):
        assert len(kwargs) == 1
        key = list(kwargs.keys())[0]
        return [A for A in self if getattr(A, key) >= kwargs[key]]

    def get(self, *args):
        return [[getattr(A, key) for key in args] for A in self]

    def get_dict(self, *args):
        return [{A.Z: getattr(A, key) for key in args} for A in self]


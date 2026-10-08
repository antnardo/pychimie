from chimie.atomes import Atome, Tableau


def denonize(x):
    return x if x is not None else 'no'


def deemptyize(x):
    return x if x != '' else 'no'


def identity(x):
    return x


def mean(L):
    return sum(L)/len(L)


class AtomeTikz(Atome):
    styles = {
        'etat': ['etat_CSTP', identity],
        'groupe': ['groupe', deemptyize],
        'bloc': ['bloc', identity],
        'origine': ['origine', identity],
        'chi': ['chi', identity],
        'radioactif': ['radioactif', identity],
        'discovery': ['discovery_year', identity],
        'metal': ['metal', deemptyize],
    }

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @property
    def coordinates(self):
        if self.bloc != 'f':
            return f'{self.colonne}-{self.ligne}'
        else:
            return f'{self.colonne-18}-{self.groupe[0].upper()}'

    @property
    def href(self):
        return f'\\href{{https://fr.wikipedia.org/wiki/{self.nom}}}{{{self.symbole}}}'

    def tikz_styles(self, style_list=None):
        L = []
        if style_list is None:
            return L
        for key in style_list:
            attr, func = self.styles[key]
            val = func(getattr(self, attr))
            L.append(f'{key}={val}')
        return L

    def node_symbole(self, attr_style=None, style=None, href=True):
        common_style = [] if style is None else [f'{style}']
        styles = ', '.join(['symbole'] + common_style + self.tikz_styles(attr_style))
        symb = self.symbole if not href else self.href
        return f"""\\draw ({self.coordinates})
    coordinate({self.symbole})
    node[{styles}] {{{symb}}};
"""

    def case(self, attr_style=None):
        styles = ', '.join(self.tikz_styles(attr_style))
        return f"\\dessinecase[{styles}]{{{self.coordinates}}}\n"

    def case_value(self, attr, style, factor=1, default_style=None):
        if attr in self.list_attributes:
            try:
                val = mean([x.value for x in getattr(self, attr)])
            except ZeroDivisionError:
                val = None
        else:
            val = getattr(self, attr)
        if val is None:
            if default_style is None:
                return f"\\dessinecase{{{self.coordinates}}}\n"
            else:
                return f"\\dessinecase[{default_style}]{{{self.coordinates}}}\n"
        else:
            val = int(val*factor)
            return f"\\dessinecase[{style}={val}]{{{self.coordinates}}}\n"

    def node_info(self, attr, style='', format=None):
        val = getattr(self, attr)
        if val is None:
            return ''
        if format is None:
            format = '{}'
        if type(val) is list:
            try:
                val = mean([v.value for v in val])
            except ZeroDivisionError:
                return ''
        return f"\\draw ({self.coordinates}) node[{style}]{{{format.format(val)}}};\n"


class Tikz:
    def __init__(self, href=True):
        self.T = Tableau(AtomeTikz)
        self.href = href
        self.s = ''

    def add_symboles(self, attr_style=None, style=None, symb_style=None):
        orig_style = ', '.join(style) if style is not None else None
        for A in self.T:
            style = orig_style if symb_style is None or A.symbole in symb_style else None
            self.s += A.node_symbole(attr_style=attr_style, style=style, href=self.href)

    def add_cases(self, attr_style=None):
        for A in self.T:
            self.s += A.case(attr_style)

    def add_info(self, attr, style, format=None):
        for A in self.T:
            self.s += A.node_info(attr, style, format=format)

    def add_cases_style(self, attr, style, norm=100, default_style=None):
        # normalisation 0-100
        if attr in self.T.list_attributes:
            L = [x for A in self.T for x in getattr(A, attr)]
        else:
            L = [getattr(A, attr) for A in self.T if getattr(A, attr) is not None]
        m = max(L)
        for A in self.T:
            self.s += A.case_value(attr, style, factor=norm/m, default_style=default_style)

    def add_legend(self, col=0, cols=None, lines=None, line=0, no_vals=None, map=None, horiz=False):
        if cols is None:
            cols = 1
        if lines is None:
            lines = 1
        if type(lines) is str:
            attr = lines
        else:
            assert type(cols) is str
            attr = cols
        if map is None:
            map = dict()
        values = []
        for A in self.T:
            val = getattr(A, attr)
            if val is not None and val != '':
                values.append(val)
        values = list(dict.fromkeys(values))  # keep order
        if no_vals is not None:
            for val in no_vals:
                try:
                    values.remove(val)
                except KeyError:
                    pass
        N = len(values)
        if horiz:
            grouped = lines // N + 1  # nb val par colonne
            cmd = 'legendhoriz'
        else:
            grouped = N // cols  # nb val par colonne
            cmd = 'legend'
        for i, val in enumerate(values):
            x, y = divmod(i, grouped)
            try:
                legend = map[val]
            except KeyError:
                legend = val.capitalize()
            self.s += f'\\{cmd}{{{attr}={val}}}{{{x+col}}}{{{y+line}}}{{{legend}}}\n'
        return N

    def export(self, filename):
        with open(filename, 'w') as f:
            f.write(self.s)
        print(filename, 'succesfully created.')
        print("excerpt:", self.s[:100].replace('\n', '    \n'))

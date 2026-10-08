from chimie.periodikz import Tikz


T = Tikz(href=False)

T.add_cases(with_style=['groupe', 'metal'])
T.add_symboles(with_style=['etat'])
# T.add_info(attr='chi', style='chi info', format='$\\chi=\\num{{{}}}$')
T.add_info(attr='M', style='M info', format='\\num{{{}}}')
T.add_info(attr='nom', style='nom')
T.add_info(attr='Z', style='Z')
T.export('simple_taux.tex')

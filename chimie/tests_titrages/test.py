import matplotlib.pyplot as plt
from mpl_utils import add_legend, add_to_legend
import numpy as np

x = np.linspace(0, 1, 100)

plt.figure()
ax = plt.gca()
lines = plt.plot(x, 2*x, label='A')
c = lines[0].get_color()
plt.plot(x, 3*x, '--', color=c)
plt.plot(x, 4*x, '.', color=c)
lg = plt.legend()
add_to_legend(lg, styles=[dict(ls='--', color=c), dict(ls='', marker='.', color=c)], label="test")


plt.show()

import pylab, numpy as np
from CoolProp.Plots import PropertyPlot
from CoolProp.CoolProp import PropsSI
import CoolProp

fig = pylab.figure(figsize=(8, 4))
ax = fig.add_axes((0, 0, 0.5, 1.0))
ax.fill(np.r_[0, 3, 3, 0, 0], np.r_[0, 0, 1, 1, 0], 'lightblue')
ax.fill(np.r_[0, 3, 3, 0, 0], np.r_[1, 1, 2, 2, 1], 'pink')
ax.plot([2.3, 2.3], [1, 2], 'k--')
ax.plot([2.3, 2.3], [0, 1], 'k')
ax.text(1.15, 0.5, 'Two-Phase', ha='center', va='center')
ax.text(1.15, 1, '1', bbox=dict(facecolor='w'), ha='center', va='center')
ax.text(2.65, 0.5, 'Superheat', ha='center', va='center')
ax.text(2.65, 1, '2', bbox=dict(facecolor='w'), ha='center', va='center')
ax.text(1.5, 1.5, 'Single-Phase', ha='center', va='center')

ax.axis('equal')
ax.axis('off')

ax2 = fig.add_axes((0.6, 0.2, 0.35, 0.7))
# Ts('R134a', axis=ax2)
# PropertyPlot('R134a', 'Ts').calc_isolines(CoolProp.iQ, num=11)
# p = PropsSI('P', 'T', 280, 'Q', 0.0, 'R134a')
# T = np.r_[np.linspace(279.95, 279.99, 2), np.linspace(280.1, 300, 50)]
# s = PropsSI('S', 'T', T, 'P', p, 'R134a') / 1000.0
# ax2.plot(s, T)
# ax2.plot([s[0], s[-1]], [280, 305], 'r')
# ax2.set_xlim(0.9, 1.9)
# ax2.set_ylim(230, 360)
fig.suptitle(r"No internal pinching for $\dot Q_{max,\varepsilon=1}$ - single-phase fluid is hotter")

pylab.show()

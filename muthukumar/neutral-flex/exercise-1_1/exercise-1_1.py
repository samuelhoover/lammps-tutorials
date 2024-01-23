import matplotlib.pyplot as plt
import numpy as np

# Exercise 1
#######################################################################

sigma = 1.0
eps = 1.0

r = np.linspace(0.001, 5, 1000)

plt.figure()

for rc in [1.12, 2.5]:
    if rc == 1.12:
        lc = 'b'
    else:
        lc = 'r'
    
    lj1 = 4 * eps * ((sigma / r) ** 12 - (sigma / rc) ** 12 + (sigma / rc) ** 6 - (sigma / r) ** 6)
    lj2 = 4 * eps * ((sigma / r) ** 12 - (sigma / r) ** 6)
    plt.plot(r, lj1, '{}-'.format(lc), linewidth=0.75, label='With second-term - $R_c$ = {}'.format(rc))
    plt.plot(r, lj2, '{}--'.format(lc), linewidth=0.75, label='W/o second-term - $R_c$ = {}'.format(rc))
    plt.plot(r, [0] * len(r), 'k', linewidth=0.5)
    plt.xlim(0, 5)
    plt.ylim(-2, 10)
    plt.xlabel('$r$')
    plt.ylabel('$U_{LJ}$')

plt.legend()
plt.tight_layout()
plt.show()

import matplotlib.pyplot as plt
import numpy as np

# use LJ units
eps = 1
eps_r = 80
sigma = 1
q = 1
r = np.linspace(1e-6, 4, 1000)
I = np.array([1e-5, 5e-4, 5e-2, 5])
kappa = np.sqrt(I) / 0.3

U_C = (sigma * q ** 2) / (eps_r * r)
U_LJ = 4 * eps * ((1 / r) ** 12 - (1 / r) ** 6)

plt.plot(r, U_C, linewidth=0.75, label='$U_C$')
for i, k in enumerate(kappa):
    plt.plot(
        r,
        ((sigma * q ** 2) / (eps_r * r)) * np.exp(-k * r),
        linewidth=0.75,
        label=rf'$U_D$, I = {I[i]} M',
    )
plt.plot(r, U_LJ, linewidth=0.75, label='$U_{LJ}$')
plt.plot([0, 0], [-1e9, 1e9], 'k', linewidth=0.25)
plt.plot([-10, 10], [0, 0], 'k', linewidth=0.25)
plt.xlabel('$r$')
plt.ylabel('$U$  $[k_BT]$')
plt.ylim(-1.5, 1.5)
plt.xlim(-0.2, 3)
plt.legend()
plt.tight_layout()
plt.savefig('salt-interactions.png', dpi=600)
plt.show()


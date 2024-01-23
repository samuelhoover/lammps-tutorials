import matplotlib.pyplot as plt
import numpy as np


# Exercise 2
#######################################################################

# Part 3
########
t_40 = np.loadtxt('n-40/run1/radius_of_gyration.dat', skiprows=1)[:, 0]
Rg2_40 = np.vstack((
    np.loadtxt('n-40/run1/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('n-40/run2/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('n-40/run3/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('n-40/run4/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('n-40/run5/radius_of_gyration.dat', skiprows=1)[:, 1],
))
mean_Rg_40 = np.sqrt(np.mean(Rg2_40, axis=0))

t_80 = np.loadtxt('n-80/run1/radius_of_gyration.dat', skiprows=1)[:, 0]
Rg2_80 = np.vstack((
    np.loadtxt('n-80/run1/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('n-80/run2/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('n-80/run3/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('n-80/run4/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('n-80/run5/radius_of_gyration.dat', skiprows=1)[:, 1],
))
mean_Rg_80 = np.sqrt(np.mean(Rg2_80, axis=0))

t_120 = np.loadtxt('n-120/run1/radius_of_gyration.dat', skiprows=1)[:, 0]
Rg2_120 = np.vstack((
    np.loadtxt('n-120/run1/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('n-120/run2/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('n-120/run3/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('n-120/run4/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('n-120/run5/radius_of_gyration.dat', skiprows=1)[:, 1],
))
mean_Rg_120 = np.sqrt(np.mean(Rg2_120, axis=0))

t_200 = np.loadtxt('n-200/run1/radius_of_gyration.dat', skiprows=1)[:, 0]
Rg2_200 = np.vstack((
    np.loadtxt('n-200/run1/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('n-200/run2/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('n-200/run3/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('n-200/run4/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('n-200/run5/radius_of_gyration.dat', skiprows=1)[:, 1],
))
mean_Rg_200 = np.sqrt(np.mean(Rg2_200, axis=0))

# plt.figure()
# plt.plot(t_40, mean_Rg_40, label='$N = 40$')
# plt.plot(t_80, mean_Rg_80, label='$N = 80$')
# plt.plot(t_120, mean_Rg_120, label='$N = 120$')
# plt.plot(t_200, mean_Rg_200, label='$N = 200$')
# plt.xlabel('Time')
# plt.ylabel('$R_g$')

# plt.legend()
# plt.tight_layout()
# plt.show()

N = np.array([40, 80, 120, 200])

Rg_40 = np.mean(mean_Rg_40[-600:])
Rg_80 = np.mean(mean_Rg_80[-600:])
Rg_120 = np.mean(mean_Rg_120[-600:])
Rg_200 = np.mean(mean_Rg_200[-600:])
Rg = np.array([Rg_40, Rg_80, Rg_120, Rg_200])

nu = np.mean(np.log(Rg[-1] / Rg[:-1]) / np.log(N[-1] / N[:-1]))

plt.figure()
plt.scatter(N, Rg, label='$\\nu$ = {:.2f}'.format(nu))
plt.xscale('log')
plt.yscale('log')
plt.xlabel('N')
plt.ylabel('$R_g$')

plt.legend()
plt.tight_layout()
plt.show()

for i, j in zip(N, Rg):
    print('N = {}: Rg = {:.3f}'.format(i, j))

print('nu = {:.2f}'.format(nu))

# Rg converges ~3e6, ~4e6, ~1.1e7, and ~1.8e7 steps for N = 40, 80, 120, 200

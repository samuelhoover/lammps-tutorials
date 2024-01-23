import matplotlib.pyplot as plt
import numpy as np


# Exercise 2
#######################################################################

# Part 2
########
t = np.loadtxt('damp-1/run1/radius_of_gyration.dat', skiprows=1)[:, 0]
Rg2_1 = np.vstack((
    np.loadtxt('damp-1/run1/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('damp-1/run2/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('damp-1/run3/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('damp-1/run4/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('damp-1/run5/radius_of_gyration.dat', skiprows=1)[:, 1],
))
mean_Rg_1 = np.sqrt(np.mean(Rg2_1, axis=0))

Rg2_10 = np.vstack((
    np.loadtxt('damp-10/run1/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('damp-10/run2/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('damp-10/run3/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('damp-10/run4/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('damp-10/run5/radius_of_gyration.dat', skiprows=1)[:, 1],
))
mean_Rg_10 = np.sqrt(np.mean(Rg2_10, axis=0))

Rg2_100 = np.vstack((
    np.loadtxt('damp-100/run1/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('damp-100/run2/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('damp-100/run3/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('damp-100/run4/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('damp-100/run5/radius_of_gyration.dat', skiprows=1)[:, 1],
))
mean_Rg_100 = np.sqrt(np.mean(Rg2_100, axis=0))

Rg2_1000 = np.vstack((
    np.loadtxt('damp-1000/run1/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('damp-1000/run2/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('damp-1000/run3/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('damp-1000/run4/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('damp-1000/run5/radius_of_gyration.dat', skiprows=1)[:, 1],
))
mean_Rg_1000 = np.sqrt(np.mean(Rg2_1000, axis=0))

plt.figure()
plt.plot(t, mean_Rg_1, linewidth=0.75, label='$damp$ = 1')
plt.plot(t, mean_Rg_10, linewidth=0.75, label='$damp$ = 10')
plt.plot(t, mean_Rg_100, linewidth=0.75, label='$damp$ = 100')
plt.plot(t, mean_Rg_1000, linewidth=0.75, label='$damp$ = 1000')
plt.xlabel('Time')
plt.ylabel('$R_g$')

plt.legend()
plt.tight_layout()
plt.show()

# Visual analysis on Rg shows that increasing damp causes slower equilibration
# but too low of a damp can cause Rg to fluctuate and never steady.
# Did not print energy and temperature however

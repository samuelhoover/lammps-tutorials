import matplotlib.pyplot as plt
import numpy as np


# Exercise 2
#######################################################################

# Part 1
########
t = np.loadtxt('run1/radius_of_gyration.dat', skiprows=1)[:, 0]
Rg2 = np.vstack((
    np.loadtxt('run1/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('run2/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('run3/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('run4/radius_of_gyration.dat', skiprows=1)[:, 1],
    np.loadtxt('run5/radius_of_gyration.dat', skiprows=1)[:, 1],
))
mean_Rg = np.sqrt(np.mean(Rg2, axis=0))

plt.figure()
plt.plot(t, mean_Rg)
plt.xlabel('Time')
plt.ylabel('$R_g$')

plt.tight_layout()
plt.show()

# Visual analysis shows system equilibrates around 2.75e6 timesteps

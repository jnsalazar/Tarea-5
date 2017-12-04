import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

data = np.genfromtxt("RadialVelocities.dat", skip_header = True, delimiter = " ")
vel_obs = data[:,1]
r_obs = data[:,0]
vel = np.genfromtxt("fit.txt")

plt.scatter(r_obs, vel_obs, label = "Observaciones", s = 2)
plt.plot(r_obs, vel, label = "Ajuste", c = "k")
plt.legend()
plt.xlabel("radio (kpc)")
plt.ylabel("Velocidad (km/s)")
plt.savefig("curva_ajuste.pdf")
plt.close()
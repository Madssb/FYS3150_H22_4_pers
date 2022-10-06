# Program to visualize particle motion
import numpy as np
import matplotlib.pyplot as plt

t, x, y, z = np.loadtxt('data.txt', unpack=True)

plt.plot(t, z)
plt.show()

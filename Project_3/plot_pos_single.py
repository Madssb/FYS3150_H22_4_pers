# Program to visualize particle motion
import numpy as np
import matplotlib.pyplot as plt

t, z_num, z_ana = np.loadtxt('data_test_single.txt', unpack=True)

plt.plot(t, z_num, label='Num')
plt.plot(t, z_ana, label='Ana')
plt.legend()
plt.show()

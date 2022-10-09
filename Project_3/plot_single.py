# Program to visualize particle motion
import numpy as np
import matplotlib.pyplot as plt

t, z_num, z_ana = np.loadtxt('data_test_single.txt', unpack=True)

plt.figure(figsize=(8, 4.5))

plt.title(f'Simulation for $t={np.ceil(t[-1])}\,\mu s$')

plt.plot(t, z_num, color='red', lw=1, label='Num')
plt.plot(t, z_ana, color='black', ls='dashed', lw=.5, label='Ana')

plt.xlabel('t [$\mu s$]')
plt.ylabel('z [$\mu m$]')
plt.legend()

plt.show()
# plt.savefig('single_particle_z_dir.pdf')

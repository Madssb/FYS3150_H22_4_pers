# File to plot figures
import numpy as np
import matplotlib.pyplot as plt

_t, x1, x2, y1, y2, z1, z2, vx1, vx2, vy1, vy2, vz1, vz2 = np.loadtxt('data_test_double.txt', unpack=True)

t = _t[::2]
N = len(t)

r1 = np.zeros((N, 3))
r1[:, 0] = x1[::2]
r1[:, 1] = y1[::2]
r1[:, 2] = z1[::2]

v1 = np.zeros((N, 3))
v1[:, 0] = vx1[::2]
v1[:, 1] = vy1[::2]
v1[:, 2] = vz1[::2]

r2 = np.zeros((N, 3))
r2[:, 0] = x2[::2]
r2[:, 1] = y2[::2]
r2[:, 2] = z2[::2]

v2 = np.zeros((N, 3))
v2[:, 0] = vx2[::2]
v2[:, 1] = vy2[::2]
v2[:, 2] = vz2[::2]

x_min = np.min(r1[:, 0]) - 1
x_max = np.max(r1[:, 0]) + 1
y_min = np.min(r1[:, 1]) - 1
y_max = np.max(r1[:, 1]) + 1

plt.title(f'Particle motion after $t={int(t[-1])}\,\mu s$')

plt.plot(r1[:, 0], r1[:, 1], color='red', lw=1, label="$P_1$")
plt.plot(r2[:, 0], r2[:, 1], color='royalblue', lw=1, label="$P_2$")

plt.xlabel('x [$\mu m$]')
plt.ylabel('y [$\mu m$]')

plt.legend()
plt.show()

# File to plot figures
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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

fig = plt.figure(figsize=(8, 6))
# fig.suptitle(f'Simulation for $t={np.ceil(t[-1])}\,\mu s$\nwith particle interactions')
fig.suptitle(f'Simulation for $t={np.ceil(t[-1])}\,\mu s$\nwithout particle interactions')

ax1 = plt.subplot2grid((2, 4), (0, 0), colspan=2)
ax2 = plt.subplot2grid((2, 4), (0, 2), colspan=2)
ax3 = plt.subplot2grid((2, 4), (1, 0), colspan=4)

ax1.set_title('Motion')
ax1.plot(r1[:, 0], r1[:, 1], color='red', lw=.5, label="$P_1$")
ax1.plot(r2[:, 0], r2[:, 1], color='royalblue', lw=.5, label="$P_2$")
ax1.set_xlabel('x [$\mu m$]')
ax1.set_ylabel('y [$\mu m$]')
ax1.axis('equal')

ax2.set_title('Phase space $(x,v_x)$')
ax2.plot(r1[:, 0], v1[:, 0], color='red', lw=.5)
ax2.plot(r2[:, 0], v2[:, 0], color='royalblue', lw=.5)
ax2.set_xlabel('x [$\mu m$]')
ax2.set_ylabel('$v_x$ [$\mu m/\mu s$]')
ax2.axis('equal')

ax3.set_title('Phase space $(z,v_z)$')
ax3.plot(r1[:, 2], v1[:, 2], color='red', lw=.5)
ax3.plot(r2[:, 2], v2[:, 2], color='royalblue', lw=.5)
ax3.set_xlabel('z [$\mu m$]')
ax3.set_ylabel('$v_z$ [$\mu m/\mu s$]')
ax3.axis('equal')

fig.tight_layout()
plt.figlegend()
# plt.savefig('double_with_interactions.pdf')
plt.savefig('double_without_interactions.pdf')

fig2 = plt.figure(figsize=(8, 6))
ax = plt.axes(projection='3d')

ax.set_title(f'Motion in $(x,y,z)$-plane for $t={np.ceil(t[-1])}\,\mu s$\nwithout particle interactions')
# ax.set_title(f'Motion in $(x,y,z)$-plane for $t={np.ceil(t[-1])}\,\mu s$\nwith particle interactions')
ax.scatter(r1[0, 0], r1[0, 1], r1[0, 2], s=16, color='green', label='$P_{1,start}$')
ax.plot3D(r1[:, 0], r1[:, 1], r1[:, 2], color='red', lw=.75, label='$P_1$')
ax.scatter(r2[0, 0], r2[0, 1], r2[0, 2], s=16, color='black', label='$P_{2,start}$')
ax.plot3D(r2[:, 0], r2[:, 1], r2[:, 2], color='royalblue', lw=.75, label='$P_2$')
ax.scatter(r1[-1, 0], r1[-1, 1], r1[-1, 2], s=16, color='red')
ax.scatter(r2[-1, 0], r2[-1, 1], r2[-1, 2], s=16, color='royalblue')
ax.set_xlabel('x [$\mu m$]')
ax.set_ylabel('y [$\mu m$]')
ax.set_zlabel('z [$\mu m$]')
ax.legend()

plt.savefig('3d_double_without_interactions.pdf')
# plt.savefig('3d_double_with_interactions.pdf')
plt.show()

# Program to visualize particle motion
import numpy as np
import matplotlib.pyplot as plt

w_z = np.sqrt(2 * 2.41e6 / (40.078 * 500**2))

t, x_num, y_num, z_num, x_ana, y_ana, z_ana = np.loadtxt('data_test_single.txt', unpack=True)
dt = t[1] - t[0]

fig = plt.figure(figsize=(8 * 1.4, 4.5 * 1.4))
fig.suptitle('Position vector $\mathbf{r}=(x,y,z)$' + f'\n$dt=${dt:.1e}$\,\mu s$')

ax1 = plt.subplot2grid((2, 4), (0, 0), colspan=2)
ax2 = plt.subplot2grid((2, 4), (0, 2), colspan=2)
ax3 = plt.subplot2grid((2, 4), (1, 1), colspan=2)

ax1.set_title('Position in x-direction')
ax1.plot(t, x_num, color='red', lw=1, label='Approx')
ax1.plot(t, x_ana, color='black', ls='dashed', lw=1, label='Exact')
ax1.set_xlabel('Time [$\mu s$]')
ax1.set_ylabel('Position [$\mu m$]')

ax2.set_title('Position in y-direction')
ax2.plot(t, y_num, color='red', lw=1)
ax2.plot(t, y_ana, color='black', ls='dashed', lw=1)
ax2.set_xlabel('Time [$\mu s$]')
ax2.set_ylabel('Position [$\mu m$]')


ax3.set_title('Position in z-direction')
ax3.plot(t, z_num, color='red', lw=1)
ax3.plot(t, z_ana, color='black', ls='dashed', lw=1)
ax3.set_xlabel('Time [$\mu s$]')
ax3.set_ylabel('Position [$\mu m$]')


plt.figlegend(loc='lower right')
fig.tight_layout()

plt.figure()
plt.title(f'Position in z-direction\n$dt=${dt:.1e}$\,\mu s$, $\omega_z=${w_z:.3e} ' + '$(\mu s)^{-1}$')
plt.plot(t, z_num, color='red', lw=1)
plt.plot(t, z_ana, color='black', ls='dashed', lw=1, label='Exact')
plt.xlabel('Time [$\mu s$]')
plt.ylabel('Position [$\mu m$]')
plt.figlegend(loc='center right')
# plt.savefig('single_zdir.pdf')
# plt.savefig('single_zdir_with_exact.pdf')

plt.show()

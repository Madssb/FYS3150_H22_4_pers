# Plot how many particles remains after a 500 mu*s simulation
import numpy as np
import matplotlib.pyplot as plt

_f1, _f2, _f3 = np.loadtxt('no_of_particles.txt', usecols=(1, 2, 3), unpack=True)
f1, f2, f3 = _f1[0], _f2[0], _f3[0]

w_V, N1, N2, N3 = np.loadtxt('no_of_particles.txt', usecols=(0, 4, 5, 6), unpack=True)

plt.figure(figsize=(8, 4.5))
# plt.title('Particles in trap after $500\,\mu s$\nwithout particle interactions')
plt.title('Particles in trap after $500\,\mu s$\nwith particle interactions')
plt.plot(w_V, N1, color='red', lw=1, label=f'$f={f1:.1f}$')
plt.plot(w_V, N2, color='royalblue', lw=1, label=f'$f={f2:.1f}$')
plt.plot(w_V, N3, color='black', lw=1, label=f'$f={f3:.1f}$')
plt.xlabel('$\omega_V$ [MHz]')
plt.ylabel('$\\frac{N}{N_{tot}}$', fontsize=16)
plt.legend()
# plt.ylim([0, 1.5])
# plt.savefig('particles_remaining_nointeraction.pdf')
# plt.savefig('particles_remaining_nointeraction_finegrain.pdf')
# plt.savefig('particles_remaining_finegrain.pdf')
plt.show()

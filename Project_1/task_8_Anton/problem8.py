import numpy as np
import matplotlib.pyplot as plt

fig1 = plt.figure()
fig2 = plt.figure()
fig3 = plt.figure()
ax1 = fig1.add_subplot()
ax2 = fig2.add_subplot()
ax3 = fig3.add_subplot()
# axin = ax.inset_axes([0.2, 0.1, 0.2, 0.2])


N = 7
n_arr = np.array([10**i for i in range(1, N+1)])
maxrelerror = np.zeros(len(n_arr))
maxabserror = np.zeros(len(n_arr))
for i, n in enumerate(n_arr):
    data = np.loadtxt(f'problem7_data_n={n}.txt')   # Collect numerical data
    # Excluding endpoints that contains 0
    x = data[:,0][1:-1]
    v = data[:,1][1:-1]
    u = 1 - (1 - np.exp(-10))*x - np.exp(-10*x)     # Analytic solution

    abserror = abs(u-v)
    relerror = abs((u - v)/u)

    ax1.plot(x, abserror, label='n = 10' + r'$^{' + f'{i+1}' + r'}$')
    ax2.plot(x, relerror, label='n = 10' + r'$^{' + f'{i+1}' + r'}$')

    maxrelerror[i] = np.max(relerror)
    maxabserror[i] = np.max(abserror)
    print(f'|\t{n}\t|\t{maxrelerror[i]:.3e}\t|')

ax3.plot(n_arr, maxrelerror, 'r')
# ax3.plot(n_arr, maxabserror, 'b')

ax1.set_yscale('log')
ax2.set_yscale('log')
ax3.set_yscale('log')
ax3.set_xscale('log')

ax1.set_title('Absolute error', weight='bold', fontsize=20)
ax2.set_title('Relative error', weight='bold', fontsize=20)
ax3.set_title('Maximum relative error', weight='bold', fontsize=20)

ax1.set_xlabel('x', weight='bold', fontsize=30)
ax2.set_xlabel('x', weight='bold', fontsize=30)
ax3.set_xlabel(r'n$_{steps}$', weight='bold', fontsize=30)

ax1.set_ylabel(r'$\Delta$', weight='bold', fontsize=30)
ax2.set_ylabel(r'$\epsilon$', weight='bold', fontsize=30)
ax3.set_ylabel(r'max($\epsilon$)', weight='bold', fontsize=30)

ax1.legend(prop={'size': 20})
ax2.legend(prop={'size': 20})

fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()

plt.show()

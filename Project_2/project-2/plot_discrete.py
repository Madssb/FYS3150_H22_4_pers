# Plot the discrete solution solved by jacobi_eigensolver in source.cpp
import numpy as np
import matplotlib.pyplot as plt

def plot_discrete(filename):
    x, v_1, v_2, v_3 = np.loadtxt(filename, unpack=True)
    n = len(x) - 1

    v_list = np.array([v_1, v_2, v_3])
    v_colors = np.array(['black', 'royalblue', 'red'])

    fig, ax = plt.subplots(3, 1, figsize=(9,6), sharex=True)

    for i in range(len(v_list)):
        ax[i].plot(x, np.zeros(len(x)), linestyle='dashed', color='dimgrey', lw=.8, label=f'$v_{int(i+1)}$=0')
        ax[i].plot(x, v_list[i], color=v_colors[i], label=f'$v_{int(i+1)}$')
        ax[i].set_ylabel(f'$v_{int(i+1)}$')
        ax[i].legend()

    ax[0].set_title(f'n={int(n)}')
    ax[-1].set_xlabel('$\hat{x}$')
    plt.tight_layout()

plot_discrete('3_eigvecs_N_is_9_DONT_MODIFY.txt')
# plt.savefig('n_is_10.pdf')

plot_discrete('3_eigvecs_N_is_99_DONT_MODIFY.txt')
# plt.savefig('n_is_100.pdf')

plt.show()

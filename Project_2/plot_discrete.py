# Plot the discrete solution solved by jacobi_eigensolver in source.cpp
import numpy as np
import matplotlib.pyplot as plt

def plot_discrete(filename_1, filename_2):
    x, v_1, v_2, v_3 = np.loadtxt(filename_1, unpack=True)
    ana_v_1, ana_v_2, ana_v_3 = np.loadtxt(filename_2, unpack=True)
    n = len(x) - 1

    v_list = np.array([v_1, v_2, v_3])
    ana_v_list = np.array([ana_v_1, ana_v_2, ana_v_3])
    v_colors = np.array(['black', 'royalblue', 'red'])

    fig, ax = plt.subplots(3, 1, figsize=(9,6), sharex=True)

    for i in range(len(v_list)):
        ax[i].plot(x, np.zeros(len(x)), linestyle='dashed', color='dimgrey', lw=.8)
        ax[i].plot(x, ana_v_list[i], color=v_colors[i], linestyle='dashed', label='Analytical')
        ax[i].plot(x, v_list[i], color=v_colors[i], label='Numerical')
        ax[i].set_ylabel(f'$v_{int(i+1)}$')
        ax[i].legend()

    ax[0].set_title(f'n={int(n)}')
    ax[-1].set_xlabel('$\hat{x}$')
    plt.tight_layout()

plot_discrete('3_eigvecs_N_is_9.txt', '3_ANA_eigvecs_N_is_9.000000.txt')
plt.savefig('n_is_10.pdf')

plot_discrete('3_eigvecs_N_is_99.txt', '3_ANA_eigvecs_N_is_99.000000.txt')
plt.savefig('n_is_100.pdf')
plt.show()

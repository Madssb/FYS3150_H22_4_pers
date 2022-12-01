# In this program we define a potential matrix
import numpy as np
import matplotlib.pyplot as plt

def potential(n_slits, out_filename, spatial_step_length=.005, V_0=1e10, thickness=.02, center_length=.05, apperture=.05, save=None):

    out_filename += '.dat'

    M = 1 / spatial_step_length + 1
    N = int(np.ceil(M - 2))
    V = np.zeros((N, N))

    h = spatial_step_length

    c_x = N // 2                        # Center index along x axis
    c_y = N // 2                        # Center index along y axis
    c_l = int(center_length * N)        # Index length of center wall
    a = int(apperture * N)              # Index length of apperture
    t = int(thickness * N)              # Index length of wall thickness

    e_l = N - int(n_slits) * a - c_l    # Index length of end wall pieces
    e_l //= 2

    for i in range(N):
        for j in range(c_y - t // 2, c_y + t // 2 + 1):

            if (i < e_l
            or  e_l + a <= i <= e_l + a + c_l
            or  e_l + 2 * a + c_l < i):

                V[i, j] = V_0

    if save == True:

        np.savetxt('../' + out_filename, V)

    else:

        plt.imshow(V)
        plt.show()

# potential(2, 'zero_potential', V_0=0, save=True)
potential(2, 'zero_potential', V_0=0)

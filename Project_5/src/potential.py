# In this program we define a potential matrix
import numpy as np

def potential(n_slits, spatial_step_length=.005, V_0=1e10, thickness=.02, center_length=.05, apperture=.05):

    M = 1 / spatial_step_length + 1
    N = int(np.ceil(M - 2))

    V = np.zeros((N, N))

    center_walls = n_slits - 1
    end_length = .5 * (center_walls * center_length + n_slits * apperture)

    x_min = .5 * (1 - thickness)
    x_max = .5 * (1 + thickness)

    center_start = .5 * (1 - center_length)
    center_stop = .5 * (1 + center_length)

    if n_slits == 2:

        for j in range(N):

            yi = (j + 1) * spatial_step_length

            for i in range(N):

                xi = (i + 1) * spatial_step_length

                if x_min <= xi <= x_max:
                    if yi <= end_length or center_start <= yi <= center_stop \
                    or end_length + n_slits * apperture + center_length <= yi:

                        V[j, i] = V_0

    print(V)

potential(2, spatial_step_length=.1, thickness=.3, center_length=.3, apperture=.1)

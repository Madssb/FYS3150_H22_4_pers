# In this program we define a potential matrix
import numpy as np
import matplotlib.pyplot as plt

def potential(n_slits, spatial_step_length=.005, V_0=1e10, thickness=.02, center_length=.05, apperture=.05):

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

    #
    # center_walls = n_slits - 1
    # end_length = .5 * (center_walls * center_length + n_slits * apperture)
    #
    # x_min = .5 * (1 - thickness)
    # x_max = .5 * (1 + thickness)
    #
    # center_start = .5 * (1 - center_length)
    # center_stop = .5 * (1 + center_length)
    #
    # if n_slits == 2:
    #
    #     for i in range(N):
    #
    #         yi = (i + 1) * spatial_step_length
    #
    #         for j in range(N):
    #
    #             xi = (j + 1) * spatial_step_length
    #
    #             if x_min < xi < x_max:
    #                 if yi < end_length or center_start < yi < center_stop \
    #                 or end_length + n_slits * apperture + center_length < yi:
    #
    #                     V[i, j] = V_0
    if N == 9:
        print(V)
        print(f'M={M}, V={N}x{N}')
    # plt.imshow(V)
    # plt.show()
    # np.savetxt(f'../{n_slits}_slits_{N}x{N}_potential.dat', V)

potential(2)

# potential(2, spatial_step_length=.1, thickness=.05, center_length=.25, apperture=.2)

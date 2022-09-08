import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

x = pa.mat()
u = pa.mat()
v = pa.mat()

x.load("x.bin")
u.load("u.bin")
v.load("v.bin")

plt.plot(x, v, color='red', label='$v(x)$')
plt.plot(x, u, lw=.9, linestyle='dashed', color='royalblue', label='$u(x)$')
plt.legend()
plt.show()

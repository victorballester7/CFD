# -*- coding:utf-8 -*-
#
# Advection-Diffusion and stability
# Emmanuel Dormy
# Numerical Methods for Fluid Dynamics, ENS 2023,
# TD2
#

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as lg
import matplotlib.pyplot as plt
# plt.ion()


# Size of the domain
L = 5
# Kinematic viscosity
nu = 1.0
# Advection flow
UO = 25.0

# Numerical parameters
dt = 2.e-3         # Time step
NP = 50            # Number of pts in the domain, including bdy pts
deltax = L / (NP - 1)  # Grid step
N = NP - 1           # Periodic BC

# Matrix construction for the discrete operator
a = (nu / deltax**2 + 0.5 * (UO / deltax))
b = -(2 * nu / deltax**2)
c = (nu / deltax**2 - 0.5 * (UO / deltax))
subd = a * np.ones(N)
d = b * np.ones(N)
supd = c * np.ones(N)

# Periodic BC (we think it as a 5-diagonal matrix)
C = [supd, subd, d, supd, subd]
offsets = [-N + 1, -1, 0, 1, N - 1]
T = sp.dia_matrix((C, offsets), shape=(N, N))

print(T.todense())

# Discrete eigenvalues via the matrix
p, v = np.linalg.eig(T.todense())

# Discrete eigenvalues via analytics
#   Because of Nyquist: kmax=N pi / L
kdeltax = 2 * np.pi * np.linspace(np.ceil(-N / 2), np.floor(N / 2), N) / N
pp = -UO * 1j * np.sin(kdeltax) / deltax + 2 * nu * \
    (np.cos(kdeltax) - 1) / deltax**2

# Eigenvalues for the continuous equation
kkk = np.linspace(0, np.pi / deltax, 1000)
ppp = -UO * 1j * kkk - nu * (kkk**2)


m = p * dt
mm = pp * dt
mmm = ppp * dt


plt.figure(1)

fig_min = (m.real).min()
fig_min = np.minimum(fig_min, -2)
fig_min = np.floor(fig_min) - 1
plt.plot(np.real(mmm), np.imag(mmm), '.g')
plt.plot(np.real(mm), np.imag(mm), 'xr')
plt.plot(np.real(m), np.imag(m), '+b')
theta = np.linspace(0, 2 * np.pi, 100)
plt.plot(-1 + np.cos(theta), np.sin(theta), 'k--')
plt.axis([fig_min, 0.1, -1.1, 1.1])
plt.show()

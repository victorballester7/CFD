# -*- coding:utf-8 -*-
#
#       1D Finite Difference approximation of the wave equation
#       with non-reflecting boundary conditions
#       (one-way wave equation: u_t = u_x)
#
#       Emmanuel Dormy (2021)
#

import numpy as np

import matplotlib.pyplot as plt

# Domain length and resolution
Len = 1
N = 300
dx = Len / (N - 2)
x = np.linspace(-dx / 2, Len + dx / 2, N)

# Time step
dt = 0.8 * dx
NT = 300

# Initialization
pold = np.zeros((N, 1))
p = np.zeros((N, 1))
pnew = np.zeros((N, 1))
pxx = np.zeros((N, 1))
s = np.zeros((N, 1))

# Initial Condition
x0 = Len / 2
delta = 1 / 100
rr2 = (x - x0)**2
p[:, 0] = 0.5 * np.exp(-rr2 / delta)
pold = p.copy()

# Time
time = 0

# Time-loop
for i in range(NT):

    print('Iteration ', i, '/', NT)

# Laplacian of p
    pxx[1:N - 1] = (p[0:N - 2] - 2 * p[1:N - 1] + p[2:N])

# Discrete wave equation
    pnew = (dt * dt / dx / dx) * pxx + 2 * p - pold

# Nonreflecting boundary conditions (Engquist-Majda)


# reflecting boundary conditions (Neumann)
    # pnew[0] = pnew[1]
    # pnew[N - 1] = pnew[N - 2]

# reflecting boundary conditions (Dirichlet)
    # pnew[0] = 0
    # pnew[N - 1] = 0

# update
    pold = p.copy()
    p = pnew.copy()
    time = time + dt

# Visualisation
    if (i % 5 == 0):

        plt.figure(1)
        plt.clf()
        plt.plot(x[1:N - 1], p[1:N - 1])
        plt.axis([0, 1, -0.6, 0.6])
        plt.title("{:10.2f}".format(time))
        plt.draw()
        plt.pause(0.001)

#   if (i%20==0):
#      plt.ginput()     # input graphique, pour pause

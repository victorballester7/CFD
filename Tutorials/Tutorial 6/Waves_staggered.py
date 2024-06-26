# -*- coding:utf-8 -*-
#
#       Wave equation on a staggered mesh
#
#       Emmanuel Dormy (2021)
#

import numpy as np


# affichage graphique
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


fig = plt.figure(1, figsize=[10, 7])
ax = fig.add_subplot(projection='3d')


# Physical parameters
xLen = 3
yLen = 1
c = 1  # Speed of sound

# Numerical parameters
Nx = 300
Ny = 100

# staggered mesh
dx = xLen / Nx
dy = yLen / (Ny - 1 / 2)

NT = 200  # Number of time steps

# Inits of arrays and variables
p = np.zeros((Nx, Ny))
u = np.zeros((Nx + 1, Ny))
v = np.zeros((Nx, Ny))
time = 0

# Spatial structure of source term
# Pressure mesh
x = np.linspace(dx / 2, xLen - dx / 2, Nx)
y = np.linspace(dy / 2, yLen, Ny)

xx, yy = np.meshgrid(x, y)

# Source term
x0 = xLen / 4
y0 = yLen / 2
delta = 1 / 50
rr = np.sqrt((xx - x0)**2 + (yy - y0)**2)
rr = np.exp(-rr / delta).transpose()

# Amplitude and time period of source term
amp = 20
T = 0.2


dt = 0.01

# CFL restriction
dt = 0.8 * min(dx, dy) / c / np.sqrt(2)


for it in range(NT):

    print('iteration ', it, ' out of ', NT)
    print('time ', time)

    # First step update velocity (with Leapfrog)
    u[1:Nx, :] = u[1:Nx, :] - dt / dx * \
        c ^ 2 * (p[2:Nx + 1, :] - p[0:Nx - 1, :])
    v[:, 1:Ny] = v[:, 1:Ny] - dt / dy * \
        c ^ 2 * (p[:, 2:Ny + 1] - p[:, 0:Ny - 1])

    # Update boundary points
    # periodic boundary conditions in x, no flux at y = 0
    u[0, :] = u[0, :] - dt / dx * c ^ 2 * (p[1, :] - p[Nx, :])
    u[Nx, :] = u[0, :]
    v[:, 0] = 0

    # Compute pressure source term
    s = rr * amp * np.sin(2 * np.pi * time / T)

    # Second step update pressure

    p[1:Nx - 1,
      1:Ny - 1] = p[1:Nx - 1,
                    1:Ny - 1] - dt * ((u[0:Nx - 1,
                                       1:Ny - 1] - u[2:Nx + 1,
                                       1:Ny - 1]) / (2 * dx) + (v[1:Nx - 1,
                                                                  0:Ny - 2] - v[1:Nx - 1,
                                                                                2:Ny + 1]) / (2 * dy)) + s[1:Nx - 1,
                                                                                                           1:Ny - 1]

    # Update boundary points
    # p = 0 at y = Ly
    p[:, Ny] = 0
    # p periodic in x
    p[0, :] = p[0, :] - dt * ((u[Nx, :] - u[1, :]) /
                              (2 * dx) + (v[0, :] - v[0, :]) / (2 * dy)) + s[0, :]
    p[Nx, :] = p[0, :]

    # TO BE COMPLETED

    time = time + dt

    # Visualisation
    if (it % 5 == 0):
        plt.figure(1)
        plt.clf()

        ax = fig.add_subplot(1, 1, 1)
        ax = fig.gca(projection='3d')
        surf = ax.plot_wireframe(xx, yy, p.transpose(), rstride=5, cstride=5)
        ax.set_zlim(-0.1, 0.1)
        ax.view_init(60, -60)

        plt.axis('image')
        plt.draw()
        plt.pause(0.1)

        plt.figure(2)
        plt.clf()
        plt.pcolor(xx, yy, p.transpose(), cmap=cm.coolwarm, antialiased=False)
        plt.axis('image')
        plt.clim(-0.1, 0.1)
        plt.draw()
        plt.pause(0.1)

        plt.figure(3)
        plt.clf()
        plt.plot(y, p[int(Nx / 4), :])
        plt.draw()
        plt.pause(0.1)

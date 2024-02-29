# -*- coding:utf-8 -*-
#
#       Advection test with a semi-lagrangial scheme
#
#       Emmanuel Dormy (2021)
#

import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as lg
import matplotlib.pyplot as plt


####
def Semilag2(u, v, q):
    """
    Second order semi-Lagrangian advection
    """

    ADVq = Semilag(u, v, q)
    aux = Semilag(-u, -v, ADVq)
    ADVq = q + (q - aux) / 2.
    ADVq = Semilag(u, v, ADVq)

    return ADVq


####
def Semilag(u, v, q):
    """
    1st order semi-Lagrangian advection
    """
    ADVq = np.zeros((NX, NY))

    # Matrices where 1 is right, 0 is left or center
    Mx2 = np.sign(np.sign(u[1:-1, 1:-1]) + 1.)
    Mx1 = 1. - Mx2

    # Matrices where 1 is up, 0 is down or center
    My2 = np.sign(np.sign(v[1:-1, 1:-1]) + 1.)
    My1 = 1. - My2

    # Matrices of absolute values for u and v
    au = abs(u[1:-1, 1:-1])
    av = abs(v[1:-1, 1:-1])

    # Matrices of coefficients respectively central, external, same x, same y
    Cc = (dx - au * dt) * (dy - av * dt) / dx / dy
    Ce = dt * dt * au * av / dx / dy
    Cmx = (dx - au * dt) * av * dt / dx / dy
    Cmy = dt * au * (dy - dt * av) / dx / dy

    # Computes the advected quantity
    ADVq[1:-1, 1:-1] = (Cc * q[1:-1, 1:-1] +
                        Ce * (Mx1 * My1 * q[2:, 2:] +
                              Mx2 * My1 * q[:-2, 2:] +
                              Mx1 * My2 * q[2:, :-2] +
                              Mx2 * My2 * q[:-2, :-2]) +
                        Cmx * (My1 * q[1:-1, 2:] +
                               My2 * q[1:-1, :-2]) +
                        Cmy * (Mx1 * q[2:, 1:-1] +
                               Mx2 * q[:-2, 1:-1]))

    return ADVq

#########################################
######        MAIN PROGRAM         ######
#########################################


# Domain Size
# aspect_ratio = LY/LX

aspect_ratio = float(1.)
LY = float(1.)
LX = LY / aspect_ratio

# Grid Size
# (incuding ghost points)

NX = int(102)
NY = int(100)

# Size of the domain without ghost points
nx = NX - 2
ny = NY - 2

# Number of time iterations
nitermax = int(10001)

# Plot every 1000 iterations
modulo = int(1000)

# ARRAYS & INITIAL CONDITIONS

# Arrays definitions
u = np.zeros((NX, NY))
v = np.zeros((NX, NY))
q = np.zeros((NX, NY))
ADVq = np.zeros((NX, NY))

# Finite difference coefficients
dx = LX / (nx)
dy = LY / (ny)

dx_2 = 1. / (dx * dx)
dy_2 = 1. / (dy * dy)

# Grid, for plotting only
x = np.linspace(dx / 2, LX - dx / 2, nx)
y = np.linspace(dy / 2, LY - dy / 2, ny)
[xx, yy] = np.meshgrid(x, y)

# Advection field
u[1:-1, 1:-1] = np.sin(np.pi * np.transpose(xx) / LX)**2 * \
    np.sin(2 * np.pi * np.transpose(yy) / LY)
v[1:-1, 1:-1] = -np.sin(2 * np.pi * np.transpose(xx) / LX) * \
    (np.sin(np.pi * np.transpose(yy) / LY)**2)

# Initial field to be advected
r = np.sqrt((np.transpose(xx) - 0.5)**2 + (np.transpose(yy) - 0.75)**2)
r = np.where(r > 0.2, 0.2, r)
r = r / 0.2
q[1:-1, 1:-1] = (1 + np.cos(np.pi * r)) / 4
q0 = q


# Time step (could be improved using the CFL condition)
dt = 0.0001
t = 0.  # Simulated time

################
# MAIN LOOP
################

for niter in range(nitermax):

    t += dt

    # semi-Lagrangian advection
    q = Semilag2(u, v, q)

    if ((niter + 1) % modulo == 0):

        # PLOT
        plotlabel = "t = %1.5f" % (t)
        plt.clf()
        plt.title(plotlabel)
        plt.pcolor(xx, yy, np.transpose(q[1:-1, 1:-1]))
        plt.axis('image')
        plt.draw()
        plt.pause(0.1)

# reverse flow: u -> -u, v -> -v

for niter in range(nitermax):

    t += dt

    # semi-Lagrangian advection
    q = Semilag2(-u, -v, q)

    if ((niter + 1) % modulo == 0):

        # PLOT
        plotlabel = "t = %1.5f" % (t)
        plt.clf()
        plt.title(plotlabel)
        plt.pcolor(xx, yy, np.transpose(q[1:-1, 1:-1]))
        plt.axis('image')
        plt.draw()
        plt.pause(0.1)

# Final plot

plotlabel = "t = %1.5f" % (t)
plt.clf()
plt.title(plotlabel)
plt.pcolor(xx, yy, np.transpose(q[1:-1, 1:-1]))
plt.axis('image')
plt.draw()


plt.figure()
plt.title("error")
plt.pcolor(xx, yy, np.transpose(q[1:-1, 1:-1] - q0[1:-1, 1:-1]))
plt.axis('image')
plt.colorbar()
plt.show()

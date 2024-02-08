# -*- coding:utf-8 -*-
#
# Convergence of Finite Differences on a 2D Poisson problem
# Emmanuel Dormy
# Numerical Methods for Fluid Dynamics, ENS 2023,
# TD1
#

from sys import exit
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as lg
import matplotlib.pyplot as plt
# plt.ion() # I don't need .ion() if I use plt.show() instead of
# plt.draw() and plt.pause()


def BuildLaPoisson():
    """
    Matrice du Laplacien 2D
    """

    # Definition of the 1D Laplace operator

    NXi = NX - 2
    NYi = NY - 2

    # AXE X
    # Diagonal terms
    dataNX = [np.ones(NXi), -2 * np.ones(NXi), np.ones(NXi)]

    # AXE Y
    # Diagonal terms
    dataNY = [np.ones(NYi), -2 * np.ones(NYi), np.ones(NYi)]

    # Their positions
    # this is to set the position of the array with respect to the diagonal
    offsets = np.array([-1, 0, 1])
    DXX = sp.dia_matrix((dataNX, offsets), shape=(NXi, NXi)) * dx_2
    DYY = sp.dia_matrix((dataNY, offsets), shape=(NYi, NYi)) * dy_2
    # print(DXX.todense())
    # print(DYY.todense())

    # 2D Laplace operator
    # kronecker product id \times DXX + DYY \times id
    LAP = sp.kron(sp.eye(NYi, NYi), DXX) + sp.kron(DYY, sp.eye(NXi, NXi))

    # print(LAP.todense())
    return LAP


def LUdecomposition(LAP):
    """
    returns the LU decomposition of a sparse matrix LAP
    """
    return lg.splu(LAP.tocsc(),)


def ResoLap(splu, RHS):
    """
    solve the system

    SPLU * x = RHS

    Args:
    --RHS: 2D array((NY,NX))
    --splu: LU decomposed matrix shape (NY*NX, NY*NX)

    Return: x = array[NY,NX]

    """
    # array 2D -> array 1D
    f2 = RHS.ravel()

    # Solving the linear system
    x = splu.solve(f2)

    return x.reshape(RHS.shape)


#########################################
# MAIN: Programme principal
#########################################


# Taille adimensionnee du domaine
# aspect_ratio = LY/LX

aspect_ratio = float(2.)
LY = float(1.)
LX = LY / aspect_ratio

# GRID RESOLUTION

# Taille des tableaux (points fantomes inclus)
NX = int(50)
NY = int(100)

# Allocate f and u
f = np.ones((NY, NX))
u = np.zeros((NY, NX))

# Elements differentiels
dx = LX / (NX - 1)
dy = LY / (NY - 1)

dx_2 = 1. / (dx * dx)
dy_2 = 1. / (dy * dy)


# Maillage
x = np.linspace(0, LX, NX)
y = np.linspace(0, LY, NY)
[xx, yy] = np.meshgrid(x, y)

f = ((np.pi / LX)**2 + (np.pi / LY)**2) * \
    np.sin(np.pi * xx / LX) * np.sin(np.pi * yy / LY)
analytic = - np.sin(np.pi * xx / LX) * np.sin(np.pi * yy / LY)
analytic = analytic.ravel()

# Matrix construction
LAPoisson = BuildLaPoisson()
LUPoisson = LUdecomposition(LAPoisson)

# Solve the linear system
u[1:-1, 1:-1] = ResoLap(LUPoisson, RHS=f[1:-1, 1:-1])

# FIGURE draw works only if plt.ion()
plt.clf()
plt.pcolor(xx, yy, u, shading='auto')
plt.colorbar()
plt.axis('image')
# plt.draw()
# plt.pause(1)
plt.show()


# exit(0)  # Stop here for the time being

# Convergence study
error = np.zeros(10)
h = np.zeros(10)

for k in range(10):
    NX = 10 * (k + 1)
    NY = NX
    f = np.zeros((NY, NX))
    u = np.zeros((NY, NX))
    dx = LX / (NX - 1)
    dy = LY / (NY - 1)
    dx_2 = 1. / (dx * dx)
    dy_2 = 1. / (dy * dy)
    x = np.linspace(0, LX, NX)
    y = np.linspace(0, LY, NY)
    [xx, yy] = np.meshgrid(x, y)
    LAPoisson = BuildLaPoisson()
    LUPoisson = LUdecomposition(LAPoisson)

    # Expression of the RHS
    f = ((np.pi / LX)**2 + (np.pi / LY)**2) * \
        np.sin(np.pi * xx / LX) * np.sin(np.pi * yy / LY)
    analytic = - np.sin(np.pi * xx / LX) * np.sin(np.pi * yy / LY)
    # analytic = analytic.ravel()

    # Solve the linear system
    u[1:-1, 1:-1] = ResoLap(LUPoisson, RHS=f[1:-1, 1:-1])

    error[k] = np.max(np.abs(analytic[:] - u[:]))
    h[k] = dx


plt.figure()
plt.loglog(h, error, 'ok-')
plt.loglog(h, (error[0] / h[0]) * h, 'k--')
plt.loglog(h, (error[0] / h[0]) * h**2, 'r--')
plt.loglog(h, (error[0] / h[0]) * h**4, 'g--')
plt.axis('tight')
plt.xlabel('dx')
plt.ylabel('error')
plt.legend(('error', '1st order', '2nd order', '4th order'))
plt.show()
# plt.draw()
# plt.pause(0.001)

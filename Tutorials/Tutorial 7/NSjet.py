# -*- coding:utf-8 -*-
#
#       Navier-Stokes flow in a driven cavity
#
#       Emmanuel Dormy (2021)
#

import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as lg


# affichage graphique
import matplotlib.pyplot as plt


def BuildLAPD():
    """
    Laplacian with Dirichlet BC

    """
    # Dropping ghost points (-2)
    NXi = nx
    NYi = ny

    # 1D Laplace operator

    # X-axis
    # Diagonal terms
    dataNXi = [np.ones(NXi), -2 * np.ones(NXi), np.ones(NXi)]

    # Boundary conditions
    dataNXi[1][0] = -3.  # Dirichlet
    dataNXi[1][-1] = -3.  # Dirichlet

    # Y-axis
    # Diagonal terms
    dataNYi = [np.ones(NYi), -2 * np.ones(NYi), np.ones(NYi)]

    # Boundary conditions
    dataNYi[1][0] = -3.  # Dirichlet
    dataNYi[1][-1] = -3.  # Dirichlet

    # Their positions
    offsets = np.array([-1, 0, 1])
    DXX = sp.dia_matrix((dataNXi, offsets), shape=(NXi, NXi)) * dx_2
    DYY = sp.dia_matrix((dataNYi, offsets), shape=(NYi, NYi)) * dy_2
    # print(DXX.todense())
    # print(DYY.todense())

    # 2D Laplace operator
    LAPD = sp.kron(sp.eye(NYi, NYi), DXX) + sp.kron(DYY, sp.eye(NXi, NXi))

    return LAPD


def BuildLAPN():
    """
    Laplacian matrix for Phi with Neumann BC

    The value is set at one point (here [0][1]) to ensure uniqueness

    """
    # Dropping ghost points (-2)
    NXi = nx
    NYi = ny

    # 1D Laplace operator

    # X-axis
    # Diagonal terms
    dataNXi = [np.ones(NXi), -2 * np.ones(NXi), np.ones(NXi)]

    # Boundary conditions
    dataNXi[1][0] = -1.  # Neumann
    dataNXi[1][-1] = -3.  # Dirichlet

    # Y-axis
    # Diagonal terms
    dataNYi = [np.ones(NYi), -2 * np.ones(NYi), np.ones(NYi)]

    # Boundary conditions
    dataNYi[1][0] = -1.  # Neumann
    dataNYi[1][-1] = -1.  # Neumann

    # Their positions
    offsets = np.array([-1, 0, 1])
    DXX = sp.dia_matrix((dataNXi, offsets), shape=(NXi, NXi)) * dx_2
    DYY = sp.dia_matrix((dataNYi, offsets), shape=(NYi, NYi)) * dy_2
    # print(DXX.todense())
    # print(DYY.todense())

    # 2D Laplace operator
    LAP = sp.kron(DXX, sp.eye(NYi, NYi)) + sp.kron(sp.eye(NXi, NXi), DYY)

    # BUILD CORRECTION MATRIX

    # Upper Diagonal terms
    dataNYNXi = [np.zeros(NYi * NXi)]
    offset = np.array([1])

    # Fix coef: 2+(-1) = 1 ==> Dirichlet at a single point
    dataNYNXi[0][1] = -1 * dx_2

    LAP0 = sp.dia_matrix((dataNYNXi, offset), shape=(NYi * NXi, NYi * NXi))

    # tmp = LAP + LAP0
    # print(LAP.todense())
    # print(LAP0.todense())
    # print(tmp.todense())

    return LAP + LAP0


def LUdecomposition(LAP):
    """
    return the Incomplete LU decomposition
    of a sparse matrix LAP
    """
    return lg.splu(LAP.tocsc(),)


def Resolve(splu, RHS):
    """
    solve the system

    SPLU * x = RHS

    Args:
    --RHS: 2D array((NY,NX))
    --splu: (Incomplete) LU decomposed matrix
            shape (NY*NX, NY*NX)

    Return: x = array[NY,NX]

    Rem1: taille matrice fonction des CL

    """
    # array 2D -> array 1D
    f2 = RHS.ravel()

    # Solving the linear system
    x = splu.solve(f2)

    return x.reshape(RHS.shape)


####
def Laplacien(x):
    """
    calcule le laplacien scalaire
    du champ scalaire x(i,j)

    pas de termes de bord car ghost points

    """
    rst = np.empty((NX, NY))

    coef0 = -2 * (dx_2 + dy_2)

    rst[1:-1, 1:-1] = ((x[2:, 1:-1] + x[:-2, 1:-1]) * dx_2 +
                       (x[1:-1, 2:] + x[1:-1, :-2]) * dy_2 +
                       (x[1:-1, 1:-1]) * coef0)
    return rst


def divergence(u, v):
    """
    divergence avec points fantomes
    ne jamais utiliser les valeurs au bord

    """
    tmp = np.empty((NX, NY))

    tmp[1:-1, 1:-1] = (
        (u[2:, 1:-1] - u[:-2, 1:-1]) / dx / 2 +
        (v[1:-1, 2:] - v[1:-1, :-2]) / dy / 2)

    return tmp


###
def VelocityGhostPoints(u, v):
    # left
    u[0, :] = -u[1, :]
    u[0, 30:60] = 2 - u[1, 30:60]
    v[0, :] = -v[1, :]
    # right
    u[-1, :] = u[-2, :]
    v[-1, :] = v[-2, :]
    # bottom
    u[:, 0] = u[:, 1]
    v[:, 0] = -v[:, 1]
    # top
    u[:, -1] = u[:, -2]
    v[:, -1] = -v[:, -2]


def PhiGhostPoints(phi):
    """
    copie les points fantomes
    tjrs Neumann

    global ==> pas de return

    """
    # left
    phi[0, :] = phi[1, :]
    # right
    phi[-1, :] = -phi[-2, :]
    # bottom
    phi[:, 0] = phi[:, 1]
    # top
    phi[:, -1] = phi[:, -2]


####

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
# MAIN: Programme principal
#########################################


# Domain Size
# aspect_ratio = LY/LX

LY = float(1.)
LX = float(2.)

# Grid Size

# (incuding ghost points)

NX = int(200)
NY = int(90)

# Number of iterations for the divergence
NI = 1


# Taille du domaine reel
nx = NX - 2
ny = NY - 2

# Nombre d'iterations
nitermax = int(10001)

# Modulo
modulo = int(1000)

# CONDITIONS INITIALES

# Valeurs initiales des vitesses
u = np.zeros((NX, NY))
v = np.zeros((NX, NY))

# Elements differentiels

dx = LX / (nx)
dy = LY / (ny)

dx_2 = 1. / (dx * dx)
dy_2 = 1. / (dy * dy)

# Grid, for plotting only
x = np.linspace(dx / 2, LX - dx / 2, nx)
y = np.linspace(dy / 2, LY - dy / 2, ny)
[xx, yy] = np.meshgrid(x, y)

# ATTENTION: dt_init calculer la CFL a chaque iteration...
dt = 0.0001

t = 0.  # total time

# parameters (Reynolds number)
Re = 100.0


# Tableaux avec points fantomes
# Matrices dans lesquelles se trouvent les extrapolations
ADVu = np.zeros((NX, NY))
ADVv = np.zeros((NX, NY))

# Definition des matrices ustar et vstar
ustar = np.zeros((NX, NY))
vstar = np.zeros((NX, NY))

# Definition de divstar
divstar = np.zeros((NX, NY))

# Definition de la pression phi
phi = np.zeros((NX, NY))
gradphix = np.zeros((NX, NY))
gradphiy = np.zeros((NX, NY))


# CONSTRUCTION des matrices et LU decomposition

# Matrix construction for projection step
LAPN = BuildLAPN()
LUPN = LUdecomposition(LAPN)


################
# MAIN LOOP
tStart = t

for niter in range(nitermax):

    t += dt

# for k in range(NI):

    # Advection semi-Lagrangienne
    ADVu = Semilag(u, v, u)
    ADVv = Semilag(u, v, v)

    # Diffusion step
    ustar = ADVu + dt * Laplacien(u) / Re
    vstar = ADVv + dt * Laplacien(v) / Re

    # Ghost points update
    VelocityGhostPoints(ustar, vstar)

    # Update divstar
    divstar = divergence(ustar, vstar)
#    divstar = divstar - np.mean(divstar[1:-1,1:-1])

    # Solving the linear system
    phi[1:-1, 1:-1] = Resolve(LUPN, RHS=divstar[1:-1, 1:-1])

    # update Pressure ghost points
    PhiGhostPoints(phi)

    # Update gradphi
    gradphix[1:-1, :] = (phi[2:, :] - phi[:-2, :]) / dx / 2
    gradphiy[:, 1:-1] = (phi[:, 2:] - phi[:, :-2]) / dy / 2

    # Project u
    u = ustar - gradphix
    v = vstar - gradphiy

    # Mise a jour des points fantomes
    # pour le champ de vitesse et T

    VelocityGhostPoints(u, v)

    if ((niter + 1) % modulo == 0):

        # logfile
        sys.stdout.write(
            '\niteration: %d -- %i %%\n'
            '\n'
            'total time     = %.2e\n'
            '\n'
            % (niter,
               float(niter) / nitermax * 100,
               t))

        # FIGURE draw works only if plt.ion()
        plotlabel = "t = %1.5f" % (t)
        plt.clf()
        plt.title(plotlabel)
        plt.quiver(xx[::4, ::4], yy[::4, ::4], np.transpose(
            u[1:-1:4, 1:-1:4]), np.transpose(v[1:-1:4, 1:-1:4]), 0.9)
        plt.axis('image')
        plt.draw()
        plt.pause(0.1)

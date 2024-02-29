# -*- coding:utf-8 -*-
#
# Periodic 2D flow using a spectral method
# Emmanuel Dormy
# Introduction to CFD, Cambridge 2022
#


import numpy as np
# from numpy.fft import ifft2, fft2
from scipy.fftpack import fft2, ifft2


# Plotting
import matplotlib.pyplot as plt


def phys_spec(UR, NX, NY, NXSp, NYSp):
    #  Physical to Spectral transform

    # Total nb of points
    NN = NX * NY

# Fourier transform
    UU = fft2(UR / NN)

# Only keep the relevant modes
    NXSp_2 = np.int32(NXSp / 2)
    NYSp_2 = np.int32(NYSp / 2)
    Unew = 1j * np.zeros([NXSp, NYSp])
    Unew[0:NXSp_2, 0:NYSp_2] = UU[0:NXSp_2, 0:NYSp_2]
    Unew[NXSp_2:NXSp, 0:NYSp_2] = UU[NX - NXSp_2:NX, 0:NYSp_2]
    Unew[0:NXSp_2, NYSp_2:NYSp] = UU[0:NXSp_2, NY - NYSp_2:NY]
    Unew[NXSp_2:NXSp, NYSp_2:NYSp] = UU[NX - NXSp_2:NX, NY - NYSp_2:NY]

    return Unew


def spec_phys(U, NX, NY, NXSp, NYSp):
    #  Spectral to Physical transform

    # Total nb of points
    NN = NX * NY

# Create a large array
    UU = 1j * np.zeros([NX, NY])

# Fill with relevant modes
    NXSp_2 = np.int32(NXSp / 2)
    NYSp_2 = np.int32(NYSp / 2)

    UU[0:NXSp_2, 0:NYSp_2] = U[0:NXSp_2, 0:NYSp_2]
    UU[NX - NXSp_2:NX, 0:NYSp_2] = U[NXSp_2:NXSp, 0:NYSp_2]
    UU[0:NXSp_2, NY - NYSp_2:NY] = U[0:NXSp_2, NYSp_2:NYSp]
    UU[NX - NXSp_2:NX, NY - NYSp_2:NY] = U[NXSp_2:NXSp, NYSp_2:NYSp]

# Inverse Fourier transform
    Unew = NN * np.real(ifft2(UU))

    return Unew


# Reynolds number
Re = 1.0e3

# Resolution in physical space
NX = 128
NY = 128

# Resolution in spectral space (with de-aliasing)
NXSp = 2 * np.int32(NX / 3)
NYSp = 2 * np.int32(NY / 3)

dt = 5e-2     # time step
NT = 1000      # number of iterations
modulo = 100   # modulo for plotting


ic = 1j

# Mesh generation
dx = 2 * np.pi / NX
dy = 2 * np.pi / NY
t = 0.

x = np.linspace(1, NX + 1, NX + 1)
y = np.linspace(1, NY + 1, NY + 1)
yy, xx = np.meshgrid(y, x)
XP = xx * dx
YP = yy * dy
X = XP[:-1, :-1]
Y = YP[:-1, :-1]


# Vortex initialisation
w = np.zeros([NX, NY])
# w = np.exp(-((X-np.pi)**2))
w = np.exp(-((X - np.pi)**2 + (Y - np.pi + np.pi / 4)**2) / (0.2)) \
    + np.exp(-((X - np.pi)**2 + (Y - np.pi - np.pi / 4)**2) / (0.2)) \
    - 0.5 * np.exp(-((X - np.pi - np.pi / 4)**2 +
                   (Y - np.pi - np.pi / 4)**2) / (0.4))

w_init = w.copy()

# 2D arrays of wavenumbers
# in the x direction
kx = ic * (np.outer((np.mod(np.linspace(1, NXSp, NXSp) -
           np.ceil(NXSp / 2 + 1), NXSp) - np.floor(NXSp / 2)), np.ones(NYSp)))

# in the y direction
ky = ic * np.outer(np.ones(NXSp),
                   (np.mod(np.linspace(1,
                                       NYSp,
                                       NYSp) - np.ceil(NYSp / 2 + 1),
                           NYSp) - np.floor(NYSp / 2)))


# Laplacian in Fourier space
k2 = kx**2 + ky**2
k2_fixed = k2

# "fixed" Laplacian for the Poisson equation
k2_fixed[0, 0] = 1

w_hat = phys_spec(w, NX, NY, NXSp, NYSp)

psi_hat = np.zeros([NXSp, NYSp])

for i in range(NT):

    # Time-stepping
    # Dissipation (with crank-nicolson)
    C = - k2 / Re * dt / 2
    w_hat = (1 + C) / (1 - C) * w_hat

    psi_hat = w_hat / k2_fixed

    # get u from psi
    u = spec_phys(ic * ky * psi_hat, NX, NY, NXSp, NYSp)
    v = spec_phys(-ic * kx * psi_hat, NX, NY, NXSp, NYSp)

    diffX_w = spec_phys(ic * kx * w_hat, NX, NY, NXSp, NYSp)
    diffY_w = spec_phys(ic * ky * w_hat, NX, NY, NXSp, NYSp)

    # advection term
    adve = u * diffX_w + v * diffY_w
    adve_hat = phys_spec(adve, NX, NY, NXSp, NYSp)

    w_hat = w_hat - adve_hat * dt

    if (i % modulo == 0):
        # Plot of the vorticity field
        w = spec_phys(w_hat, NX, NY, NXSp, NYSp)
        plt.figure(1)
        plt.set_cmap('jet')
        plt.clf()
        t = i * dt
        plotlabel = "t = %1.5f" % (t)
        plt.title(plotlabel)
        plt.pcolormesh(XP, YP, w, shading='flat')
        plt.colorbar()
        plt.axis('image')
        plt.draw()
        plt.pause(0.1)

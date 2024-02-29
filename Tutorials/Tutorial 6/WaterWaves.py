# -*- coding:utf-8 -*-
#
#       1D Water waves at arbirarily depth (spectral)
#
#       Emmanuel Dormy (2021)
#
   
import numpy as np
import matplotlib.pyplot as plt


g=9.8
h0=0.001
#h0=1.0
#h0=1000.0

NSp=256
NX=2*NSp

NT=100000

dt=1e-3

def phys_spec(UR,NX,NSp):
   Unew=np.fft.rfft(UR/NX,NX)
   Unew.resize(NSp)
   return Unew

def spec_phys(U,NX,NSp):
   Unew=NX*np.fft.irfft(U,NX)
   return Unew



xtmp=np.linspace(0,2*np.pi,NX+1)

x=xtmp[0:NX]

eta=np.exp(-20*(x-np.pi)**2)
psi=0*x

eta_hat=phys_spec(eta,NX,NSp)
psi_hat=phys_spec(psi,NX,NSp)


k=np.linspace(0,NSp-1,NSp)
# Deep water
eig=k
# Arb. depth
### TO BE MODIFIED

for it in range(NT):

  eta_hat=eta_hat+dt*eig*psi_hat
  psi_hat=psi_hat-dt*g*eta_hat

    # Visualisation
  if (it%100==0):
      plt.figure(1)
      eta=spec_phys(eta_hat,NX,NSp)
      plt.clf()
      plt.plot(x,eta)
      plt.title("{:10.2f}".format(it*dt))
#      plt.title(it*dt)
      plt.axis([0, 2*np.pi, -1, 1])
      plt.draw()
      plt.pause(0.001)

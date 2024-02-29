# -*- coding:utf-8 -*-
#
# Precision of a Chebychev discretisation
# Emmanuel Dormy
# Introduction to CFD, Cambridge 2022
#

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as lg
from scipy.sparse import csr_matrix

###### Plotting
import matplotlib.pyplot as plt
plt.ion()
                 

def CDmat(n):
    """
    This routine computes compact differentiation matrices
    """


    dataA = [-3*np.ones(n), 0*np.ones(n), 3*np.ones(n)]
# suppress 1st and last lines
    dataA[2][1]  = 0.
    dataA[1][0]  = 0.
    dataA[1][-1] = 0.
    dataA[0][-2] = 0.
    
    offsets = np.array([-1,0,1])                    
    A = sp.dia_matrix((dataA,offsets), shape=(n,n))
    A_cl=np.zeros([n,n])
    A_cl[0,0:4]=[-17/6, 3/2, 3/2, -1/6]
    A_cl[-1,-4:]=[1/6, -3/2, -3/2, 17/6]

    A=A+sp.csr_matrix(A_cl)

    dataB = [np.ones(n), 4*np.ones(n), np.ones(n)]
# modify 1st and last lines
    dataB[2][1]  = 3.
    dataB[1][0]  = 1.
    dataB[1][-1] = 1.
    dataB[0][-2] = 3.
    
    offsets = np.array([-1,0,1])                    
    B = sp.dia_matrix((dataB,offsets), shape=(n,n))

    return A,B

                

def Cmat(n):
    """
    This routine computes Chebychev differentiation matrices

 N       = number of modes
 D0      = zero'th derivative matrix
 D1      = first derivative matrix
 D2      = second derivative matrix
 D4      = fourth derivative matrix
    """


# initialize

# create D0

    vec=np.linspace(0,n,n+1)
    D0=np.cos(np.outer(vec,vec)*np.pi/n);

# create higher derivative matrices

    lv=np.size(vec) # lv = n+1;
    
    D1=np.zeros([lv,lv])
    D2=np.zeros([lv,lv])
    D3=np.zeros([lv,lv])
    D4=np.zeros([lv,lv])
    D1[:,1]=D0[:,0]
    D1[:,2]=4*D0[:,1]
    D2[:,2]=4*D0[:,0]
    for i in range(n-2): # We compute Tj, which is in column j+1
      j=i+2
      jj=j+1
      D1[:,j+1]=2*jj*D0[:,j]+jj*D1[:,j-1]/(jj-2);
      D2[:,j+1]=2*jj*D1[:,j]+jj*D2[:,j-1]/(jj-2);
      D3[:,j+1]=2*jj*D2[:,j]+jj*D3[:,j-1]/(jj-2);
      D4[:,j+1]=2*jj*D3[:,j]+jj*D4[:,j-1]/(jj-2);
    return D0,D1,D2,D4



#########################################
###### MAIN: Programme principal
#########################################



kmax=20
   
# CD approach
nnn=np.zeros(kmax)
err_CD=np.zeros(kmax)
   
for k in range(kmax):
    n=5*(k+1)
       
    nnn[k]=n

# Compact Finite Differences

    dx = 2/(n-1)
    x  = np.linspace(-1,1,n)

    u=np.sin(np.pi*x)
    du=np.pi*np.cos(np.pi*x)
    d2u=-np.pi**2*np.sin(np.pi*x)

    CDmats=CDmat(n)
    D1_A=CDmats[0]
    D1_B=CDmats[1]


    du_CD=lg.spsolve(dx*D1_B,D1_A.dot(u))
    d2u_CD=lg.spsolve(dx*D1_B,D1_A.dot(du_CD))    
    err_CD[k]=np.max(abs(d2u-d2u_CD))

plt.figure(1)
plt.loglog(nnn,err_CD,'o')
plt.grid(which='both')    
plt.xlabel('Number of points')
plt.ylabel('Error')
plt.title('Compact Scheme')


# Spectral approaches (Chebyshev)

nnn2=np.zeros(kmax)
err_Cheb=np.zeros(kmax)

for k in range(kmax):
    n=5+k
    nnn2[k]=n
       
# Chebyshev
    nn=n-1
    vec=np.linspace(0,nn,nn+1)
    x_Cheb=np.cos(np.pi*vec/nn)

    u=np.sin(np.pi*x_Cheb)
    du=np.pi*np.cos(np.pi*x_Cheb)
    d2u=-np.pi**2*np.sin(np.pi*x_Cheb)
    CHEBmats=Cmat(nn)
    D0=CHEBmats[0]
    D1=CHEBmats[1]
    D2=CHEBmats[2]
    D4=CHEBmats[3]

    au=np.linalg.solve(D0,u)
    d2u_Cheb=D2.dot(au)
    err_Cheb[k]=np.max(abs(d2u-d2u_Cheb))

plt.figure(2)
plt.loglog(nnn,err_Cheb,'o')
plt.grid(which='both')    
plt.xlabel('Number of points')
plt.ylabel('Error')
plt.title('Chebyshev')

plt.draw()

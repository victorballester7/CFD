# -*- coding:utf-8 -*-
#
#       Stokes flow in a driven cavity
#       (using an iterative pressure correction code)
#
#       Emmanuel Dormy (2021)
#

import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as lg


###### Plotting the results
import matplotlib.pyplot as plt
plt.ion()
                 


def BuildLAPD():
    """
    Laplacian with Dirichlet BC

    """
    ### Dropping ghost points (-2)
    NXi = nx
    NYi = ny

    ###### 1D Laplace operator

    ###### X-axis
    ### Diagonal terms
    dataNXi = [np.ones(NXi), -2*np.ones(NXi), np.ones(NXi)]   
    
    ### Boundary conditions
    dataNXi[1][0]  = -3.  # Dirichlet
    dataNXi[1][-1] = -3.  # Dirichlet

    ###### Y-axis
    ### Diagonal terms
    dataNYi = [np.ones(NYi), -2*np.ones(NYi), np.ones(NYi)] 
   
    ### Boundary conditions
    dataNYi[1][0]  = -3.  # Dirichlet
    dataNYi[1][-1] = -3.  # Dirichlet

    ###### Their positions
    offsets = np.array([-1,0,1])                    
    DXX = sp.dia_matrix((dataNXi,offsets), shape=(NXi,NXi)) * dx_2
    DYY = sp.dia_matrix((dataNYi,offsets), shape=(NYi,NYi)) * dy_2
    #print(DXX.todense())
    #print(DYY.todense())
    
    ####### 2D Laplace operator (Dirichlet)
    LAPD = sp.kron(DXX, sp.eye(NYi,NYi)) + sp.kron(sp.eye(NXi,NXi), DYY)
    
    return LAPD 

def BuildLAPN():
    """
    Laplacian matrix for Phi with Neumann BC

    The value is set at one point (here [0][1]) to ensure uniqueness
    
    """
    ### Dropping ghost points (-2)
    NXi = nx
    NYi = ny

    ###### 1D Laplace operator

    ###### X-axis
    ### Diagonal terms
    dataNXi = [np.ones(NXi), -2*np.ones(NXi), np.ones(NXi)]   
    
    ### Boundary conditions
    dataNXi[1][0]  = -1.  # Neumann
    dataNXi[1][-1] = -1.  # Neumann

    ###### Y-axis
    ### Diagonal terms
    dataNYi = [np.ones(NYi), -2*np.ones(NYi), np.ones(NYi)] 
   
    ### Boundary conditions
    dataNYi[1][0]  = -1.  # Neumann
    dataNYi[1][-1] = -1.  # Neumann

    ###### Their positions
    offsets = np.array([-1,0,1])                    
    DXX = sp.dia_matrix((dataNXi,offsets), shape=(NXi,NXi)) * dx_2
    DYY = sp.dia_matrix((dataNYi,offsets), shape=(NYi,NYi)) * dy_2
    #print(DXX.todense())
    #print(DYY.todense())
    
    ####### 2D Laplace operator
    LAP = sp.kron(DXX,sp.eye(NYi,NYi)) + sp.kron(sp.eye(NXi,NXi),DYY)
    
    ####### BUILD CORRECTION MATRIX

    ### Upper Diagonal terms
    dataNYNXi = [np.zeros(NYi*NXi)]
    offset = np.array([1])

    ### Fix coef: 2+(-1) = 1 ==> Dirichlet at a single point
    dataNYNXi[0][1] = -1 * dy_2

    LAP0 = sp.dia_matrix((dataNYNXi,offset), shape=(NYi*NXi,NYi*NXi))
  
    return LAP + LAP0



def LUdecomposition(LAP):
    """
    return the Incomplete LU decomposition 
    of a sparse matrix LAP
    """
    return  lg.splu(LAP.tocsc(),)


def Resolve(splu,RHS):
    """
    solves the system

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
def Laplacian(x):
    """
    computes the Laplacian of the 
    scalar field x(i,j)
    
    no boundary terms because of ghost points

    """
    rst = np.empty((NX,NY))

    coef0 = -2*(dx_2 + dy_2) 
    

    rst[1:-1,1:-1] = ( (x[2:, 1:-1] + x[:-2, 1:-1])*dx_2 +
                       (x[1:-1, 2:] + x[1:-1, :-2])*dy_2 +  
                       (x[1:-1, 1:-1])*coef0 )    
    return rst

def divergence(u,v):
    """
    computes the divergence of the 
    vector field (u,v)

    """
    tmp = np.empty((NX,NY))
    
    tmp[1:-1,1:-1] = (
        (u[2:, 1:-1] - u[:-2, 1:-1])/dx/2 +
        (v[1:-1, 2:] - v[1:-1, :-2])/dy/2 )
        
    return tmp


###
def VelocityGhostPoints(u,v):
    """
        No-slip (Dirichlet) boundary conditions

        global ==> no need for return
    """
    ### bottom
    u[:,  0] = -u[:,  1] 
    v[:,  0] = -v[:,  1] 
    ### top     
    u[:, -1] = 2-u[:, -2] # Try to understand that
    v[:, -1] = -v[:, -2] 
    ### right     
    u[0,  :] = -u[1,  :] 
    v[0,  :] = -v[1,  :] 
    ### left    
    u[-1, :] = -u[-2, :] 
    v[-1, :] = -v[-2, :] 

        
def PhiGhostPoints(phi):
    """
    Neumann boundary conditions

    global ==> no need for return 

    """
    ### bottom               
    phi[:,  0] = phi[:,  1]
    ### top            
    phi[:, -1] = phi[:, -2]
    ### right
    phi[0,  :] = phi[1,  :]
    ### left             
    phi[-1, :] = phi[-2, :]



#########################################
######         MAIN program        ######
#########################################


###### Geometry of the Domain ######

### aspect_ratio = LY/LX  

aspect_ratio = float(1.)
LY = float(1.)
LX = LY/aspect_ratio

###### Grid Size

### (incuding ghost points)

NX = int(30)
NY = int(32)

# Number of iterations for the divergence
NI=1

### Real size of the domain
nx = NX-2
ny = NY-2

### Number of iterations
nitermax = int(10001)

### Modulo
modulo = int(1000)



###### INITIAL CONDITIONS #######

##### Initial velocities
u = np.zeros((NX,NY)) 
v = np.zeros((NX,NY))

##### Differential elements 
dx = LX/(nx)
dy = LY/(ny)

dx_2 = 1./(dx*dx)
dy_2 = 1./(dy*dy)

### Grid, for plotting only
x = np.linspace(dx/2,LX-dx/2,nx) 
y = np.linspace(dy/2,LY-dy/2,ny)
[xx,yy] = np.meshgrid(x,y) 

### dt (could be recomputed using the CFL condition) 
dt = 0.0001

t = 0. # time



###### Arrays with ghost points

### matrices ustar and vstar
ustar = np.zeros((NX,NY))
vstar = np.zeros((NX,NY))

### divergence array
divstar = np.zeros((NX,NY))

### Pressure phi
phi      = np.zeros((NX,NY))
gradphix = np.zeros((NX,NY))
gradphiy = np.zeros((NX,NY))


###### Building the system matrices and LU decomposition

### Matrix construction for projection step
LAPN = BuildLAPN() 
LUPN = LUdecomposition(LAPN)



#####################
###### MAIN TIME LOOP
##################### 
tStart = t

for niter in range(nitermax):
        
    t += dt

    ###### Diffusion step

    ustar = u + dt*Laplacian(u) 
    vstar = v + dt*Laplacian(v) 
        
    ###### Ghost points update
    VelocityGhostPoints(ustar,vstar)

    ### Update divstar 
    divstar = divergence(ustar,vstar)
    divstar = divstar - np.mean(divstar[1:-1,1:-1])
    
    ### Solving the linear system
    phi[1:-1,1:-1] = Resolve(LUPN, RHS=divstar[1:-1,1:-1])

    ### update Pressure ghost points 
    PhiGhostPoints(phi)

    ### Update gradphi
    gradphix[1:-1, :] = (phi[2:, :] - phi[:-2, :])/dx/2
    gradphiy[:, 1:-1] = (phi[:, 2:] - phi[:, :-2])/dy/2

    ### Project u
    u = ustar - gradphix
    v = vstar - gradphiy

    ###### Updating velocity ghost points

    VelocityGhostPoints(u,v)

    # Plot
    if ((niter+1)%modulo==0):

        ###### logfile
        sys.stdout.write(
            '\niteration: %d -- %i %%\n'
            '\n'
            'total time     = %.2e\n'
            '\n'
            %(niter,                    
              float(niter)/nitermax*100,
              t))
        
        ###### FIGURE draw works only if plt.ion()
        plotlabel = "t = %1.5f" %(t)
        plt.clf()
        plt.title(plotlabel)
        plt.quiver(xx,yy,np.transpose(u[1:-1,1:-1]),np.transpose(v[1:-1,1:-1]),4)
        plt.axis('image')
        plt.draw()
        plt.pause(0.001)
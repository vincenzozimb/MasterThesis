# -*- coding: utf-8 -*-
""" 
mainTRG.py
---------------------------------------------------------------------
Script file for initializing a partition function (here the classical \
Ising model) before passing to a Tensor Renormalization Group (TRG) \
coarse-graining routine. The magnetization, free energy and internal \
energy are then computed and compared against the exact results.

by Glen Evenbly (c) for www.tensors.net, (v1.1) - last modified 26/1/2019
"""

#### Preamble
import numpy as np
from numpy import linalg as LA

from ncon import ncon
from doTRG import doTRG 

##### Example 1: square lattice classical Ising model #####
###########################################################

##### Set bond dimensions and temperature
chiM = 32
chiH = 32
chiV = 32
relTemp = 0.98 # temp relative to the crit temp, in [0,inf]
numlevels = 16 # number of coarse-grainings
OPTS_midsteps = 10 # iterations in isometry optimization

##### Define partition function (classical Ising)
Tc = (2/np.log(1+np.sqrt(2)))
Tval  = relTemp*Tc
betaval = 1/Tval
Jtemp = np.zeros((2,2,2)); Jtemp[0,0,0]=1; Jtemp[1,1,1]=1
Etemp = np.array([[np.exp(betaval),np.exp(-betaval)],[np.exp(-betaval),np.exp(betaval)]])

Ainit = ncon([Jtemp,Jtemp,Jtemp,Jtemp,Etemp,Etemp,Etemp,Etemp],
             [[-1,8,1],[-2,2,3],[-3,4,5],[-4,6,7],[1,2],[3,4],[5,6],[7,8]])
Xloc = (1/np.sqrt(2))*np.array([[1,1],[1,-1]])
Ainit = ncon([Ainit,Xloc,Xloc,Xloc,Xloc],[[1,2,3,4],[1,-1],[2,-2],[3,-3],[4,-4]])
Atemp = ncon([Ainit,Ainit,Ainit,Ainit],[[-2,-3,3,1],[3,-4,-5,2],[-1,1,4,-8],[4,2,-6,-7]]).reshape(4,4,4,4)
Anorm = np.ones(numlevels+1)
Anorm[0] = LA.norm(Atemp.flatten())
A = [0 for x in range(numlevels+1)]
A[0] = Atemp/Anorm[0]

##### Do iterations of TRG
SPerrs = np.zeros((numlevels,3))
qC = [0 for x in range(numlevels)]
vC = [0 for x in range(numlevels)]
wC = [0 for x in range(numlevels)]
for k in range(numlevels):
    A[k+1], qC[k], vC[k], wC[k], Anorm[k+1], SPerrs[k,:] = doTRG(A[k],chiM,chiH,chiV, midsteps = OPTS_midsteps)
    print('RGstep: %d, Truncation Errors: %e, %e, %e' % (k,SPerrs[k,0],SPerrs[k,1],SPerrs[k,2]))


"""
evalTRG: evaluates outputs from TNR algorithm to compute several \
quantities for the Ising model (free energy, internal energy)
"""
def evalTRG(numlevels, Tval, A, vC, wC, Anorm):

    ##### Evaluate free energy density as a function of system size
    NumSpins = 2**(2*np.int64(np.array(range(1,21)))+1)
    FreeEnergy = np.zeros(numlevels);
    for k in range(1,numlevels+1): 
        Hgauge = ncon([vC[k-1],vC[k-1]],[[1,2,-1],[2,1,-2]])
        Vgauge = ncon([wC[k-1],wC[k-1]],[[1,2,-1],[2,1,-2]])
        FreeEnergy[k-1] = -Tval*(sum((4**np.int64(np.array(range(k,-1,-1))))*np.log(Anorm[:(k+1)])) +
                  np.log(ncon([A[k],Hgauge,Vgauge],[[1,3,2,4],[1,2],[3,4]]))) / NumSpins[k]

    ##### Evaluate local density matrix, internal energy, spontaneous mag
    maglambda = 1/(np.sinh(2/Tval)**2) # mag field (quantum Ising) from temp (classical Ising)
    rhoone = [0 for x in range(numlevels+1)]
    rhoone[numlevels] = np.eye(wC[numlevels-1].shape[2]) / wC[numlevels-1].shape[2]
    for k in range(numlevels-1,-1,-1):
        rhoone[k] = ncon([rhoone[k+1],wC[k],wC[k]],[[1,2],[-1,3,1],[-2,3,2]])
    
    sX = np.array([[0, 1], [1, 0]])
    sY = np.array([[0, -1j], [1j, 0]])
    sZ = np.array([[1, 0], [0,-1]])
    sI = np.array([[1, 0], [0,1]])
    hloc = (-np.kron(sX,sX) - (maglambda/2)*(np.kron(sZ,sI) + np.kron(sI,sZ))).reshape(4,4)
    IntEnergy = np.trace(hloc @ rhoone[0])
    
    sXcg = [0 for x in range(numlevels+1)]
    sXcg[0] = np.kron(sX,sI)
    for k in range(numlevels):
        sXcg[k+1] = ncon([sXcg[k],wC[k],wC[k]],[[1,2],[1,3,-1],[2,3,-2]])
    
    dtemp = LA.eigvalsh(sXcg[numlevels])
    ExpectX = max(abs(dtemp))
    
    ##### Compare with exact results (thermodynamic limit)
    N = 1000000;
    x = np.linspace(0,np.pi,N+1)
    y = np.log(np.cosh(2*betaval)*np.cosh(2*betaval) + (1/maglambda)*np.sqrt(
            1+maglambda**2 -2*maglambda*np.cos(2*x)))
    FreeExact = -Tval*((np.log(2)/2) + 0.25*sum(y[1:(N+1)] + y[:N])/N)
    RelFreeErr = abs((FreeEnergy[numlevels-1]-FreeExact)/FreeExact)
    
    x = np.linspace(0,2*np.pi,N+1)
    y = np.sqrt((maglambda-1)**2 +4*maglambda*np.sin(x/2)**2)
    IntExact = -0.5*sum(y[1:(N+1)] +  y[:N])/N
    RelIntErr = abs((IntEnergy-IntExact)/IntExact)
    
    ExpectXExact = 0
    if maglambda < 1:
        ExpectXExact = (1 - np.sinh(2*betaval)**(-4))**(1/8)
    
    return RelFreeErr, RelIntErr, ExpectX, ExpectXExact

#----------------------------------------------------
# Evaluate obserables of interest
RelFreeErr, RelIntErr, ExpectX, ExpectXExact = evalTRG(numlevels, Tval, A, vC, wC, Anorm)
print('----------------- Final Results -----------------')
print('Mag: %f, ExactMag: %f, Err_FreeEn: %e, Err_IntEn: %e, RGsteps: %d, Bond dims: %d, %d, %d' % (
        ExpectX,ExpectXExact,RelFreeErr,RelIntErr,numlevels,chiM,chiV,chiH))
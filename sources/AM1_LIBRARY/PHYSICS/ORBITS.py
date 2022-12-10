# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 18:29:45 2022

@author: serg_
"""
from AM1_LIBRARY.NUMERIC.NewtonSolve import Jac, newton
from numpy import array, sqrt, zeros
from numpy.linalg import norm, eig
from scipy.optimize import newton, fsolve

def Kepler(U, t): 

    x = U[0]; y = U[1]; dxdt = U[2]; dydt = U[3]
    d = ( x**2  +y**2 )**1.5

    return  array( [ dxdt, dydt, -x/d, -y/d ] ) 

def CR3BP(U, t, mu): #Circular Restricted 3 Body Problem
    
    x = U[0]; y = U[1]
    vx = U[2]; vy = U[3]
    
    r1 = sqrt( (x+mu)**2 + y**2 )
    r2 = sqrt( (x-1+mu)**2 + y**2 )
    
    dxdt = vx
    dydt = vy
    
    dvxdt = 2*vy+x-( (1-mu)*(x+mu) ) / (r1**3) - mu*(x+mu-1)/(r2**3)
    dvydt = -2*vx + y -( (1-mu) / (r1**3) + mu/(r2**3) )*y
    
    return array([ dxdt, dydt, dvxdt, dvydt])

def LG_Points( U0 ,N ,mu ):
    
    LGP = zeros([5,2])

    def F(Y):
        X = zeros(4)
        X[0:2] = Y
        X[2:4] = 0
        FX = CR3BP(X, 0 , mu)
        return FX[2:4]
   
    for i in range(N):
        LGP[i,:] = fsolve(F, U0[i,0:2])
        
    return LGP

def LP_Stab( U0, mu ):

    def F(Y):
        return CR3BP(Y, 0 , mu)

    A = Jac(F, U0)
    values, vectors = eig(A)

    return values
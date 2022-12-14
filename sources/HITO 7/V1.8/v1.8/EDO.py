# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 19:52:48 2022
@author: serg_
"""

from numpy import zeros, float64

"""
"""

def CauchyProblem( F, t, U0, TempScheme): 
    
    #N means the number of time nodes and Nv to time & velocity number of coordinates

     N, Nv=  len(t)-1, len(U0) 
     U = zeros( (N+1, Nv), dtype=float64 ) 
     U[0,:] = U0

     for n in range(N): #el ultimo es N-2, ya que calculo en N-1 que es el ultimo, al empezar en 0!!
        U[n+1,:] = TempScheme( U[n, :], t[n+1] - t[n], t[n],  F ) 

     return U

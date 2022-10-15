# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 19:39:39 2022

@author: serg_
"""
from scipy.optimize import fsolve

"""
----------------------------------------------------------
TEMP_SCHEMES

This Module includes different temp schemes to integrate numerical problems

-Inputs:
    
    U: State Vector at integration time tn
    dt: time step
    t: integration time instant tn
    F(U,t): Function dU/dt = F(U,t)
    
-Output: the return
    
    U State Vector at desired step tn + dt
----------------------------------------------------------
"""

def Euler_Scheme(U, dt, t, F):
        
    return U + dt * F(U, t)



def RK4_Scheme(U, dt, t, F):

    k1 = F(U,t)
    k2 = F(U + k1 * dt/2, t + dt/2)
    k3 = F(U + k2 * dt/2, t + dt/2)
    k4 = F(U + k3 * dt, t + dt)

    return U + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)



def InvEuler_Scheme(U, dt, t, F):
    
    def InvEuler_Eq(X):     
          return X - U - dt * F(X, t)

    return fsolve(InvEuler_Eq, U ) 



def CrankNicolson_Scheme(U, dt, t, F ): 

    def CN_Eq(X): 
         
         return  X - U - (F(X,t+dt) + F(U, t))*dt/2
  
    return fsolve( CN_Eq, U )




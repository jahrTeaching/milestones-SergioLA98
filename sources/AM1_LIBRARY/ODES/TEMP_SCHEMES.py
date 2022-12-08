# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 19:39:39 2022

@author: serg_
"""
from AM1_LIBRARY.NUMERIC.NewtonSolve import newton

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
    Euler_Scheme.__name__ = "Euler"    
    return U + dt * F(U, t)



def RK4_Scheme(U, dt, t, F):

    RK4_Scheme.__name__ = "Runge-Kutta 4"
    k1 = F(U,t)
    k2 = F(U + k1 * dt/2, t + dt/2)
    k3 = F(U + k2 * dt/2, t + dt/2)
    k4 = F(U + k3 * dt, t + dt)

    return U + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)



def InvEuler_Scheme(U, dt, t, F):
    
    InvEuler_Scheme.__name__ = "Inverse Euler"
    def InvEuler_Eq(X):     
          return X - U - dt * F(X, t)

    return newton(InvEuler_Eq, U ) 



def CrankNicolson_Scheme(U, dt, t, F ): 

    CrankNicolson_Scheme.__name__ = "Crank Nicolson"
    def CN_Eq(X): 
         
         return  X - U - (F(X,t+dt) + F(U, t))*dt/2
  
    return newton( CN_Eq, U )


def LeapFrog (U, dt, t, F ):
     
     LeapFrog.__name__ = "LeapFrog"
     
     if t == 0:
         U = U + dt*F(U, t) #Euler Inizialitation
     else:
        p = int(len(U)/2)
        
        U_aux = U #puntero mismo ID
        aux = F (U_aux, t)

        U_aux [p :] += aux[p :]*dt/2.0 #no cambia ID, actualiza
        U_aux [: p] += aux[: p]*dt

        aux = F (U_aux, t)
        
        U_aux [p :] +=  aux[p :]*dt/2.0
        
        #U = U_aux innecesario, Uaux tiene el mismo id, es un alias
        # se va actualizando
        
     return U



from scipy.optimize import newton
# from Functions.Physics import Kepler

""""
Functions for the temporal schemes
    Inputs:
        U = vector value at tn
        dt = time step
        F = function dU/dt
        t = tn
    
    Output:
        U = vector value at (tn+dt)
"""""

def Euler(U, dt, t, F_E, Nb, M): 

    return (U + dt*F_E(U,t, Nb, M))

def Crank_Nicolson(U, dt, t, F, Nb, M): 
    def f_CN(X): 
         
         return  X - (U + dt/2*F(U, t, Nb, M)) - dt/2*F(X, t+dt, Nb, M)
 
    return newton(f_CN, U)

def RK4(U, dt, t, F, Nb, M):

    k1 = F(U, t, Nb, M)
    k2 = F(U + dt*k1/2, t + dt/2, Nb, M)
    k3 = F(U + dt*k2/2, t + dt/2, Nb, M)
    k4 = F(U + dt*k3, t + dt, Nb, M)
    
    return U + dt*(k1 + 2*k2 + 2*k3 + k4)/6

def Inverse_Euler(U, dt, t, F, Nb, M):
    def f_I(X):
        return X - U - dt*F(X,t, Nb, M)

    return newton(f_I, U)

def LeapFrog(U, dt, t, F, Nb, M):
    p = int(len(U)/2)

    V = U
    A = F(V,t, Nb, M)

    V[p:] += A[p:]*dt/2.0
    V[:p] += A[:p]*dt

    A = F(V,t, Nb, M)

    V[p:] += A[p:]*dt/2.0

    return V

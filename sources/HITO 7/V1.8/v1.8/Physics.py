import numpy as np

def Kepler(U,t):
    
    F0 = U[2] # dxdt
    F1 = U[3] # dydt
    F2 = -U[0]/(U[0]**2 + U[1]**2)**(3/2) # -x/(x^2+y^2)3/2
    F3 = -U[1]/(U[0]**2 + U[1]**2)**(3/2) # -y/(x^2+y^2)3/2

    return np.array([F0, F1, F2, F3])

def Linear_Oscilator(U, t):
    r = U[0]
    drdt = U[1]
    return np.array([drdt, -r])
import numpy as np
from Cauchy_Problem import Cauchy_problem
from Temporal_Schemes import Euler, Crank_Nicolson, RK4, Inverse_Euler, LeapFrog
from fNbody import F_NBody
from initialConditions import initial_cond
from Physics import Kepler

print('Enter number of bodies:')
Nb = int(input())#1024*72*2 # nº cuerpos tiene que ser 1024 (nº threads por bloque) * 72 (nº bloques o SM) * A (nº entero a elegir)
Nc = 3 # nº coord

t0 = 0 # timepo inicial
print('Enter time [s]:')
tf = int(input()) # tiempo final
N = 10*tf # number of time steps

t = np.linspace(t0, tf, N)
U0,M = initial_cond(Nb,Nc)

U = Cauchy_problem(RK4, F_NBody, t, U0, Nb, M)

U_s = np.reshape( U, (N, Nb, 2, Nc) ) 
r = np.reshape( U_s[:, :, 0, :], (N, Nb, Nc) )
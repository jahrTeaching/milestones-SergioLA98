import numpy as np

def Cauchy_problem(Temporal_Scheme, F, t, U0, Nb, M): 

     N = len(t) - 1
     N0 = len(U0)
     U = np.zeros((N+1, N0)) 

     U[0,:] = U0

     for i in range(N):
        #print(i)
        U[i+1,:] = Temporal_Scheme(U[i,:], t[i+1] - t[i], t[i], F, Nb, M) 
        #print(U[i+1,:])

     return U

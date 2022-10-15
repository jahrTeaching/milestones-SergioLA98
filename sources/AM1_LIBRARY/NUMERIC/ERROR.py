# -*- coding: utf-8 -*-
from numpy import linspace, zeros
from numpy.linalg import norm
from AM1_LIBRARY.ODES.EDO import CauchyProblem





def ErrorNum( Kepler, Scheme, t, U0):

    N = len(t)
    tf = t[N-1]
    E = zeros([N, len(U0)])
    t1 = t #N puntos, va de 0 a N-1 (tf), N puntos
    t2 = linspace(0,tf,2*N) #2a malla con doble de puntos
    # recordamos que el 0 cuenta y es [0,N) = [0, N-1]
    U1 = CauchyProblem(Kepler, t1, U0, Scheme)
    U2 = CauchyProblem(Kepler, t2, U0, Scheme)
    
    k = 10 #numero de puntos que pinto en la graf
    for i in range(0,N):
        E[i,:] = norm((U2[2*i,:] - U1[i,:]))
    
    return E


def ConvergenceRate( U, t, N):
    
    
    
    return
# -*- coding: utf-8 -*-
from numpy import linspace, zeros, log10, array, round_, reshape, float64
from numpy.linalg import norm 
from AM1_LIBRARY.ODES.EDO import CauchyProblem
from sklearn.linear_model import LinearRegression


def ConvergenceRate( Kepler, Scheme, t, U0, m):
    

    N = len(t)
    N2 = 2*len(t)
    tf = t[N-1]
    logE = zeros(m)
    logN = zeros(m)
    t1 = t #N puntos, va de 0 a N-1 (tf), N puntos
    t2 = linspace(0,tf,N2) #2a malla con doble de puntos
    # recordamos que el 0 cuenta y es [0,N) = [0, N-1]
    U1 = CauchyProblem(Kepler, t1, U0, Scheme)
    

    for i in range(m):
        
        
        U2 = CauchyProblem(Kepler, t2, U0, Scheme)
        logE[i] = log10(norm(U2[-1,:]-U1[-1,:]))
        logN[i] = log10(N)
        
        if i == (m-1):
            break
        
        t1 = t2
        U1 = U2
        N = N2
        N2 = 2*N2
        t2 = linspace(0,tf,N2)
       
    for j in range(m): 
             if abs(logE[j]) > 12 :  break #evito tramo final si aumenta
    j = min(j, m-1) 
    reg = LinearRegression().fit(logN[0:j+1].reshape((-1, 1)),logE[0:j+1]) 
    order = round_(abs(reg.coef_),1)
    mError = reg.coef_
    logN_lineal = logN[0:j+1]
    logE_lineal = reg.predict(logN[0:j+1].reshape((-1, 1)))
    logE_total = logE[:] - log10(1 - 1 / (2**order))

    #x = logN[0:j+1];  y = logE[0:j+1]
    #A = vstack( [ x, ones(len(x)) ] ).T #Numpy regression example
    #mError, cError = lstsq(A, y, rcond=None)[0]   #m pendiente y c coef
    #order = abs(mError)

    return logE, logN, logE_lineal, logN_lineal, order, mError, logE_total
    

    #logE = logE - log10( 1 - 1./2**order) #resto lo desplazado
    
    #return logE, logN, order, cError, mError



def ErrorNum( Kepler, Scheme, t, U0, order):

    N = len(t)
    tf = t[N-1]
    E = zeros([N, len(U0)])
    t1 = t #N puntos, va de 0 a N-1 (tf), N puntos
    t2 = linspace(0,tf,2*N) #2a malla con doble de puntos
    # recordamos que el 0 cuenta y es [0,N) = [0, N-1]
    U1 = CauchyProblem(Kepler, t1, U0, Scheme)
    U2 = CauchyProblem(Kepler, t2, U0, Scheme)
    
    
    for i in range(0,N):
        E[i,:] = (U2[2*i,:] - U1[i,:])/( 1 - 1./(2**order))
    
    return E

def StabilityRegion(Scheme):
    #PARA METODOS UNIPASO
    N = 100
    x, y = linspace(-5, 5, 100), linspace(-5, 5, 100)
    rho =  zeros( (N, N),  dtype=float64)

    for i in range(N): 
      for j in range(N):

          w = complex(x[i], y[j])
          r = Scheme( 1., 1., 0., lambda u, t: w*u ) #lambda es una vable q puede ser lo q sea luego
          # Scheme( U, dt, t, F )
          rho[i, j] = abs(r)
          
    return rho, x, y
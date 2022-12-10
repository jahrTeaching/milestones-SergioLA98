from numpy import zeros, reshape, array
from numpy.linalg import norm


def F_NBody(U, t, Nb, Nc):
    
    """
----------------------------------------------------------
N_BODY

This Module includes different temp schemes to integrate numerical problems

-Inputs:
    
    U: State Vector at integration time tn
    r,v: position and velocity
    dr/dt = v
    d = r(j) - r(i)
    dv/dt = sum(j)[Gm(j)d/abs(d)^3]
    F(U,t): Function dU/dt = F(U,t)
    
-Output: the return
    
    U State Vector at desired step tn + dt
----------------------------------------------------------
"""
    F =  zeros(len(U))
    Us  = reshape(U,(Nb, 2, Nc)) #automatic update by pointers
    dUs = reshape(F,(Nb, 2, Nc))

    r = reshape(Us[:,0,:],(Nb, Nc))
    v = reshape(Us[:,1,:],(Nb, Nc))
    dr = reshape(dUs[:,0,:], (Nb,Nc))
    dv = reshape(dUs[:,1,:], (Nb,Nc))
    
    for i in range(Nb):
        dr[i,:] = v[i,:]
        for j in range(Nb):
            if j != i:
                d = r[j,:] - r[i,:]
                dv[i,:] = dv[i,:] + d[:]/norm(d)**3 

    return F
    
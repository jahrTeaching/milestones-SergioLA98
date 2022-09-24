

#MILESTONE 2: INTEGRATE KEPLER WITH EULER, CRANCK-NICOLSON AND RK4 METHOD AND DISCUSS


from numpy import array, zeros, linspace
import matplotlib.pyplot as plt
    







def prueba():

    N = 200
    U0 = array( [ 1, 0, 0, 1 ])
    t = linspace(0, 20, num = N)

    U = CauchyProblem( Kepler, t, U0, Euler)

    plt.plot(U[:,0] ,U[:,1])
    plt.show()


def Kepler(U, t): 

    x = U[0]; y = U[1]; dxdt = U[2]; dydt = U[3]
    d = ( x**2  +y**2 )**1.5

    return  array( [ dxdt, dydt, -x/d, -y/d ] ) 


def Euler(U, dt, t, F): 

    return U + dt * F(U, t)


def CauchyProblem( F, t, U0, TemporalScheme ):

    N, Nv=  len(t)-1, len(U0)
    U = array( zeros([N+1, Nv] ) )

    U[0,:] = U0

    for n in range(N):

       U[n+1,:] = TemporalScheme( U[n, :], t[n+1] - t[n], t[n],  F ) 

    return U
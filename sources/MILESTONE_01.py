
#MILESTONE 1. CREATE EULER, CRANK NICOLSON, RK4 SCRIPT TO INTEGRATE KEPLER ORBIT

import matplotlib.pyplot as plt
from matplotlib import rc #LaTex Font
from numpy import array, zeros, linspace, sqrt
from scipy.optimize import fsolve


rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True); plt.rc('font', family='serif')


def MIL01():


    #TIME and TIME STEP
    tf = 20 
    dt = [ 0.001]#, 0.02, 0.05, 0.001]
    
    #INITIAL CONDITIONS U0 MATRIX
    U0 = array( [1, 0, 0, 1])

    #SOLVING MATRIX

    U_Euler = {} #Dictionary to allocate different size matrix for different time steps and same total time
    U_CrankNic = {}
    U_RK4 = {}

    for in_dt in dt: #initialize different inputs of the dictionaries

        N = int( tf/in_dt )+ 1 #Time Nodes
        t = linspace( 0, tf, N) #vector tiempos
        U_Euler [in_dt] = zeros( ( len(U0), N ))
        U_CrankNic [in_dt] = zeros( ( len(U0), N ))
        U_RK4 [in_dt] = zeros( ( len(U0), N ))

        U_Euler [ in_dt] [ :, 0] = U0
        U_CrankNic [ in_dt] [ :, 0] = U0
        U_RK4 [ in_dt] [ :, 0] = U0

    for i in range ( 1, N): #Initial Node is in U0

        U_Euler [ in_dt] [ :, i] = Euler( U_Euler [ in_dt] [ :, i-1], t, in_dt, F_Kepler)
        U_CrankNic [ in_dt] [ :, i] = CrankNicolson( U_CrankNic [ in_dt] [ :, i-1], t, in_dt, F_Kepler )
        U_RK4 [ in_dt] [ :, i]     = RK4( U_RK4 [ in_dt] [ :, i-1], t, in_dt, F_Kepler )


        ## PLOTTING DIFFERENT SCHEMES

     #EULER
    x = linspace(-1,1,1000)
    y = sqrt(1-x**2)

    colours = ['red', 'orange', 'blue', 'purple']

    fig, ax = plt.subplots(1,1, figsize=(11,11), constrained_layout='true')
    ax.set_xlim(-1.5,1.5)
    ax.set_ylim(-1.5,1.5)
    ax.set_title('Kepler solution for Euler Scheme', fontsize=30)
    ax.grid()
    ax.set_xlabel(r'$x/r$',fontsize=15)
    ax.set_ylabel(r'$y/r$',fontsize=15)

    for i in range(len(dt)):
    
        ax.plot( U_Euler[dt[i]][0,:], U_Euler[dt[i]][1,:], c=colours[i], label=r'$\Delta t$ = '+str(dt[i]))

    ax.plot(x,y,'k'); ax.plot(x,-y,'k')
    ax.plot(0,0,'k-o', markersize=12)
    ax.legend(loc=0, fancybox=False, edgecolor="black", ncol = 1, fontsize=16)

    plt.show()

    #CRANK-NICOLSON
    fig, ax = plt.subplots(1,1, figsize=(11,11), constrained_layout='true')
    ax.set_xlim(-1.1,1.1)
    ax.set_ylim(-1.1,1.1)
    ax.set_title('Kepler solution for C-N Scheme', fontsize=30)
    ax.grid()
    ax.set_xlabel(r'$x/r$',fontsize=15)
    ax.set_ylabel(r'$y/r$',fontsize=15)

    for i in range(len(dt)):
    
        ax.plot( U_CrankNic[dt[i]][0,:], U_CrankNic[dt[i]][1,:], c=colours[i], label=r'$\Delta t$ = '+str(dt[i]))

    ax.plot(x,y,'k'); ax.plot(x,-y,'k')
    ax.plot(0,0,'k-o', markersize=12)
    ax.legend(loc=0, fancybox=False, edgecolor="black", ncol = 1, fontsize=16)

    plt.show()

    #RK4
    fig, ax = plt.subplots(1,1, figsize=(11,11), constrained_layout='true')
    ax.set_xlim(-1.1,1.1)
    ax.set_ylim(-1.1,1.1)
    ax.set_title('Kepler solution for RK4 Scheme', fontsize=30)
    ax.grid()
    ax.set_xlabel(r'$x/r$',fontsize=15)
    ax.set_ylabel(r'$y/r$',fontsize=15)

    for i in range(len(dt)):
    
        ax.plot( U_RK4[dt[i]][0,:], U_RK4[dt[i]][1,:], c=colours[i], label=r'$\Delta t$ = '+str(dt[i]))

    ax.plot(x,y,'k'); ax.plot(x,-y,'k')
    ax.plot(0,0,'k-o', markersize=12)
    ax.legend(loc=0, fancybox=False, edgecolor="black", ncol = 1, fontsize=16)

    plt.show()





def F_Kepler( U, t): #U( x, y, x', y' ) dU/dt = F(U,t) = ( x', y', x'', y'' )

    F1 = U[2]     
    F2 = U[3]    
    F3 = -U[0]/(U[0]**2 + U[1]**2)**(3/2)
    F4 = -U[1]/(U[0]**2 + U[1]**2)**(3/2)

    return array( [ F1, F2, F3, F4] )


def Euler( U, t, dt, F_Kepler):
    return U + dt * F_Kepler( U, t)


def CrankNicolson( U, t, dt, F_Kepler): #U(n+1) = U(n) + 0.5 * dt * ( F(Un+1,tn+1) + F(Un,tn))
   
   def CrankNicSol(X):

        return X - U - dt/2 * ( F_Kepler(U,t) + F_Kepler(X,t+dt) )

   return fsolve( CrankNicSol, U ) 


def RK4( U, t, dt, F_Kepler): #U(n+1) = U(n) + (k1 + 2*k2 + 2*k3 + k4)/6

    k1 = F_Kepler( U, t )
    k2 = F_Kepler( U + dt * k1/2, t + dt/2 )    
    k3 = F_Kepler( U + dt * k2/2, t + dt/2 )    
    k4 = F_Kepler( U + dt *k3, t + dt   )

    return U + dt * (k1 + 2 * k2 + 2 * k3 + k4)/6
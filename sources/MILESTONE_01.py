
#MILESTONE 1. CREATE EULER, CRANK NICOLSON, RK4 SCRIPT TO INTEGRATE KEPLER ORBIT

import matplotlib.pyplot as plt
from numpy import array, zeros, linspace, sqrt

def Euler():


    #TIME and TIME STEP
    tf = 20 
    dt = [ 0.1, 0.2, 0.02]

    #INITIAL CONDITIONS U0 MATRIX
    U0 = array( [1, 0, 0, 1])

    #SOLVING MATRIX

    U_Euler = {} #Dictionary to allocate different size matrix for different time steps and same total time
    U_Euler = array( zeros( ( 4, N)))
    U_Euler[ :, 0] = U0

    for i in range ( 1, N):
        F = array([U_Euler[ 2, i-1], U_Euler[ 3, i-1],
                   -U_Euler[0,i-1]/(U_Euler[0,i-1]**2+U_Euler[1,i-1]**2)**(3/2), 
                   -U_Euler[1,i-1]/(U_Euler[0,i-1]**2+U_Euler[1,i-1]**2)**(3/2)])
        U_Euler[ :, i] = U_Euler[ :, i-1] + dt * F


    plt.figure(1)
    plt.plot(U_Euler[0,:],U_Euler[1,:])
    plt.show()

    #PLOTTING EULER SCHEME
    x = linspace(-1,1,1000)
    y = sqrt(1-x**2)

    colours = ['purple', 'orange', 'blue', 'red']

    fig, ax = plt.subplots(1,1, figsize=(11,11), constrained_layout='true')
    ax.set_xlim(-1.85,1.85)
    ax.set_ylim(-1.85,1.85)
    ax.set_title('Kepler solution for Euler Scheme', fontsize=30)
    ax.grid()
    ax.set_xlabel(r'$x$',fontsize=15)
    ax.set_ylabel(r'$y$',fontsize=15)

    #for i in range(len(dt)):
    
    ax.plot( U_Euler[0,:], U_Euler[1,:], c=colours[1], label=r'$\Delta t$ = '+str(dt))

    ax.plot(x,y,'k'); ax.plot(x,-y,'k')
    ax.plot(0,0,'k-o', markersize=12)
    ax.legend(loc=0, fancybox=False, edgecolor="black", ncol = 1, fontsize=16)

    plt.show()



def Crank_Nicolson():

    #NUMBER OF TIME STEP
    N = 200 #Different time steps

    #TIME
    tf = 20 
    dt = tf / N

    #INITIAL CONDITIONS U0 MATRIX
    U0 = array( [1, 0, 0, 1])

    #SOLVING MATRIX
    U = array( zeros( ( 4, N)))
    U[ :, 0] = U0

    for i in range ( 1, N):
        F = array([U[ 2, i-1], U[ 3, i-1],
                   -U[0,i-1]/(U[0,i-1]**2+U[1,i-1]**2)**(3/2), 
                   -U[1,i-1]/(U[0,i-1]**2+U[1,i-1]**2)**(3/2)])
        U[ :, i] = U[ :, i-1] + dt * F


    plt.figure(1)
    plt.plot(U[0,:],U[1,:])
    plt.show()

    #PLOTTING EULER SCHEME
    x = linspace(-1,1,1000)
    y = sqrt(1-x**2)

    colours = ['purple', 'orange', 'blue', 'red']

    fig, ax = plt.subplots(1,1, figsize=(11,11), constrained_layout='true')
    ax.set_xlim(-1.85,1.85)
    ax.set_ylim(-1.85,1.85)
    ax.set_title('Kepler solution for Euler Scheme', fontsize=30)
    ax.grid()
    ax.set_xlabel(r'$x$',fontsize=15)
    ax.set_ylabel(r'$y$',fontsize=15)

    #for i in range(len(dt)):
    
    ax.plot( U[0,:], U[1,:], c=colours[1], label=r'$\Delta t$ = '+str(dt))

    ax.plot(x,y,'k'); ax.plot(x,-y,'k')
    ax.plot(0,0,'k-o', markersize=12)
    ax.legend(loc=0, fancybox=False, edgecolor="black", ncol = 1, fontsize=16)

    plt.show()

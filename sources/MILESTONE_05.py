from AM1_LIBRARY.ODES.TEMP_SCHEMES import RK4_Scheme, LeapFrog
from AM1_LIBRARY.ODES.EDO import CauchyProblem
from AM1_LIBRARY.PHYSICS.NBODY import F_NBody

from numpy import linspace, size, zeros, reshape
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt



def Mil5( tf, N, Nb, Nc):
    
    #Nb numero de cuerpos, de 2 a 4
    #Nc numero de coordenadas
        
    def Init_NBody(Nb, Nc): #esta en el scope al estar dentro de Mil5
        
        U0 = zeros(2*Nb*Nc)
        U_0 = reshape(U0, (Nb, 2, Nc))
        r0 = reshape(U_0[:,0,:], (Nb, Nc)) #actualiza solo al ser puntero 
        v0 = reshape(U_0[:,1,:], (Nb, Nc)) #actualiza solo al ser puntero 
        
        if Nb == 2:
            #Cuerpo 1
            r0[0,:] = [2., 2., 0.]
            v0[0,:] = [-0.4, 0., 0.]
            #Cuerpo 2
            r0[1,:] = [-2., 2., 0.]
            v0[1,:] = [0., -0.4, 0.]
            
        elif Nb == 3:
            #Cuerpo 1
            r0[0,:] = [2., 2., 0.]
            v0[0,:] = [-0.4, 0., 0.]
            #Cuerpo 2
            r0[1,:] = [-2., 2., 0.]
            v0[1,:] = [0., -0.4, 0.]
            #Cuerpo 3 
            r0[2, :] = [ -2., -2., 0. ] 
            v0[2, :] = [ 0.4, 0., 0. ] 
            
        else:
            #Cuerpo 1
            r0[0,:] = [2., 2., 0.]
            v0[0,:] = [-0.5, 0., 0.]
            #Cuerpo 2
            r0[0,:] = [-2., 2., 0.]
            v0[0,:] = [0., -0.5, 0.]
            #Cuerpo 3 
            r0[2, :] = [ -2., -2., 0. ] 
            v0[2, :] = [ 0.5, 0., 0. ]  
            #Cuerpo 4 
            r0[3, :] = [ 2., -2., 0. ] 
            v0[3, :] = [ 0., 0.5, 0. ]
            
        return U0
    
    def F(U,t): #Wrapper para la F de Cauchy al tener más argumentos
        return F_NBody(U, t, Nb, Nc)
    
    U0 = Init_NBody(Nb, Nc)
    #Scheme = RK4_Scheme
    Scheme = LeapFrog
    for i in range(size(N)):
        dt = tf/(N[i]-1.)
        t = linspace(0, tf, N[i]) #N puntos el linspace
        U = CauchyProblem(F, t, U0, Scheme)
        #print(size(U)) 
        U_rec = reshape(U,(N[i], Nb, 2, Nc))
        r_rec = reshape(U_rec[:,:,0,:],(N[i],Nb,Nc))
        color = ['b','r','g','y']
        
        for j in range(Nb):
            plt.plot( r_rec[:,j,0],r_rec[:,j,1], color=color[j], label = 'Orbit' + str(j+1)) #xy trayectoria, es la r
            
        plt.axis('tight')
        plt.xlim([-5,5])
        plt.ylim([-5,5])
        plt.title(f'N = {Nb} body problem. XY Projection')
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.legend(loc="best")
        plt.show()
        
        fig = plt.figure(2)
        ax = fig.add_subplot(projection='3d')
        ax.plot(r_rec[:, 0, 0], r_rec[:, 0, 1], r_rec[:, 0, 2], color[0], label = 'Orbit 1')
        ax.plot(r_rec[:, 1, 0], r_rec[:, 1, 1], r_rec[:, 1, 2], color[1], label = 'Orbit 2')
        ax.plot(r_rec[:, 2, 0], r_rec[:, 2, 1], r_rec[:, 2, 2], color[2], label = 'Orbit 3')
        ax.plot(r_rec[:, 3, 0], r_rec[:, 3, 1], r_rec[:, 3, 2], color[3], label = 'Orbit 4')
        plt.title(f'N = {Nb} body problem. 3D')
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.axis('tight')
        plt.legend(loc = "best")
        ax.set_xlim3d(-5, 5) # viewrange for z-axis should be [-4,4] 
        ax.set_ylim3d(-5,5) # viewrange for y-axis should be [-2,2] 
        ax.set_zlim3d(-5, 5) # viewrange for x-axis should be [-2,2] 
        #plt.grid(True)
        plt.show()  
        
    
from AM1_LIBRARY.PHYSICS.OSCILATOR import LinearOscilator
from AM1_LIBRARY.ODES.TEMP_SCHEMES import Euler_Scheme, InvEuler_Scheme, LeapFrog, CrankNicolson_Scheme, RK4_Scheme
from AM1_LIBRARY.ODES.EDO import CauchyProblem
from AM1_LIBRARY.NUMERIC.ERROR import StabilityRegion

from numpy import zeros, array, linspace, size, transpose
import matplotlib.pyplot as plt

def Mil4( tf, N, U0):
    
    schemes = [ Euler_Scheme, RK4_Scheme, CrankNicolson_Scheme, InvEuler_Scheme, LeapFrog ]
    colors = ['b--','g--','r--']
    
    
    for i in range(size(schemes)):
        
        for j in range(size(N)):
            dt = tf/(N[j]-1)
            t = linspace(0, tf, N[j]) #da N puntos el linspace
            U = CauchyProblem( LinearOscilator, t, U0, schemes[i])
            
            plt.title(f'Oscilator with {schemes[i].__name__}')
            plt.xlabel("Tiempo (s)")
            plt.ylabel("X",rotation = 0)
            plt.grid()
            plt.plot(t,U[:,0], colors[j], label = 'dt =' + str(dt) + ' s')
            plt.legend(loc ='lower left') 
            
        plt.savefig('MILESTONE 4/' + 'Oscilator with'+ schemes[i].__name__+'.png')
        #plt.show()
        plt.close()
        
        
        rho,x,y = StabilityRegion(schemes[i])
        
        plt.title(f'Stability Region of {schemes[i].__name__}')
        plt.xlabel("Re")
        plt.ylabel("Im",rotation = 0)
        #plt.axis('equal')
        plt.grid()
    
        if schemes[i] == LeapFrog:
            Im = linspace(-5,5,100)
            Re = zeros(100)
            plt.plot(Re, Im, color = 'y')
            #plt.show()
        else:
            plt.contour( x, y, transpose(rho), linspace(0, 1, 20) )
            plt.contourf( x, y, transpose(rho), levels = [0,1], colors = 'whitesmoke')
        
        color2 = ['b','g','r']
        for j in range(size(N)):
            dt = tf/(N[j]-1)
            plt.plot([0,0], [dt,-dt], 'o', color = color2[j], label = 'dt =' + str(dt) + ' s')
            
        plt.legend()
        plt.xlim([-5,5])
        plt.ylim([-5,5])
        plt.savefig('MILESTONE 4/' + 'SR of'+ schemes[i].__name__+'.png')
        #plt.show()
        plt.close()
        
   
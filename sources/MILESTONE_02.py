from AM1_LIBRARY.ODES.EDO import CauchyProblem 
from AM1_LIBRARY.ODES.TEMP_SCHEMES import Euler_Scheme, CrankNicolson_Scheme
from AM1_LIBRARY.ODES.TEMP_SCHEMES import RK4_Scheme, InvEuler_Scheme
from AM1_LIBRARY.PHYSICS.ORBITS import Kepler 
import matplotlib.pyplot as plt
#from numpy import where

 
from numpy import linspace


def Simulation(tf, N, U0): 
   
    t = linspace(0, tf, N)
    schemes = [ Euler_Scheme, RK4_Scheme, CrankNicolson_Scheme, InvEuler_Scheme ]
    LegendScheme = ["Euler Scheme", "RK4 Scheme", "Cranck Nicolson Scheme", "InvEuler_Scheme"]
    
    for method in range(len(schemes)-0):
       U =  CauchyProblem( Kepler, t, U0, schemes[method]) 
       print(method)
       fig, ax = plt.subplots(1,1, figsize=(11,11), constrained_layout='true')
       ax.set_xlim(-1.85,1.85)
       ax.set_ylim(-1.85,1.85)
       ax.set_title("Kepler for "+ LegendScheme[method], fontsize=30)
       ax.grid()
       ax.set_xlabel(r'$x/r$',fontsize=15)
       ax.set_ylabel(r'$y/r$',fontsize=15)
       plt.plot(U[:,0] , U[:,1])
       plt.show()
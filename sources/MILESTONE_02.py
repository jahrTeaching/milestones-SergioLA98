from AM1_LIBRARY.ODES.EDO import CauchyProblem 
from AM1_LIBRARY.ODES.TEMP_SCHEMES import Euler_Scheme, CrankNicolson_Scheme
from AM1_LIBRARY.ODES.TEMP_SCHEMES import RK4_Scheme, InvEuler_Scheme
from AM1_LIBRARY.PHYSICS.ORBITS import Kepler 
import matplotlib.pyplot as plt
#from numpy import where

 
from numpy import linspace


def SimulationMIL2(tf, N, U0): 

    t = linspace(0, tf, N) #va de 0 a N-1 -> N puntos
    schemes = [ Euler_Scheme, RK4_Scheme, CrankNicolson_Scheme, InvEuler_Scheme ]
    LegendScheme = ["Euler Scheme", "RK4 Scheme", "Cranck Nicolson Scheme", "InvEuler Scheme"]
    delta = round(tf/(N-1),4)
#    for indic in range(len(schemes)):
    for indic in range(len(schemes)): #(a,b] en in range, enumerate
       method = schemes[indic] 
       U =  CauchyProblem( Kepler, t, U0, method) 
       fig, ax = plt.subplots(1,1, figsize=(11,11), constrained_layout='true')
       ax.set_xlim(-1.85,1.85)
       ax.set_ylim(-1.85,1.85)
       ax.set_title("Kepler for "+ LegendScheme[indic], fontsize=30)
       ax.grid()
       ax.set_xlabel(r'$x/r$',fontsize=15)
       ax.set_ylabel(r'$y/r$',fontsize=15)
       plt.plot(U[:,0] , U[:,1],label=r'$\Delta t$ = ' + str(delta))
       ax.legend(loc=0, fancybox=False, edgecolor="black", ncol = 1, fontsize=16)
       plt.show( )
from AM1_LIBRARY.PHYSICS.OSCILATOR import LinearOscilator
from AM1_LIBRARY.ODES.TEMP_SCHEMES import Euler_Scheme, InvEuler_Scheme, LeapFrog, CrankNicolson_Scheme, RK4_Scheme
from AM1_LIBRARY.ODES.EDO import CauchyProblem

from numpy import zeros, array, linspace, size
import matplotlib.pyplot as plt

def Mil4( tf, N, U0):
    
    schemes = [ Euler_Scheme, RK4_Scheme, CrankNicolson_Scheme, InvEuler_Scheme, LeapFrog ]
    colors = ['b--','g--','r--']
    
    
    for i in range(size(schemes)):
        
        for j in range(size(N)):
            
            t = linspace(0, tf, N[j]) #da N puntos el linspace

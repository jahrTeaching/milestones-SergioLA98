# -*- coding: utf-8 -*-

from numpy import array, linspace, log10
import matplotlib.pyplot as plt

from AM1_LIBRARY.NUMERIC.ERROR import ErrorNum, ConvergenceRate
from AM1_LIBRARY.ODES.EDO import CauchyProblem
from AM1_LIBRARY.ODES.TEMP_SCHEMES import Euler_Scheme, CrankNicolson_Scheme
from AM1_LIBRARY.ODES.TEMP_SCHEMES import InvEuler_Scheme, RK4_Scheme
from AM1_LIBRARY.PHYSICS.ORBITS import Kepler



def ErrorNumerico3( tf, N, U0):
    
    schemes = [ Euler_Scheme, RK4_Scheme, CrankNicolson_Scheme, InvEuler_Scheme ]
    t = linspace(0, tf, N) #[0, N) hasta N-1 pero N puntos
    deltat = round(tf/(N-1),4)
    for Scheme in schemes:
        Error = ErrorNum(Kepler, Scheme, t, U0)
        plt.plot(t,Error[:,0])
        plt.show()
        



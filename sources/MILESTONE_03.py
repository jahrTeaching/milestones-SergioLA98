# -*- coding: utf-8 -*-

from numpy import array, linspace, log10, round_
import matplotlib.pyplot as plt

from AM1_LIBRARY.NUMERIC.ERROR import ErrorNum, ConvergenceRate
from AM1_LIBRARY.ODES.EDO import CauchyProblem
from AM1_LIBRARY.ODES.TEMP_SCHEMES import Euler_Scheme, CrankNicolson_Scheme
from AM1_LIBRARY.ODES.TEMP_SCHEMES import InvEuler_Scheme, RK4_Scheme
from AM1_LIBRARY.PHYSICS.ORBITS import Kepler



def ErrorNumerico3( tf, N, U0):
    
    
    schemes = [ Euler_Scheme, RK4_Scheme, CrankNicolson_Scheme, InvEuler_Scheme ]
    GrafSchemes = ['Euler Scheme', 'RK4 Scheme', 'Cranck-Nicolson Scheme', 'Inverse Euler Scheme']
    t = linspace(0, tf, N) #[0, N) hasta N-1 pero N puntos
    deltat = round(tf/(N-1),4)
    i = 0
    for Scheme in schemes:
        
        [logE, logN, logE_lineal, logN_lineal, order, mError, logE_total] = ConvergenceRate(Kepler, Scheme, t, U0, 8)
        #[logE, logN, order, cError, mError] = ConvergenceRate(Kepler, Scheme, t, U0, 8)
        #print('Order =', order, 'm=',mError, 'c=',cError)
        print('Order =', order, 'm=',mError)

        plt.plot(logN, logE, 'o', label = 'log(|U2-U1|)', markersize = 6)
        plt.plot(logN, logE_total, 'k', label='log(|E|)')
        #plt.plot(logN, mError*logN + cError, 'k', label='Fitted line')
        plt.legend(loc='upper right')
        plt.grid()
        plt.xlabel('log(N)')
        plt.ylabel('log(Error)')
        plt.title('Error of order' + ' ' + f'{order}' + ' '+ GrafSchemes[i])
        plt.savefig('MILESTONE 3/'+ 'Both error of ' + GrafSchemes[i]+ ' ' + str(deltat)+'.png')
        plt.show()



        plt.plot(logN, logE, 'o', label = 'Error', markersize = 6)
        plt.plot(logN_lineal, logE_lineal, 'k', label='Fitted line')
        #plt.plot(logN, mError*logN + cError, 'k', label='Fitted line')
        plt.legend(loc='upper right')
        plt.grid()
        plt.xlabel('log(N)')
        plt.ylabel('log(Error)')
        plt.title('Convergence rate of order' + ' ' + f'{order}' + ' '+ GrafSchemes[i])
        plt.savefig('MILESTONE 3/'+ 'Convergence rate of ' + GrafSchemes[i]+ ' ' + str(deltat)+'.png')
        plt.show()
        order = round_(abs(mError))
        
        Error = ErrorNum(Kepler, Scheme, t, U0, order = int(round_(abs(mError))))
        ErrorNorm = (Error[:,0]**2+Error[:,1]**2)**(1/2)
        plt.plot(t,ErrorNorm, label = GrafSchemes[i] + ' Error')
        plt.legend(loc = 'upper left')
        plt.xlabel('t (s)')
        plt.ylabel('Error')
        plt.grid()
        plt.title('Error norm of ' + GrafSchemes[i] + ', deltat = ' + str(deltat))
        plt.savefig('MILESTONE 3/'+ 'Error norm of ' + GrafSchemes[i]+ ' ' + str(deltat)+'.png')
        plt.show()
        i = i + 1


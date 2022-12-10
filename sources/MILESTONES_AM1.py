
from numpy import array 

#import MILESTONE_01
#import MILESTONE_02
#import MILESTONE_03
#import MILESTONE_04
#import MILESTONE_05
import MILESTONE_06

#MILESTONE_01.prueba()

#MILESTONE_02.SimulationMIL2(  tf = 20, N = 201, U0 = array( [ 1., 0., 0., 1. ] ) )

#MILESTONE_03.ErrorNumerico3( tf = 5, N = 101, U0 = array( [ 1., 0., 0., 1. ] ) )

#MILESTONE_04.Mil4( tf = 20, N = [ 201, 2001, 20001], U0 = array([1,0]))

#MILESTONE_05.Mil5( tf = 30, N = [20001], Nb = 4, Nc = 3)

MILESTONE_06.Mil6(mu=1.2151e-2, tf=350, N=20001, epsorder = 1e-3)
#mu = 3.0039e-7 #Tierra-Sol
#mu = 1.2151e-2 #Tierra-Luna
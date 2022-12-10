from numpy import zeros, linspace, array
import matplotlib.pyplot as plt

from AM1_LIBRARY.PHYSICS.ORBITS import Kepler
from AM1_LIBRARY.ODES.TEMP_SCHEMES import Embedded_RK
from AM1_LIBRARY.ODES.EDO import CauchyProblem



U0 = array( [ 1., 0., 0., 1. ] )
t = linspace(0, 20, 2001) #da N puntos el linspace
U = CauchyProblem( Kepler, t, U0, Embedded_RK)

plt.plot(U[:,0] , U[:,1])
plt.show()

'''

U0 = zeros(6,5) #6coordinates and 5 points
U0[:,0] = [ 0.8, 0.6, 0., 0., 0., 0.  ]   # Lagrange points
U0[:,1] = [ 0.8, -0.6, 0., 0., 0., 0.  ]
U0[:,2] = [ -0.1, 0.0, 0., 0., 0., 0.  ]
U0[:,3] = [ 0.1, 0.0, 0., 0., 0., 0.  ]
U0[:,4] = [ 1.1, 0.0, 0., 0., 0., 0.  ]

'''

from numpy import zeros, linspace, array
import matplotlib.pyplot as plt

from AM1_LIBRARY.PHYSICS.ORBITS import CR3BP, LG_Points
from AM1_LIBRARY.ODES.TEMP_SCHEMES import Embedded_RK
from AM1_LIBRARY.ODES.EDO import CauchyProblem


#mu = 3.0039e-7 #Tierra-Sol
mu = 1.2151e-2 #Tierra-Luna

U0_LG = zeros([5,4]) #5 Puntos Lagrange, 2 coord. De normal era vector 1,4 lo que mando a Cauchy
U0_LG[0,:] = [0.8, 0.6, 0, 0]
U0_LG[1,:] = [0.8, -0.6, 0, 0]
U0_LG[2,:] = [-0.1, 0, 0, 0]
U0_LG[3,:] = [0.1, 0, 0, 0]
U0_LG[4,:] = [1.01, 0, 0, 0]
t = linspace(0, 20, 2001) #da N puntos el linspace

LP = LG_Points(U0_LG, 5, mu)
print(LP)
def F(U,t):
    return CR3BP(U,t,mu)  #wrapped ya que Schemes usa F(U,t) y aqui doy mu tb

#U = CauchyProblem( F, t, U0, Embedded_RK)

#plt.plot(U[:,0] , U[:,1])
#plt.show()

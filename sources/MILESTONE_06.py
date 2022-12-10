from numpy import zeros, linspace, array, around, size
from random import random
import matplotlib.pyplot as plt

from AM1_LIBRARY.PHYSICS.ORBITS import CR3BP, LG_Points, LP_Stab
from AM1_LIBRARY.ODES.TEMP_SCHEMES import Embedded_RK, Euler_Scheme
from AM1_LIBRARY.ODES.EDO import CauchyProblem


#mu = 3.0039e-7 #Tierra-Sol
mu = 1.2151e-2 #Tierra-Luna
tf = 1000
U0_LG = zeros([5,4]) #5 Puntos Lagrange, 2 coord. De normal era vector 1,4 lo que mando a Cauchy
U0_LG[0,:] = [0.8, 0.6, 0, 0]
U0_LG[1,:] = [0.8, -0.6, 0, 0]
U0_LG[2,:] = [-0.1, 0, 0, 0]
U0_LG[3,:] = [0.1, 0, 0, 0]
U0_LG[4,:] = [1.01, 0, 0, 0]
t = linspace(0, tf, 2001) #da N puntos el linspace

LP = LG_Points(U0_LG, 5, mu)
print(LP)
plt.plot(LP[:,0],LP[:,1],'o')
plt.show()
#LP_Stab_VAL = zeros([4])
#Autoval_LP = LP_Stab(U0_LG, mu) #estabilidad
#print(Autoval_LP)

def F(U,t):
    return CR3BP(U,t,mu)  #wrapped ya que Schemes usa F(U,t) y aqui doy mu tb

#U = CauchyProblem( F, t, U0, Embedded_RK)

#plt.plot(U[:,0] , U[:,1])
#plt.show()
U0_LP_Sel = zeros(4)
U0_LP_SelStab = zeros(4)
eps = 1e-3*random()
LP_List = array([1,2,3,4,5])
for k in range(5):

    sel = k + 1

    if sel == 5:
        label = 'L2'
    elif sel == 4:
        label = 'L1'
    elif sel == 3:
        label = 'L3'
    elif sel == 2:
        label = 'L5'
    elif sel == 1:
        label = 'L4'
    
    U0_LP_Sel[0:2] = LP[sel-1,:] + eps
    U0_LP_Sel[2:4] = eps

    U0_LP_SelStab[0:2] = LP[sel-1,:]
    U0_LP_SelStab[2:4] = 0

    Autoval_LP = LP_Stab(U0_LP_SelStab, mu) #estabilidad
    print(around(Autoval_LP.real,8))

    methods = [Euler_Scheme]#[Embedded_RK]

    for j in range (size(methods)):

        U_LP = CauchyProblem(F, t, U0_LP_Sel, methods[j])

        fig, (ax1, ax2) = plt.subplots(1, 2)
        ax1.plot(U_LP[:,0], U_LP[:,1],'-',color = "r")
        ax1.plot(-mu, 0, 'o', color = "g")
        ax1.plot(1-mu, 0, 'o', color = "b")
        for i in range(5):
            ax1.plot(LP[i,0], LP[i,1] , 'o', color = "k")

        ax2.plot(U_LP[:,0], U_LP[:,1],'-',color = "r")
        ax2.plot(LP[sel - 1,0], LP[sel - 1,1] , 'o', color = "k")

        ax1.set_xlim(-2,2)
        ax1.set_ylim(-2,2)
        ax1.set_title("Vista del sistema orbital")
        ax2.set_title("Vista del punto de Lagrange")
        ax2.set_xlim(LP[sel - 1,0]-0.02,LP[sel - 1,0]+0.02)
        ax2.set_ylim(LP[sel - 1,1]-0.02,LP[sel - 1,1]+0.02)
        fig.suptitle(f"Tierra-Luna - CR3BP ({methods[j].__name__}) - Órbita alrededor de {label} con t = {tf}s" )
        for ax in fig.get_axes():
            ax.set(xlabel='x', ylabel='y')
            ax.grid()
            
        #manager = plt.get_current_fig_manager()
        #manager.full_screen_toggle()
        figure = plt.gcf()                  
        #figure.set_size_inches(16, 8)       
        #plt.savefig('Plots/Hito 6/ CR3BP ' + label +' '+ methods[j].__name__ +'.png', bbox_inches = 'tight')
        #plt.close('all')
        plt.show()

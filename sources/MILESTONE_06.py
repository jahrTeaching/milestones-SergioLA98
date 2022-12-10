from numpy import zeros, linspace, array, around, size
from random import random
import matplotlib.pyplot as plt

from AM1_LIBRARY.PHYSICS.ORBITS import CR3BP, LG_Points, LP_Stab
from AM1_LIBRARY.ODES.TEMP_SCHEMES import Embedded_RK, Euler_Scheme, InvEuler_Scheme, CrankNicolson_Scheme
from AM1_LIBRARY.ODES.EDO import CauchyProblem


#mu = 3.0039e-7 #Tierra-Sol
#mu = 1.2151e-2 #Tierra-Luna
#tf = 500
#N = 20001

def Mil6(mu, tf, N, epsorder):
    eps = epsorder*random()
    dt = tf/(N-1)
    methods = [Embedded_RK]
    #methods = [CrankNicolson_Scheme]
    
    
    U0_LG = zeros([5,4]) #5 Puntos Lagrange, 2 coord. De normal era vector 1,4 lo que mando a Cauchy
    U0_LG[0,:] = [0.8, 0.6, 0, 0]
    U0_LG[1,:] = [0.8, -0.6, 0, 0]
    U0_LG[2,:] = [-0.1, 0, 0, 0]
    U0_LG[3,:] = [0.1, 0, 0, 0]
    U0_LG[4,:] = [1.01, 0, 0, 0]
    t = linspace(0, tf, N) #da N puntos el linspace

    LPAUX = LG_Points(U0_LG, 5, mu)
    LP = zeros([5,2])
    LP[0,:] = LPAUX[3,:] #Reordeno
    LP[1,:] = LPAUX[4,:] 
    LP[2,:] = LPAUX[2,:] 
    LP[3,:] = LPAUX[0,:] 
    LP[4,:] = LPAUX[1,:] 
    labelPTot = ['L1','L2','L3','L4','L5'] #ordenados
    ShapeLP = ["<",">","d","^","v"]
    ColorLP = ["yellow","cyan","violet","sienna","lightcoral"]
    print(LP)
    for i in range(5):
        plt.plot(LP[i,0],LP[i,1],ShapeLP[i],color = ColorLP[i],label=labelPTot[i])
    plt.plot(-mu, 0, 'o', color = "g", label = 'Tierra')
    plt.plot(1-mu, 0, 'o', color = "b", label = 'Luna')
    plt.grid()
    plt.title("Puntos de Lagrange del CR3BP Tierra-Luna")
    plt.legend(loc = 'upper left',bbox_to_anchor=(1., 0.95))
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

    #LP_List = array([1,2,3,4,5])
    for k in range(5):

        sel = k + 1 #ya los he ordenado

        if sel == 1:
            labelP = 'L1'
        elif sel == 2:
            labelP = 'L2'
        elif sel == 3:
            labelP = 'L3'
        elif sel == 4:
            labelP = 'L4'
        elif sel == 5:
            labelP = 'L5'
        
        U0_LP_Sel[0:2] = LP[sel-1,:] + eps
        U0_LP_Sel[2:4] = eps

        U0_LP_SelStab[0:2] = LP[sel-1,:]
        U0_LP_SelStab[2:4] = 0

        Autoval_LP = LP_Stab(U0_LP_SelStab, mu) #estabilidad
        print(around(Autoval_LP.real,8))

        #methods = [Euler_Scheme]
        
        

        for j in range (size(methods)):

            U_LP = CauchyProblem(F, t, U0_LP_Sel, methods[j]) #ordenados tb

            #fig, (ax1, ax2) = plt.subplots(1, 2)
            plt.plot(U_LP[:,0], U_LP[:,1],'-',color = "k", label = 'Orbit')
            plt.plot(-mu, 0, 'o', color = "g", label = 'Tierra')
            plt.plot(1-mu, 0, 'o', color = "b", label = 'Luna')
            for i in range(5):
                plt.plot(LP[i,0],LP[i,1],ShapeLP[i],color = ColorLP[i],label=labelPTot[i])
            plt.xlim(-2,2)
            plt.ylim(-2,2)
            plt.title(f"Simulación CR3BP Tierra-Luna con esquema {methods[j].__name__}. Órbita en {labelP}. t = {tf}s, dt = {dt}. Vista completa" )    
            plt.legend(loc = 'upper left',bbox_to_anchor=(1., 0.95))
            plt.xlabel('X')
            plt.ylabel('Y')
            plt.grid()
            plt.show()
                    
            plt.plot(U_LP[:,0], U_LP[:,1],'-',color = "k", label = "Orbit" )
            plt.plot(LP[sel - 1,0], LP[sel - 1,1] , ShapeLP[sel-1],color = ColorLP[sel-1], label = labelPTot[sel-1])
            plt.title(f"Simulación CR3BP Tierra-Luna con esquema {methods[j].__name__}. Detalle de órbita en {labelP}. t = {tf}s, dt = {dt}" )
            plt.xlabel('X')
            plt.ylabel('Y')
            plt.legend(loc = 'upper right',bbox_to_anchor=(1, 0.5))
            plt.grid()   
            plt.xlim(LP[sel - 1,0]-0.2,LP[sel - 1,0]+0.2)
            plt.ylim(LP[sel - 1,1]-0.2,LP[sel - 1,1]+0.2)
            plt.legend(loc = 'upper left',bbox_to_anchor=(1., 0.95))
            plt.show()
                
            #manager = plt.get_current_fig_manager()
            #manager.full_screen_toggle()
            #figure = plt.gcf()                  
            #figure.set_size_inches(16, 8)       
            #plt.savefig('Plots/Hito 6/ CR3BP ' + label +' '+ methods[j].__name__ +'.png', bbox_inches = 'tight')
            #plt.close('all')
            

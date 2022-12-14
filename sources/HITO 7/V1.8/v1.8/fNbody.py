import numpy as np
import math
from numba import cuda

G = 6.67e-11 # Nm^2/kg
Nc = 3 # nº coord
# U --> vector estado
# r,v --> posicion, velocidad 
# dr(i)/dt = v(i)
# dv(i)/dt = sum(j) [G*m(j)*(r(j) - r(i)) / |r(j) - r(i)|^3]

# PARALELIZACION EMPLEANDO CUDA
@cuda.jit
def f(M_d,r_d_x,r_d_y,r_d_z,v_d_x,v_d_y,v_d_z,dr_d_x,dr_d_y,dr_d_z,dv_d_x,dv_d_y,dv_d_z):
  
  i = cuda.grid(1)
  size = len(r_d_x)

  if i < size:
    dr_d_x[i] = v_d_x[i]
    dr_d_y[i] = v_d_y[i]
    dr_d_z[i] = v_d_z[i]

    for j in range(size):
      if j!=i:
        dist_x = r_d_x[j] - r_d_x[i]
        dist_y = r_d_y[j] - r_d_y[i]
        dist_z = r_d_z[j] - r_d_z[i]

        norm_dist = math.sqrt(dist_x**2+dist_y**2+dist_z**2)
        
        # Descomentar si se quieren emplear las masas de los cuerpos para la variación de la aceleración
        dv_d_x[i] += (dist_x)/(norm_dist**3)*M_d[j]*G
        dv_d_y[i] += (dist_y)/(norm_dist**3)*M_d[j]*G
        dv_d_z[i] += (dist_z)/(norm_dist**3)*M_d[j]*G 


def F_NBody(U,t,Nb, Nc, M):

    '''
    Inputs:
            U = array con valores de posicion y velocidad para cada instante de tiempo
            t = tiempo
    Outputs: 
            F = array con valores de velocidad y aceleracion actualizadas con el tiempo
    '''
    
    U_sol = np.reshape(U,(Nb,2,Nc)) #Reshape de la U: dividir en arrays para cada cuerpo (i), cada cuerpo tiene en fila 0 la r(i) y en fila 1 la v(i)
    F = np.zeros(len(U), dtype=np.float64)
    dU_sol = np.reshape(F, (Nb, 2, Nc)) # Los valores introducidos en dU_sol serán insertados en F

    # RESHAPE DE POSICION Y VELOCIDAD
    r = np.reshape(U_sol[:, 0, :], (Nb, Nc)) # Guarda en array posiciones cada uno de los r(i) --> POSICIONES
    v = np.reshape(U_sol[:, 1, :], (Nb, Nc)) # Guarda en array velocidades cada uno de los v(i) --> VELOCIDADES
    # Separacion por coordenadas
    r_x = r[:,0]
    r_y = r[:,1]
    r_z = r[:,2]

    v_x = v[:,0]
    v_y = v[:,1]
    v_z = v[:,2]

    # RESHAPE DE VELOCIDAD Y ACELERACION
    drdt = np.reshape(dU_sol[:, 0, :], (Nb, Nc)) # Guarda en array derivada posiciones cada uno de los dr(i)/dt --> VELOCIDADES
    dvdt = np.reshape(dU_sol[:, 1, :], (Nb, Nc)) # Guarda en array derivada velocidades cada uno de los dv(i)/dt --> ACELERACIONES

    drdt_x = drdt[:,0]
    drdt_y = drdt[:,1]
    drdt_z = drdt[:,2]

    dvdt_x = dvdt[:,0]
    dvdt_y = dvdt[:,1]
    dvdt_z = dvdt[:,2]

    r_d_x = cuda.to_device(np.ascontiguousarray(r_x, dtype=np.float64))
    v_d_x = cuda.to_device(np.ascontiguousarray(v_x, dtype=np.float64))
    dr_d_x = cuda.to_device(np.ascontiguousarray(drdt_x, dtype=np.float64))
    dv_d_x = cuda.to_device(np.ascontiguousarray(dvdt_x, dtype=np.float64))

    r_d_y = cuda.to_device(np.ascontiguousarray(r_y, dtype=np.float64))
    v_d_y = cuda.to_device(np.ascontiguousarray(v_y, dtype=np.float64))
    dr_d_y = cuda.to_device(np.ascontiguousarray(drdt_y, dtype=np.float64))
    dv_d_y = cuda.to_device(np.ascontiguousarray(dvdt_y, dtype=np.float64))

    r_d_z = cuda.to_device(np.ascontiguousarray(r_z, dtype=np.float64))
    v_d_z = cuda.to_device(np.ascontiguousarray(v_z, dtype=np.float64))
    dr_d_z = cuda.to_device(np.ascontiguousarray(drdt_z, dtype=np.float64))
    dv_d_z = cuda.to_device(np.ascontiguousarray(dvdt_z, dtype=np.float64))

    M_d = cuda.to_device(np.ascontiguousarray(M, dtype=np.float64))
    # PARALELIZACION PARA CADA CUERPO
    nthreads = 1024
    nblocks = (len(r_d_x)//nthreads) + 1
    #print('Number of blocks: ',nblocks)
    f[nblocks,nthreads](M_d,r_d_x,r_d_y,r_d_z,v_d_x,v_d_y,v_d_z,dr_d_x,dr_d_y,dr_d_z,dv_d_x,dv_d_y,dv_d_z)

    drdt[:,0] = dr_d_x.copy_to_host()
    drdt[:,1] = dr_d_y.copy_to_host()
    drdt[:,2] = dr_d_z.copy_to_host()
    dvdt[:,0] = dv_d_x.copy_to_host()
    dvdt[:,1] = dv_d_y.copy_to_host()
    dvdt[:,2] = dv_d_z.copy_to_host()

    return F

import numpy as np
import math
from numba import cuda
from scipy.optimize import newton
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

from EDO import CauchyProblem
#Cauchy(F,t,U0,Scheme)
from TEMP_SCHEMES import Euler_Scheme, CrankNicolson_Scheme, RK4_Scheme, InvEuler_Scheme, LeapFrog, Embedded_RK
#Scheme(U,dt,t,F)
from fNbody import F_NBody #(U,t,Nb,Nc,M)
from initialConditions import initial_cond
from Physics import Kepler

#-------------------------------------
gpu = cuda.get_current_device()
print("name = %s" % gpu.name)
print("maxThreadsPerBlock = %s" % str(gpu.MAX_THREADS_PER_BLOCK))
print("maxBlockDimX = %s" % str(gpu.MAX_BLOCK_DIM_X))
print("maxBlockDimY = %s" % str(gpu.MAX_BLOCK_DIM_Y))
print("maxBlockDimZ = %s" % str(gpu.MAX_BLOCK_DIM_Z))
print("maxGridDimX = %s" % str(gpu.MAX_GRID_DIM_X))
print("maxGridDimY = %s" % str(gpu.MAX_GRID_DIM_Y))
print("maxGridDimZ = %s" % str(gpu.MAX_GRID_DIM_Z))
print("maxSharedMemoryPerBlock = %s" % str(gpu.MAX_SHARED_MEMORY_PER_BLOCK))
print("asyncEngineCount = %s" % str(gpu.ASYNC_ENGINE_COUNT))
print("canMapHostMemory = %s" % str(gpu.CAN_MAP_HOST_MEMORY))
print("multiProcessorCount = %s" % str(gpu.MULTIPROCESSOR_COUNT))
print("maxThreadsPerMultiprocessor = %s" % str(gpu.MAX_THREADS_PER_MULTI_PROCESSOR))
print("warpSize = %s" % str(gpu.WARP_SIZE))
print("unifiedAddressing = %s" % str(gpu.UNIFIED_ADDRESSING))
print("pciBusID = %s" % str(gpu.PCI_BUS_ID))
print("pciDeviceID = %s" % str(gpu.PCI_DEVICE_ID))

cuda.detect()
#------------------------------------------

print('Enter number of bodies:')
Nb = int(input()) # 1024*72*2 # nº cuerpos tiene que ser 1024 (nº threads por bloque) * 72 (nº bloques o SM) * A (nº entero a elegir)
Nc = 3 # nº coord

t0 = 0 # timepo inicial
print('Enter time [s]:')
tf = int(input()) # tiempo final
N = 10*tf # number of time steps
t = np.linspace(t0, tf, N)

U0, M = initial_cond(Nb,Nc) #ok

def F(U,t):#wrapp de la F para no cambiar Cauchy ni los esquemas
  return F_NBody(U, t, Nb, Nc)

U = CauchyProblem(F, t, U0, Embedded_RK)
U_rec = np.reshape(U,(N, Nb, 2, Nc))
r_rec = np.reshape(U_rec[:,:,0,:],(N,Nb,Nc))

#---------------------------------------------------

# 3D PLOT ANITAMED
def update_lines(num, walks, lines):
    for line, walk in zip(lines, walks):
        # NOTE: there is no .set_data() for 3 dim data...
        line.set_data(walk[:num, :2].T)
        line.set_3d_properties(walk[:num, 2])
    return lines

def paint(N,Nc,r,index):
  line=np.zeros((N,Nc))
  for i in range(N):
    line[i,:]=r[i,index,:]
  return line

# Data: Nb lines as (num_steps, 3) arrays
datas = [paint(N,Nc,r_rec,index) for index in range(Nb)]
#print ("datas: ",datas)

# Attaching 3D axis to the figure
fig = plt.figure()
ax = fig.add_subplot(projection="3d")

# Create lines initially without data
lines = [ax.plot([], [], [])[0] for _ in datas]

# Setting the axes properties
ax.set(xlim3d=(-10, 10), xlabel='X')
ax.set(ylim3d=(-10, 10), ylabel='Y')
ax.set(zlim3d=(-10, 10), zlabel='Z')

# Creating the Animation object
line_ani = animation.FuncAnimation(fig, update_lines, N//2, fargs=(datas, lines), interval=2500)

f = r"/content/drive/MyDrive/" + path + "/animation2.mp4" 
writervideo = animation.FFMpegWriter(fps=60) 
line_ani.save(f, writer=writervideo, dpi=300)
#plt.show()

#---------------------------------------------------

# 3D PLOT FIJO
colors = ['b','r','g','m','y','c']
fig = plt.figure(2)
ax = fig.add_subplot(projection='3d')
for i in range(Nb):
  ax.plot(r_rec[:, i, 0], r_rec[:, i, 1], r_rec[:, i, 2])#, colors[i])

plt.title(f'N = {Nb} body problem: 3D projection')
plt.xlabel("x")
plt.ylabel("y")

#plt.grid(True)
plt.show()
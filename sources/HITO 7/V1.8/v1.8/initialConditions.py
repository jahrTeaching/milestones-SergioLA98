import numpy as np

def random_vector(rand=np.random):
    '''
    Devuelve un vector unitario de coordenadas XYZ respecto ejes referencia
    '''
    phi = rand.uniform(0,2*np.pi)
    theta = rand.uniform(0,2*np.pi)
    x = np.cos(theta)*np.sin(phi)
    y = np.sin(theta)*np.sin(phi)
    z = np.cos(phi)
    vec = np.array([x,y,z])
    return vec

def random_mass(rand=np.random):
    '''
    Devuelve la masa del cuerpo [kg]
    '''
    return rand.uniform(4e5,2e30)

def random_dist(rand=np.random):
    '''
    Devuelve la distancia del cuerpo respecto a ejes referencia [m]
    '''
    return rand.uniform(46e9,6e12) # [m]

def random_veloc(rand=np.random):
    '''
    Devuelve la velocidad del cuerpo respecto a ejes referencia [m/s]
    '''
    return rand.uniform(1e3,1e5) # [m/s]

def initial_cond(Nb,Nc):
    # Descomentar el random_dist y random_veloc para aplicarle valores reales
    M = np.zeros(Nb)
    U0 = np.zeros(Nb*2*Nc)
    U_0 = np.reshape(U0, (Nb, 2, Nc) )
    r0 = np.reshape(U_0[:, 0, :], (Nb, Nc))
    v0 = np.reshape(U_0[:, 1, :], (Nb, Nc))

    for k in range(Nb):
        r0[k,:] = random_vector()*random_dist()
        v0[k,:] = random_vector()*random_veloc()
        M[k] = random_mass()

    return U0,M

def solarsystem():

    # SUN, MERCURY, VENUS, EARTH, MOON, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO 
    M = np.array([1988500., 0.330, 4.87, 5.97, 0.073, 0.642, 1898., 568., 86.8, 102., 0.0130])*1e+24

    U0 = np.zeros(11*2*3)
    U_0 = np.reshape(U0, (11, 2, 3) ) # punteros
    r0 = np.reshape(U_0[:, 0, :], (11, 3))
    v0 = np.reshape(U_0[:, 1, :], (11, 3))

    with open('Initialposition_solarsystem.txt') as data:
        substrings = data.read().split()
        p = [np.float64(substring) for substring in substrings]

    p = np.reshape(p, (11, 3))
    r0[:,:] = p[:,:]
    
    with open('Initialvelocity_solarsystem.txt') as data:
        substrings = data.read().split()
        v = [np.float(substring) for substring in substrings]

    v = np.reshape(v, (11, 3))
    v0[:,:] = v[:,:]

    return U0,M
from numpy import array

def LinearOscilator(U,t):
    r = U[0]
    drdt = U[1]
    return array([drdt, -r])
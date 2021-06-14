import numpy as np
from scipy.integrate import RK45


def Lorenz96(x0, F):
    x = np.zeros(x0.shape[0]+3)
    # 周期的なxの表現をつくる
    x[2:-1] = x0 # index 2がもとのxのindex 0, index -2がもとのxのN
    x[:2] = x0[-2:] # index 0,1 がもとのxのindex N-1, N
    x[-1] = x0[0] # index -1 がもとのxのindex 0
    
    dxdt = (x[3:] - x[:-3]) * x[1:-2] - x[2:-1] + F
    return dxdt

def L96(x, F):
    """
    Lorenz 96 model with constant forcing.
    Cited by "https://en.wikipedia.org/wiki/Lorenz_96_model"
    
    Parameters
    ----------
    x : 
        variables
    N : int
        number of sites
    """
    # Setting up vector
    x_ = np.zeros(x.shape)
    N = len(x)
    # Loops over indices (with operations and Python underflow indexing handling edge cases)
    for i in range(N):
        x_[i] = (x[(i + 1) % N] - x[i - 2]) * x[i - 1] - x[i] + F
    return x_
import numpy as np


def RungeKutta4(initial, time, model, F):
    """
    Parameters
    ----------
    F : 
        Forcing constant, variable used in L96
    kwargs : 
        modelに渡す x 以外の変数を指定
    """
    dt = time[1] - time[0]
    states = [initial]
    x = initial
    for t in time[:-1]:
        k1 = model(x, F)
        x1 = x + k1 * dt/2
        k2 = model(x1, F)
        x2 = x + k2 * dt/2
        k3 = model(x2, F)
        x3 = x + k3 * dt
        k4 = model(x3, F)
        x = x + (k1 + 2*k2 + 2*k3 + k4) * dt / 6
        states.append(x)
    states = np.stack(states)
    return states

def RK4(x, dt, model, F):
    """
    Parameters
    ----------
    F : 
        Forcing constant, variable used in L96
    """
    k1 = model(x, F)
    x1 = x + k1 * dt/2
    k2 = model(x1, F)
    x2 = x + k2 * dt/2
    k3 = model(x2, F)
    x3 = x + k3 * dt
    k4 = model(x3, F)
    return x + (k1 + 2*k2 + 2*k3 + k4) * dt / 6
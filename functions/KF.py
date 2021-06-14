import numpy as np
from OdeInt import RK4

I = np.identity(N, dtype=np.float)
H = np.identity(N, dtype=np.float)
R = np.identity(N, dtype=np.float)

def tangentM(x_a, delta=1e-8, dt=0.05, N=40):
    """
    M（RK4)の一次線形近似計算
    
    Parameter
    ---------
    delta : 
        Mの接線形近似のパラメタ
    x_a : 
    dt :
        Mの接線形近似のパラメタ
    N :
        観測点の個数
    """
    TLM = np.zeros((N, N), dtype=np.float_)
    for i in range(N):
        e_i = np.zeros(N, dtype=np.float_)
        e_i[i] = 1
        TLM[:, i] = (RK4(x_a + delta * e_i, dt, L96, 8) - RK4(x_a, dt, L96, 8)) / delta
    return TLM

def one_shotKF(x_a, P_a, y_o, alpha=1.1, delta=1e-8, dt=0.05, I=I, H=H, R=R):
    """
    KFの一回単位のシミュレーション。
    
    Parameter
    ---------
    x_a : 
    P_a :
    y_o :
    alpha :
        共分散膨張率
    delta : 
        Mの接線形近似のパラメタ
    dt : 0.05
        RK4のtime step。
    I :
    H :
    R :
    
    Return
    ------
    x_a, P_a
    """
    M = tangentM(x_a, delta)
    x_f = RK4(x_a, dt, model=L96, F=8)
    P_f = M * P_a * M.T
    P_f = P_f * alpha
    K = P_f * H.T * np.linalg.inv(H * P_f * H.T + R)
    x_a = x_f + K @ (y_o - H @ x_f)
    P_a = (I - K @ H) @ P_f
    return x_a, P_a
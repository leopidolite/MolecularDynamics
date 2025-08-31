import numpy as np 
import matplotlib.pyplot as plt
 
N = 300
rho = 0.25
L = np.sqrt(N/rho)
L2 = L/2
def periodic(sep):
    return np.where(sep > L / 2, sep - L,
           np.where(sep < -L / 2, sep + L, sep)) 

def unwrap_positions(x_mod):
    dx = np.diff(x_mod, axis=0)                    
    dx = np.where(dx >  L/2, dx - L, dx)
    dx = np.where(dx < -L/2, dx + L, dx)
    x_unwrapped = np.zeros_like(x_mod)
    x_unwrapped[0] = x_mod[0]                     
    x_unwrapped[1:] = x_unwrapped[0] + np.cumsum(dx, axis=0)
    return x_unwrapped

def precompute_neighbors(XX, YY, cage):
    T, N = XX.shape
    all_nbs = []
    c2 = cage*cage
    for t0 in range(T):
        nbs_t0 = []
        for i in range(N):
            dx = periodic(XX[t0, i] - XX[t0, :])
            dy = periodic(YY[t0, i] - YY[t0, :])
            dist2 = dx*dx + dy*dy
            nbs = np.where((dist2 <= c2) & (dist2 > 0))[0]
            nbs_t0.append(nbs)
        all_nbs.append(nbs_t0)
    return all_nbs

def MSD_cr_fast(X, Y, XX, YY, cage=2.0):
    T, N = X.shape
    msd = np.zeros(T-1, dtype=float)
    all_nbs = precompute_neighbors(XX, YY, cage)

    for dt in range(1, 1000):
        accum = 0.0
        count = 0
        for t0 in range(100):
            dX = X[t0+dt] - X[t0]
            dY = Y[t0+dt] - Y[t0]
            for i in range(N):
                nbs = all_nbs[t0][i]
                if nbs.size == 0:
                    continue
                nx = dX[nbs].mean()
                ny = dY[nbs].mean()
                dx_cr = dX[i] - nx
                dy_cr = dY[i] - ny
                accum += dx_cr*dx_cr + dy_cr*dy_cr
                count += 1
        msd[dt-1] = accum / count if count else np.nan
    return msd


temps = ['0.45','0.41', '0.37, 0.35', '0.25', '0.17']
for temp in temps:
    x_COM = np.load('/users/lli190/MolecularDynamics/300npy/300_'+temp+'T_x_COM.npy')
    y_COM = np.load('/users/lli190/MolecularDynamics/300npy/300_'+temp+'T_y_COM.npy')
    x_ = unwrap_positions(x_COM)[0:1200, :]
    y_ = unwrap_positions(y_COM)[0:1200, :]
    x_COM_ = x_COM[:, :]
    y_COM_ = y_COM[:, :]
    MSD = MSD_cr_fast(x_,y_,x_COM_, y_COM_)
    np.save('/users/lli190/MolecularDynamics/MSD_CR/MSD10_'+temp+'.npy', MSD)

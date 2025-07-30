import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft

N =108
rho = 0.25
L = np.sqrt(N/rho)

temps = ['047', '045', '043', '041', '039', '037', '035', '033', '031', '029', '027', '025', '023']

def load_file(filename, N):
    data = np.fromfile(filename, dtype=np.float64)
    assert data.size % N == 0, "File size does not match expected number of molecules"
    num_frames = data.size // N
    return data.reshape((num_frames, N))  # shape: (frames, molecules)

def full_analysis(temps):
    T_ = []
    PE_m = []
    s6 = []
    s5 = []
    s4 = []
    s3 = []
    s2 = []
    s6SD = []
    s5SD = []
    s4SD = []
    s3SD = []
    s2SD = []
    full_s6 = []

    orientational_timecorrs = []
    MSDs = []
    relaxation_times = []

    for temp in temps:
        prefix = '108long_' + temp
        read_folder_name = '/users/lli190/scratch/MD_analyzed_trajectories/saved_outputs/'

        xname= read_folder_name + prefix + '_x_COM.bin'
        yname= read_folder_name + prefix + '_y_COM.bin'
        thetaname = read_folder_name + prefix + '_theta.bin'
        with open('/users/lli190/scratch/MD_full_trajectories/Outputs/'+ prefix +'_KE.txt', 'r') as file: 
            KE = [float(line.rstrip()) for line in file]
        with open('/users/lli190/scratch/MD_full_trajectories/Outputs/'+ prefix +'_PE.txt', 'r') as file: 
            PE = [float(line.rstrip()) for line in file]

        ######## Load Data
        x_COM = load_file(xname, N)
        y_COM = load_file(yname, N)
        theta = load_file(thetaname, N)
        KE = np.array(KE)
        PE = np.array(PE)

        ######## Compute Temperature 
        T = np.mean(KE[0:])*2/3/N
        T_.append(T)
        PE_mean = np.mean(PE[0:])/N
        PE_m.append(PE_mean)

        ######## Compute Order Parameters
        def get_s_m(m):
            theta_j = theta
            s_m = np.abs(np.sum((np.exp(1j * m * theta_j).real), axis=1) / N)
            return s_m
        ### Full-run order parameters for export 
        def full_s_m(m):
            sm = get_s_m(m)
            smmean = np.mean(sm)
            s_m_var = np.sum((sm - smmean)**2)/(sm.size-1)
            s_m_SD = (s_m_var**0.5)
            return smmean, s_m_SD
        def get_s_m(m, theta):
            theta_j = theta
            s_m = np.abs(np.sum((np.exp(1j * m * theta_j).real), axis=1) / N)
            return s_m

        ### Full-run order parameters for export 
        def full_s_m(m, theta):
            sm = get_s_m(m, theta[0:])
            smmean = np.mean(sm)
            s_m_var = np.sum((sm - smmean)**2)/(sm.size-1)
            s_m_SD = (s_m_var**0.5)
            return smmean, s_m_SD, sm
        m6, m6SD, full_s_6 = full_s_m(6, theta)
        m5, m5SD,_ = full_s_m(5, theta)
        m4, m4SD,_ = full_s_m(4, theta)
        m3, m3SD,_ = full_s_m(3, theta)
        m2, m2SD,_ = full_s_m(2, theta)
        s6.append(m6)
        s6SD.append(m6SD)
        s5.append(m5)
        s5SD.append(m5SD)
        s4.append(m4)
        s4SD.append(m4SD)
        s3.append(m3)
        s3SD.append(m3SD)
        s2.append(m2)
        s2SD.append(m2SD)
        full_s6.append(full_s_6)

        ##### Time Correlation Functions
        def orientational_correlation(m):
            length = theta[0:,0].size
            C_total = np.zeros(length)
            for i in range(N):
                A = theta[0:,i]
                A = np.append(np.exp(1j*m*A), np.zeros(length))
                C_ = fft(A)
                C = ifft(C_*np.conj(C_)).real[0:length]
                norm = np.arange(C.size, 0, -1)
                C = C / norm
                C_total += C

            return C_total/N
        corr_m6 = orientational_correlation(6)
        orientational_timecorrs.append(corr_m6)

        time = np.linspace(0, corr_m6.size*0.001*10, corr_m6.size)
        def simpsons_integral(t, y):
            t = np.asarray(t)
            y = np.asarray(y)
            n = len(t) - 1

            h = t[1] - t[0]
            return (h / 3) * (
                y[0] + y[-1] +
                4 * np.sum(y[1:n:2]) +
                2 * np.sum(y[2:n-1:2])
            )
        tau_relax = simpsons_integral(time[:-2], corr_m6[:-2] - np.mean(corr_m6[700000:800000]))
        relaxation_times.append(tau_relax)

        ### Mean-Squared displacement
        def unwrap_positions(x_mod):
            dx = np.diff(x_mod, axis=0)                    
            dx = np.where(dx >  L/2, dx - L, dx)
            dx = np.where(dx < -L/2, dx + L, dx)
            x_unwrapped = np.zeros_like(x_mod)
            x_unwrapped[0] = x_mod[0]                     
            x_unwrapped[1:] = x_unwrapped[0] + np.cumsum(dx, axis=0)
            return x_unwrapped
        def autocorrFFT(x):
            N = len(x)
            F = fft(x, n=2*N)
            PSD = F * np.conj(F)
            res = ifft(PSD).real[:N]
            norm = N - np.arange(N)
            return res /norm

        def msd_fft(r):
            N = len(r)
            D = np.square(r).sum(axis=1)
            D = np.append(D, 0)
            S2 = sum(autocorrFFT(r[:, i]) for i in range(r.shape[1]))
            Q = 2 * D.sum()
            S1 = np.zeros(N)
            for m in range(N):
                Q -= D[m-1] + D[N-m]
                S1[m] = Q / (N - m)
            return S1 - 2 * S2
        def MSD(x, y):    
            msd_total = np.zeros(1000000)
            for i in range(x.shape[1]):
                r_i = np.stack((unwrap_positions(x[0:, i]), unwrap_positions(y[0:, i])), axis=1)  # shape (T, 2)
                msd_total += msd_fft(r_i)

            msd_avg = msd_total / x.shape[1]
            return msd_avg
        msd_t = MSD(x_COM, y_COM)
        MSDs.append(msd_t)

    ################################################MODIFY THIS#######################
    savepath = '/users/lli190/scratch/MD_analyzed_trajectories/saved_outputs/'
    ################################################MODIFY THIS#######################

    np.savez('/users/lli190/MolecularDynamics/108_analyzed.npz',
            T_=T_,
            PE_m=PE_m,
            s6=s6, s5=s5, s4=s4, s3=s3, s2=s2,
            s6SD=s6SD, s5SD=s5SD, s4SD=s4SD, s3SD=s3SD, s2SD=s2SD,
            orientational_timecorrs=orientational_timecorrs,
            MSDs=MSDs, full_s6 = full_s6,
            relaxation_times=relaxation_times)
    

full_analysis(temps)

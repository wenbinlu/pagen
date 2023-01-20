from scipy.interpolate import interp1d, RectBivariateSpline
import numpy as np
from math import log10, floor


def read_F_eta_nu(fdir, fname):   # for F2 and F3
    f = open(fdir + fname + '.txt', 'r')
    # print('read data from ' + fdir + fname + '.txt')
    f.seek(0)  # go to the beginning
    data = f.read().split('\n')  # all rows
    Neta_nu = len(data)
    logF_arr = np.zeros(Neta_nu)
    for i in range(Neta_nu):
        logF_arr[i] = data[i]
    f.close()
    return logF_arr


def read_IK_logeta_logthe(fdir, fname):   # include I1-I8, K1-K6
    f = open(fdir + fname + '.txt', 'r')
    f.seek(0)  # go to the beginning
    data = f.read().split('\n')  # all rows
    N1, N2 = len(data), len(data[0].split('\t'))
    logIK_arr = np.zeros((N1, N2))
    for i in range(N1):
        row = data[i].split('\t')
        logIK_arr[i, :] = row[:]
    f.close()
    return logIK_arr


def intp_F_logeta_nu(n, eta_nu_arr, Fdata):
    return interp1d(eta_nu_arr, Fdata[n-2][:])


def intp_I_logeta_logthe(n, logeta_arr, logthe_arr, Idata):
    return RectBivariateSpline(logeta_arr, logthe_arr, Idata[n-1][:],
                               kx=3, ky=3, s=0)   # s is smoothing parameter


def intp_K_logeta_logthe(n, logeta_arr, logthe_arr, Kdata):
    return RectBivariateSpline(logeta_arr, logthe_arr, Kdata[n-1][:],
                               kx=3, ky=3, s=0)   # s is smoothing parameter


# read only 2 tables
def read_data_K14(Ye, logYe_arr, logeta_arr, logthe_arr, Kdata, fdir):
    n_list = [1, 4]   # only these four tables are loaded
    # new prescription (linear interpolation)
    logYe = log10(Ye)
    dlogYe = logYe_arr[1] - logYe_arr[0]
    iYe1 = int(floor((logYe-logYe_arr[0])/dlogYe))
    iYe2 = int(iYe1 + 1)
    Ntab = len(n_list)
    Nlogeta, Nlogthe = len(logeta_arr), len(logthe_arr)
    Kdata1 = np.zeros((Ntab, Nlogeta, Nlogthe))    # for iYe1
    Kdata2 = np.zeros((Ntab, Nlogeta, Nlogthe))    # for iYe2
    logYe1 = logYe_arr[iYe1]
    for i in range(Ntab):
        n = n_list[i]
        fname = 'K%d_logYe_m%dp%04d'\
            % (n, floor(-logYe1), round((-logYe1) % 1*1e4))
        Kdata1[i, :] = read_IK_logeta_logthe(fdir, fname)
    logYe2 = logYe_arr[iYe2]
    for i in range(Ntab):
        n = n_list[i]
        fname = 'K%d_logYe_m%dp%04d'\
            % (n, floor(-logYe2), round((-logYe2) % 1*1e4))
        Kdata2[i, :] = read_IK_logeta_logthe(fdir, fname)
    # linear interpolation between Kdata1 and Kdata2
    for i in range(Ntab):
        n = n_list[i]
        slope = (Kdata2[i, :] - Kdata1[i, :])/dlogYe
        Kdata[n-1][:] = Kdata1[i, :] + slope*(logYe - logYe1)
    # del Kdata1, Kdata2, slope
    return None


# read only 4 tables
def read_data_K1245(Ye, logYe_arr, logeta_arr, logthe_arr, Kdata, fdir):
    n_list = [1, 2, 4, 5]   # only these four tables are loaded
    # new prescription (linear interpolation)
    logYe = log10(Ye)
    dlogYe = logYe_arr[1] - logYe_arr[0]
    iYe1 = int(floor((logYe-logYe_arr[0])/dlogYe))
    iYe2 = int(iYe1 + 1)
    Ntab = len(n_list)
    Nlogeta, Nlogthe = len(logeta_arr), len(logthe_arr)
    Kdata1 = np.zeros((Ntab, Nlogeta, Nlogthe))    # for iYe1
    Kdata2 = np.zeros((Ntab, Nlogeta, Nlogthe))    # for iYe2
    logYe1 = logYe_arr[iYe1]
    for i in range(Ntab):
        n = n_list[i]
        fname = 'K%d_logYe_m%dp%04d'\
            % (n, floor(-logYe1), round((-logYe1) % 1*1e4))
        Kdata1[i, :] = read_IK_logeta_logthe(fdir, fname)
    logYe2 = logYe_arr[iYe2]
    for i in range(Ntab):
        n = n_list[i]
        fname = 'K%d_logYe_m%dp%04d'\
            % (n, floor(-logYe2), round((-logYe2) % 1*1e4))
        Kdata2[i, :] = read_IK_logeta_logthe(fdir, fname)
    # linear interpolation between Kdata1 and Kdata2
    for i in range(Ntab):
        n = n_list[i]
        slope = (Kdata2[i, :] - Kdata1[i, :])/dlogYe
        Kdata[n-1][:] = Kdata1[i, :] + slope*(logYe - logYe1)
    # del Kdata1, Kdata2, slope
    return None


# read all 6 tables
def read_data_K123456(Ye, logYe_arr, logeta_arr, logthe_arr, Kdata, fdir):
    n_list = [1, 2, 3, 4, 5, 6]
    # new prescription (linear interpolation)
    logYe = log10(Ye)
    dlogYe = logYe_arr[1] - logYe_arr[0]
    iYe1 = int(floor((logYe-logYe_arr[0])/dlogYe))
    iYe2 = int(iYe1 + 1)
    Ntab = len(n_list)
    Nlogeta, Nlogthe = len(logeta_arr), len(logthe_arr)
    Kdata1 = np.zeros((Ntab, Nlogeta, Nlogthe))    # for iYe1
    Kdata2 = np.zeros((Ntab, Nlogeta, Nlogthe))    # for iYe2
    logYe1 = logYe_arr[iYe1]
    for i in range(Ntab):
        n = n_list[i]
        fname = 'K%d_logYe_m%dp%04d'\
            % (n, floor(-logYe1), round((-logYe1) % 1*1e4))
        Kdata1[i, :] = read_IK_logeta_logthe(fdir, fname)
    logYe2 = logYe_arr[iYe2]
    for i in range(Ntab):
        n = n_list[i]
        fname = 'K%d_logYe_m%dp%04d'\
            % (n, floor(-logYe2), round((-logYe2) % 1*1e4))
        Kdata2[i, :] = read_IK_logeta_logthe(fdir, fname)
    # linear interpolation between Kdata1 and Kdata2
    for i in range(Ntab):
        n = n_list[i]
        slope = (Kdata2[i, :] - Kdata1[i, :])/dlogYe
        Kdata[n-1][:] = Kdata1[i, :] + slope*(logYe - logYe1)
    # del Kdata1, Kdata2, slope
    return None

from parameters import *
from read_intp_funcs import read_F_eta_nu, read_IK_logeta_logthe, \
    intp_F_logeta_nu, intp_I_logeta_logthe
from data_init import *

for n in [1, 2, 3, 4, 5, 6, 7, 8]:
    fname = 'I%d' % n
    Idata[n-1][:] = read_IK_logeta_logthe(fdir, fname)

for n in [2, 3]:
    fname = 'F%d' % n
    Fdata[n-2][:] = read_F_eta_nu(fdir, fname)


# ----- interpolate all data
intp_I1 = intp_I_logeta_logthe(1, logeta_arr, logthe_arr, Idata)  # pair pressure
intp_I2 = intp_I_logeta_logthe(2, logeta_arr, logthe_arr, Idata)  # pair mass density
intp_I3 = intp_I_logeta_logthe(3, logeta_arr, logthe_arr, Idata)  # pair entropy
intp_I4 = intp_I_logeta_logthe(4, logeta_arr, logthe_arr, Idata)  # q_nu
intp_I5 = intp_I_logeta_logthe(5, logeta_arr, logthe_arr, Idata)  # q_bnu
intp_I6 = intp_I_logeta_logthe(6, logeta_arr, logthe_arr, Idata)  # lam_ep
intp_I7 = intp_I_logeta_logthe(7, logeta_arr, logthe_arr, Idata)  # lam_en
intp_I8 = intp_I_logeta_logthe(8, logeta_arr, logthe_arr, Idata)  # pair thermal energy density

intp_F2 = intp_F_logeta_nu(2, eta_nu_arr, Fdata)    # n_nu or n_bnu
intp_F3 = intp_F_logeta_nu(3, eta_nu_arr, Fdata)    # U_nu or U_bnu

# ---- below is outdated not used
# read_data_K123456(Ye0)     # load K123456 data for the initial Ye0
# intp_K1 = intp_K_logeta_logthe(1)   # kap_nu
# intp_K2 = intp_K_logeta_logthe(2)   # kap_t_nu
# intp_K3 = intp_K_logeta_logthe(3)   # kap_nu_P
# intp_K4 = intp_K_logeta_logthe(4)   # kap_bnu
# intp_K5 = intp_K_logeta_logthe(5)   # kap_t_bnu
# intp_K6 = intp_K_logeta_logthe(6)  # kap_bnu_P
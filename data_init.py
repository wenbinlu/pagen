import numpy as np

# grid of (logeta, logthe), fixed for all FD_integs
Nlogeta, Nlogthe = 300, 100
logeta_arr = np.linspace(-4.5, 1.5, Nlogeta)     # eta = (mu-1)/the
# note: make sure eta_max > x_max (in 'generate_numeric_data')
logthe_arr = np.linspace(-0.2, 2.0, Nlogthe)      # the = kT/mec^2
NlogYe = 200    # only used for FYe_eta_logthe
logYe_arr = np.linspace(-1.7, -0.28, NlogYe)

Neta_nu = 10000  # for nu or bnu
eta_nu_arr = np.linspace(-100, 100, Neta_nu)  # +/-(mu_nu/the)

I1 = np.zeros((Nlogeta, Nlogthe))   # Ppair*V/mec2
I2 = np.zeros((Nlogeta, Nlogthe))   # Ye*rho*V/mp
I3 = np.zeros((Nlogeta, Nlogthe))   # spair*rho*V/k_B
I4 = np.zeros((Nlogeta, Nlogthe))   # q_nu,thin
I5 = np.zeros((Nlogeta, Nlogthe))   # q_bnu,thin
I6 = np.zeros((Nlogeta, Nlogthe))   # lam_nu,thin
I7 = np.zeros((Nlogeta, Nlogthe))   # lam_bnu,thin
# add more 2D arrays
I8 = np.zeros((Nlogeta, Nlogthe))   # U_pairs*V/mec2
Idata = [I1, I2, I3, I4, I5, I6, I7, I8]

F2 = np.zeros(Neta_nu)  # n_nu*V/the**3
F3 = np.zeros(Neta_nu)  # U_nu*V/(the**4*mec2)
Fdata = [F2, F3]

# ---- used for for each Ye
K1 = np.zeros((Nlogeta, Nlogthe))  # kap_nu*mp/sigma0  (Rosseland mean for energy flux)
K2 = np.zeros((Nlogeta, Nlogthe))  # tkap_nu*mp/sigma0  (Rosseland mean for number flux)
K3 = np.zeros((Nlogeta, Nlogthe))  # kap_nu_P*mp/sigma0 (Planck mean)
K4 = np.zeros((Nlogeta, Nlogthe))  # kap_bnu*mp/sigma0  (Rosseland mean for energy flux)
K5 = np.zeros((Nlogeta, Nlogthe))  # tkap_bnu*mp/sigma0  (Rosseland mean for number flux)
K6 = np.zeros((Nlogeta, Nlogthe))  # kap_bnu_P*mp/sigma0 (Planck mean)
Kdata = [K1, K2, K3, K4, K5, K6]

# ---- below is not used

# read inner disk data
# print('read the inner disk neutrino emission data')
# fname = 'inner_disk_emission'
# f = open(fdir + fname + '.txt', 'r')
# f.seek(0)  # go to the beginning
# nbegin = 1
# data_all = f.readlines()[nbegin:]  # a list of all rows, skip the first nbegin rows
# NMdot = len(data_all)
# logMdot_arr = np.zeros(NMdot)
# logLnu_rin = np.zeros(NMdot)
# eps_nu_rin = np.zeros(NMdot)
# logLbnu_rin = np.zeros(NMdot)
# eps_bnu_rin = np.zeros(NMdot)
# for i in range(NMdot):
#     row = data_all[i].strip().split()
#     [logMdot_arr[i], logLnu_rin[i], eps_nu_rin[i], logLbnu_rin[i], eps_bnu_rin[i]]\
#         = [float(element) for element in row]
# f.close()

# interpolate to obtain the inner disk neutrino emission
# intp_logLnu = interp1d(logMdot_arr, logLnu_rin, fill_value='extrapolate')
# intp_eps_nu = interp1d(logMdot_arr, eps_nu_rin, fill_value='extrapolate')
# intp_logLbnu = interp1d(logMdot_arr, logLbnu_rin, fill_value='extrapolate')
# intp_eps_bnu = interp1d(logMdot_arr, eps_bnu_rin, fill_value='extrapolate')
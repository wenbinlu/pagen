from other_funcs import phi_p
from math import log10, exp
from data_init import *
from parameters import *
from scipy.interpolate import interp1d

F2 = np.zeros(Neta_nu)  # n_nu
F3 = np.zeros(Neta_nu)  # U_nu
F4 = np.zeros(Neta_nu)  # for <E_nu^2>, cross-section
F5 = np.zeros(Neta_nu)  # for <E_nu^3>, energy absorption rate
F6 = np.zeros(Neta_nu)  # for <E_nu^4>, correction to bnup energy absorption rate

# maximum particle energy considered (for pairs and neutrinos)
xmax = 5. + eta_nu_arr.max()
# obtain the function x(t) precisely using a fine grid
Nt = 100000
x_arr_fine = np.empty(Nt)
t_arr_fine = np.linspace(0, 1, Nt, endpoint=False)
dt_fine = t_arr_fine[1] - t_arr_fine[0]
t_arr_fine += dt_fine/2.  # mid point
temp = 0.
for i in range(Nt):
    t = t_arr_fine[i]
    temp += phi_p(t)*dt_fine
    x_arr_fine[i] = temp
x_t_max = x_arr_fine.max()
A = xmax/x_t_max
# print('A=%.7e' % A)
t_arr_fine += dt_fine/2.  # the right boundary of each bin
x_arr_fine *= A
xt_intp = interp1d(t_arr_fine, x_arr_fine)    # interpolation

# note: since a fixed grid in t is used for all integrals, we should simply store
#       all x(t) and phi_p(t) as two arrays
Nx = 200    # number of steps for the numerical FD integrals
x_arr = np.zeros(Nx)
phip_arr = np.zeros(Nx)
t_arr = np.linspace(0, 1, Nx, endpoint=False)
dt = t_arr[1]-t_arr[0]
t_arr += dt/2.  # mid-point
for i in range(Nx):
    t = t_arr[i]
    phip_arr[i] = phi_p(t)
    x_arr[i] = xt_intp(t)
del x_arr_fine, t_arr_fine  # no longer need these fine grid arrays
# the following integrals need to use Nx, dt, x_arr, phip_arr, A, x_t_max

#########################################
print('generating F2, F3, F4, F5, F6 data')
for i in range(Neta_nu):
    eta_nu = eta_nu_arr[i]
    Feta_nu2, Feta_nu3, Feta_nu4, Feta_nu5, Feta_nu6 = 0., 0., 0., 0., 0.
    for ix in range(Nx):
        x = x_arr[ix]
        phip = phip_arr[ix]
        a = x**2/(exp(x-eta_nu) + 1)*phip
        Feta_nu2 += a
        Feta_nu3 += a*x
        Feta_nu4 += a*x**2
        Feta_nu5 += a*x**3
        Feta_nu6 += a*x**4
    F2[i] = log10(Feta_nu2*A*dt)
    F3[i] = log10(Feta_nu3*A*dt)
    F4[i] = log10(Feta_nu4*A*dt)
    F5[i] = log10(Feta_nu5*A*dt)
    F6[i] = log10(Feta_nu6*A*dt)
print('writing F2, F3, F4, F5, F6 into files')
Fdata = [F2, F3, F4, F5, F6]
for n in [2, 3, 4, 5, 6]:
    savename = 'F%d' % n
    wrt_content = Fdata[n-2]
    f = open(fdir + savename + '.txt', 'w')
    f.seek(0)  # go to the beginning
    for i in range(Neta_nu):
        if i == 0:
            f.write('%.8e' % (wrt_content[i]))
        else:
            f.write('\n%.8e' % (wrt_content[i]))
    f.close()

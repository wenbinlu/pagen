from math import sqrt, log10, log, pi, exp, floor
from other_funcs import phi_p, f_mass
from data_init import *
from parameters import *
from read_intp_funcs import read_F_eta_nu
from scipy.interpolate import interp1d
import multiprocessing
from multiprocessing import Process

# NOTE: need to run 'gen_F_data.py' first (unless 'FN.txt' are already in fdir)

# maximum particle energy considered (for pairs and neutrinos)
xmax = max(15. + 10**(logthe_arr.max()), 40.)
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
Nx = 100    # no steps for the numerical integrals  (150 for better performance)
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

print('generating I1--I8 data')
k_arr = np.array([0.5, 1.5, 2.5, 3.5, 4.5])   # for all the Fk's
Nk = len(k_arr)
Fk_e_arr = np.zeros(Nk)  # FD integrals for electrons
Fk_inc_e_arr = np.zeros(Nk)  # incomplete FD integrals for electrons
Fk_p_arr = np.zeros(Nk)  # for positrons
Fk_inc_p_arr = np.zeros(Nk)  # for positrons

percent = 0
for i in range(Nlogeta):
    eta = 10**logeta_arr[i]
    if floor(i*100./Nlogeta) > percent:
        print('%d percent' % percent)
        percent += 10
    for j in range(Nlogthe):
        the = 10**logthe_arr[j]
        eta_p = -eta - 2./the   # for positrons
        for n in range(Nk):
            k = k_arr[n]
            Fk_e, Fk_inc_e = 0., 0.
            Fk_p, Fk_inc_p = 0., 0.
            x0 = (Q-1.)/the
            for ix in range(Nx):
                x = x_arr[ix]   # full integrals
                phip = phip_arr[ix]
                a = x**k*sqrt(1+0.5*the*x)*phip
                Fk_e += a/(exp(x-eta)+1)
                Fk_p += a/(exp(x-eta_p)+1)
                x1 = x + x0     # incomplete integrals
                a1 = x1**k*sqrt(1+0.5*the*x1)*phip
                Fk_inc_e += a1/(exp(x1-eta)+1)
                Fk_inc_p += a1/(exp(x1-eta_p)+1)
            # store the data
            Fk_e_arr[n] = Fk_e*A*dt
            Fk_inc_e_arr[n] = Fk_inc_e*A*dt
            Fk_p_arr[n] = Fk_p*A*dt
            Fk_inc_p_arr[n] = Fk_inc_p*A*dt
        # now all full/incomplete FD integrals have been obtained
        I1[i, j] = log10(1./3*2**1.5*the**2.5
                         * (Fk_e_arr[1] + Fk_p_arr[1]
                            + 0.5*the*(Fk_e_arr[2] + Fk_p_arr[2])
                            )
                         )
        I2[i, j] = log10(sqrt(2)*the**1.5
                         * (Fk_e_arr[0] - Fk_p_arr[0]
                            + the*(Fk_e_arr[1] - Fk_p_arr[1])
                            )
                         )
        I8[i, j] = log10(sqrt(2)*the**2.5
                         * (Fk_e_arr[1] + Fk_p_arr[1]
                            + the*(Fk_e_arr[2] + Fk_p_arr[2])
                            )
                         )  # only thermal, not including rest mass
        # I8[i, j] = log10(sqrt(2)*the**1.5*(Fk_e_arr[0] + 2*the*Fk_e_arr[1]
        #                                    + Fk_p_arr[0] + 2*the*Fk_p_arr[1]
        #                                    + the**2*(Fk_e_arr[2]+Fk_p_arr[2])
        #                                    )
        #                  )    # including rest mass
        I3[i, j] = log10(sqrt(2)*the**1.5
                         * (4./3*the*Fk_e_arr[2]
                            + (5./3-eta*the)*Fk_e_arr[1]
                            - eta*Fk_e_arr[0]
                            + 4./3*the*Fk_p_arr[2]
                            + (5./3-eta_p*the)*Fk_p_arr[1]
                            - eta_p*Fk_p_arr[0]
                            )
                         )
        I4[i, j] = log10(sqrt(2)*the**1.5
                         * (the**4*Fk_inc_e_arr[4]
                            - (3.*Q-4)*the**3*Fk_inc_e_arr[3]
                            + 3.*(Q-1)*(Q-2)*the**2*Fk_inc_e_arr[2]
                            + (4.-Q)*(Q-1)**2*the*Fk_inc_e_arr[1]
                            - (Q-1)**3*Fk_inc_e_arr[0]
                            )
                         )
        I5[i, j] = log10(sqrt(2)*the**1.5
                         * (the**4*Fk_p_arr[4]
                            + (3.*Q+4)*the**3*Fk_p_arr[3]
                            + 3.*(Q+1)*(Q+2)*the**2*Fk_p_arr[2]
                            + (4.+Q)*(Q+1)**2*the*Fk_p_arr[1]
                            + (Q+1)**3*Fk_p_arr[0]
                            )
                         )
        I6[i, j] = log10(sqrt(2)*the**1.5
                         * (the**3*Fk_inc_e_arr[3]
                            - (2.*Q-3)*the**2*Fk_inc_e_arr[2]
                            + (Q-1.)*(Q-3)*the*Fk_inc_e_arr[1]
                            + (Q-1.)**2*Fk_inc_e_arr[0]
                            )
                         )
        I7[i, j] = log10(sqrt(2)*the**1.5
                         * (the**3*Fk_p_arr[3]
                            + (2.*Q+3)*the**2*Fk_p_arr[2]
                            + (Q+1.)*(Q+3)*the*Fk_p_arr[1]
                            + (Q+1.)**2*Fk_p_arr[0]
                            )
                         )

print('writing I1--I8 data into files')
# the same logeta for each row
# the same logthe for each column
Idata = [I1, I2, I3, I4, I5, I6, I7, I8]
for n in [1, 2, 3, 4, 5, 6, 7, 8]:
    savename = 'I%d' % n
    wrt_content = Idata[n-1]
    f = open(fdir + savename + '.txt', 'w')
    f.seek(0)  # go to the beginning
    for i in range(Nlogeta):
        if i >= 1:
            f.write('\n')
        for j in range(Nlogthe):
            if j == 0:
                f.write('%.8e' % wrt_content[i, j])
            else:
                f.write('\t%.8e' % wrt_content[i, j])
    f.close()

# exit()  # stop here if we only need to generate I1-I8 data

# nu/bnu data for number and energy density
F2 = read_F_eta_nu(fdir, 'F2')
F3 = read_F_eta_nu(fdir, 'F3')
print('interpolate F2 and F3 data for KN calculations')
# interpolate (eta_nu_arr, F1) and (eta_nu_arr, F2)
F2_intp = interp1d(eta_nu_arr, F2)
F3_intp = interp1d(eta_nu_arr, F3)


# generate K1, K2, K3 opacity data for neutrinos
def gen_K123(nlist, Nlen):   # Nlen = len(nlist), not used
    K1_wrt = np.zeros((Nlogeta, Nlogthe))
    K2_wrt = np.zeros((Nlogeta, Nlogthe))
    K3_wrt = np.zeros((Nlogeta, Nlogthe))
    for n in nlist:
        Ye = 10**logYe_arr[n]
        for i in range(Nlogeta):
            eta = 10**logeta_arr[i]
            for j in range(Nlogthe):
                the = 10**(logthe_arr[j])
                mu = eta*the + 1.
                # need baryonic mass fractions Xp, Xn, Xa
                Xa_guess = 0.01
                rho = 10**I2[i, j]*cst.mp/(Ye*Vunit)
                Xp, Xn, Xa = f_mass(Ye, rho, the, Xa_guess)
                # chemical potential of neutrinos
                eta_nu = eta - (Q-1.)/the - log(Xn/Xp)
                if eta_nu > eta_nu_arr.max() or eta_nu < eta_nu_arr.min():  # test
                    print('eta_nu=%.3e outside the range!' % eta_nu)
                    print('eta=%.3e, the=%.3e, Xp=%.3e, Xn=%.3e, Ye=%.3e' %
                          (eta, the, Xp, Xn, Ye))
                    exit()
                Fnu_flux, Fnu_num, Fnu_flux_P = 0., 0., 0.  # 'Fnu_flux_P for Planck-mean'
                b = 0.30*Xp + 0.36*Xn + 0.21*Xa
                for ix in range(Nx):
                    x = x_arr[ix]
                    phip = phip_arr[ix]
                    a = x**2/(exp(x-eta_nu) + 1)*phip
                    temp = (x*the + Q - mu)/the
                    if temp > log(1. + 1./epsilon_small):
                        f_block = 1.    # no blocking effect
                    elif temp < log(epsilon_small):  # Taylor expansion
                        f_block = exp(temp)  # extremely strong blocking
                    else:
                        f_block = 1 - 1./(exp(temp) + 1)    # in the middle
                    c = 1.41*Xn*f_block*(x + Q/the)*sqrt((x+Q/the)**2 - 1/the**2)
                    d = a/(b*x**2 + c)
                    Fnu_flux += x*d
                    Fnu_num += d
                    Fnu_flux_P += x*a*(b*x**2 + c)
                Fnu_flux *= A*dt
                Fnu_num *= A*dt
                Fnu_flux_P *= A*dt
                K1_wrt[i, j] = log10(the**2*10**F3_intp(eta_nu)/Fnu_flux)
                K2_wrt[i, j] = log10(the**2*10**F2_intp(eta_nu)/Fnu_num)
                K3_wrt[i, j] = log10(the**2*Fnu_flux_P/10**F2_intp(eta_nu))
        # write data into files
        Kdata_wrt = [K1_wrt, K2_wrt, K3_wrt]
        for m in [1, 2, 3]:
            logYe = log10(Ye)
            savename = 'K%d_logYe_m%dp%04d'\
                       % (m, floor(-logYe), round((-logYe) % 1*1e4))
            wrt_content = Kdata_wrt[m-1]
            f = open(fdir + savename + '.txt', 'w')
            f.seek(0)  # go to the beginning
            for i in range(Nlogeta):
                if i >= 1:
                    f.write('\n')
                for j in range(Nlogthe):
                    if j == 0:
                        f.write('%.8e' % wrt_content[i, j])
                    else:
                        f.write('\t%.8e' % wrt_content[i, j])
            f.close()
    print('finished nlist = [%d, %d]' % (nlist[0], nlist[-1]))


# generate K4, K5, K6 opacity data for anti-neutrinos
def gen_K456(nlist, Nlen):  # Nlen = len(nlist), not used
    K4_wrt = np.zeros((Nlogeta, Nlogthe))
    K5_wrt = np.zeros((Nlogeta, Nlogthe))
    K6_wrt = np.zeros((Nlogeta, Nlogthe))
    for n in nlist:
        Ye = 10**logYe_arr[n]
        for i in range(Nlogeta):
            eta = 10**logeta_arr[i]
            for j in range(Nlogthe):
                the = 10**logthe_arr[j]
                mu = eta*the + 1.   # electron chemical potential
                # need baryonic mass fractions Xp, Xn, Xa
                Xa_guess = 0.01
                rho = 10**I2[i, j]*cst.mp/(Ye*Vunit)
                Xp, Xn, Xa = f_mass(Ye, rho, the, Xa_guess)
                # chemical potential of neutrinos
                eta_bnu = -(eta - (Q-1.)/the - log(Xn/Xp))
                b = 0.30*Xp + 0.36*Xn + 0.21*Xa     # scattering term
                x0 = (Q+1.)/the  # the break point
                A0 = x0/x_t_max
                # for the first segment of integral: 0 < x < x0 (no absorption)
                Fbnu_flux1, Fbnu_num1, Fbnu_flux1_P = 0., 0., 0.
                for ix in range(Nx):
                    x = x_arr[ix]
                    phip = phip_arr[ix]
                    a = 1./(exp(x-eta_bnu) + 1)*phip
                    Fbnu_flux1 += x*a/b
                    Fbnu_num1 = a/b
                    Fbnu_flux1_P += x**3 * a * b * x**2
                Fbnu_flux1 *= A0*dt
                Fbnu_num1 *= A0*dt
                Fbnu_flux1_P *= A0*dt
                # for the second segment of integral: x0 < x < xmax(~50)
                Fbnu_flux2, Fbnu_num2, Fbnu_flux2_P = 0., 0., 0.
                if x0 < xmax:   # this segment exists
                    for ix in range(Nx):
                        x = x_arr[ix] + x0
                        phip = phip_arr[ix]
                        a = 1./(exp(x-eta_bnu) + 1)*phip
                        d = (x-Q/the)**2 - 1/the**2
                        if d < epsilon_small:
                            d = 2*x_arr[ix]  # linear expansion for small (x-x0)
                        y = 0.51099*the*x   # energy of bnu in MeV
                        gfit = y**(-0.07056 + 0.02018*log(y) - 0.001953*(log(y))**3)
                        temp = (x*the - Q + mu)/the
                        if temp > log(1. + 1./epsilon_small):
                            f_block = 1.    # no blocking effect
                        elif temp < log(epsilon_small):  # Taylor expansion
                            f_block = exp(temp)  # extremely strong blocking
                        else:
                            f_block = 1 - 1./(exp(temp) + 1)    # in the middle
                        c = 1.48*gfit*Xp*f_block*(x-Q/the)*sqrt(d)/x**2   # absorption
                        Fbnu_flux2 += x*a/(b + c)
                        Fbnu_num2 += a/(b + c)
                        Fbnu_flux2_P += x**3 * a * (b + c) * x**2
                    Fbnu_flux2 *= A*dt
                    Fbnu_num2 *= A*dt
                    Fbnu_flux2_P *= A*dt
                Fbnu_flux = Fbnu_flux1 + Fbnu_flux2
                Fbnu_num = Fbnu_num1 + Fbnu_num2
                Fbnu_flux_P = Fbnu_flux1_P + Fbnu_flux2_P
                K4_wrt[i, j] = log10(the**2*10**F3_intp(eta_bnu)/Fbnu_flux)
                K5_wrt[i, j] = log10(the**2*10**F2_intp(eta_bnu)/Fbnu_num)
                K6_wrt[i, j] = log10(the**2*Fbnu_flux_P/10**F3_intp(eta_bnu))
        # write data into files
        Kdata_wrt = [K4_wrt, K5_wrt, K6_wrt]
        for m in [4, 5, 6]:
            logYe = log10(Ye)
            savename = 'K%d_logYe_m%dp%04d'\
                       % (m, floor(-logYe), round((-logYe) % 1*1e4))
            wrt_content = Kdata_wrt[m-4]
            f = open(fdir + savename + '.txt', 'w')
            f.seek(0)  # go to the beginning
            for i in range(Nlogeta):
                if i >= 1:
                    f.write('\n')
                for j in range(Nlogthe):
                    if j == 0:
                        f.write('%.8e' % wrt_content[i, j])
                    else:
                        f.write('\t%.8e' % wrt_content[i, j])
            f.close()
    print('finished nlist = [%d, %d]' % (nlist[0], nlist[-1]))


Ncpu = multiprocessing.cpu_count()
print('number of CPUs', Ncpu)
# divide the task into Ncpu chunks
nlist_chunks = np.array_split(range(NlogYe), Ncpu)

print('generating K1, K2, K3 data')
procs = [Process(target=gen_K123,
                 args=(nlist_chunks[n], len(nlist_chunks[n])))
         for n in range(Ncpu)]
for p in procs:
    p.start()
for p in procs:
    p.join()


print('generating K4, K5, K6 data')
procs = [Process(target=gen_K456,
                 args=(nlist_chunks[n], len(nlist_chunks[n])))
         for n in range(Ncpu)]
for p in procs:
    p.start()
for p in procs:
    p.join()


exit()
# old single-core script no longer used
#########################################
print('generating K1, K2, K3 data for neutrinos')
for n in range(NlogYe):
    Ye = 10**logYe_arr[n]
    percent = 0
    for i in range(Nlogeta):
        eta = 10**logeta_arr[i]
        if floor(i*100./Nlogeta) > percent:
            print('%d percent' % percent)
            percent += 10
        for j in range(Nlogthe):
            the = 10**(logthe_arr[j])
            mu = eta*the + 1.
            # need baryonic mass fractions Xp, Xn, Xa
            Xa_guess = 0.01
            rho = 10**I2[i, j]*cst.mp/(Ye*Vunit)
            Xp, Xn, Xa = f_mass(Ye, rho, the, Xa_guess)
            # chemical potential of neutrinos
            eta_nu = eta - (Q-1.)/the - log(Xn/Xp)
            if eta_nu > eta_nu_arr.max() or eta_nu < eta_nu_arr.min():  # test
                print('eta_nu=%.3e outside the range!' % eta_nu)
                print('eta=%.3e, the=%.3e, Xp=%.3e, Xn=%.3e, Ye=%.3e' %
                      (eta, the, Xp, Xn, Ye))
                exit()
            Fnu_flux, Fnu_num, Fnu_flux_P = 0., 0., 0.  # 'Fnu_flux_P for Planck-mean'
            b = 0.30*Xp + 0.36*Xn + 0.21*Xa
            for ix in range(Nx):
                x = x_arr[ix]
                phip = phip_arr[ix]
                a = x**2/(exp(x-eta_nu) + 1)*phip
                temp = (x*the + Q - mu)/the
                if temp > log(1. + 1./epsilon_small):
                    f_block = 1.    # no blocking effect
                elif temp < log(epsilon_small):  # Taylor expansion
                    f_block = exp(temp)  # extremely strong blocking
                else:
                    f_block = 1 - 1./(exp(temp) + 1)    # in the middle
                c = 1.41*Xn*f_block*(x + Q/the)*sqrt((x+Q/the)**2 - 1/the**2)
                d = a/(b*x**2 + c)
                Fnu_flux += x*d
                Fnu_num += d
                Fnu_flux_P += x*a*(b*x**2 + c)
            Fnu_flux *= A*dt
            Fnu_num *= A*dt
            Fnu_flux_P *= A*dt
            K1[i, j] = log10(the**2*10**F3_intp(eta_nu)/Fnu_flux)
            K2[i, j] = log10(the**2*10**F2_intp(eta_nu)/Fnu_num)
            K3[i, j] = log10(the**2*Fnu_flux_P/10**F2_intp(eta_nu))
    print('writing K1, K2, K3 into files for Ye=%.3e' % Ye)
    Kdata = [K1, K2, K3]
    for m in [1, 2, 3]:
        logYe = log10(Ye)
        savename = 'K%d_logYe_m%dp%04d'\
                   % (m, floor(-logYe), round((-logYe) % 1*1e4))
        wrt_content = Kdata[m-1]
        f = open(fdir + savename + '.txt', 'w')
        f.seek(0)  # go to the beginning
        for i in range(Nlogeta):
            if i >= 1:
                f.write('\n')
            for j in range(Nlogthe):
                if j == 0:
                    f.write('%.8e' % wrt_content[i, j])
                else:
                    f.write('\t%.8e' % wrt_content[i, j])
        f.close()

print('generating K4, K5, K6 data for anti-neutrinos')
for n in range(NlogYe):
    Ye = 10**logYe_arr[n]
    percent = 0
    for i in range(Nlogeta):
        eta = 10**logeta_arr[i]
        if floor(i*100./Nlogeta) > percent:
            print('%d percent' % percent)
            percent += 10
        for j in range(Nlogthe):
            the = 10**(logthe_arr[j])
            mu = eta*the + 1.   # electron chemical potential
            # need baryonic mass fractions Xp, Xn, Xa
            Xa_guess = 0.01
            rho = 10**I2[i, j]*cst.mp/(Ye*Vunit)
            Xp, Xn, Xa = f_mass(Ye, rho, the, Xa_guess)
            # chemical potential of neutrinos
            eta_bnu = -(eta - (Q-1.)/the - log(Xn/Xp))
            b = 0.30*Xp + 0.36*Xn + 0.21*Xa     # scattering term
            x0 = (Q+1.)/the  # the break point
            A0 = x0/x_t_max
            # for the first segment of integral: 0 < x < x0 (no absorption)
            Fbnu_flux1, Fbnu_num1, Fbnu_flux1_P = 0., 0., 0.
            for ix in range(Nx):
                x = x_arr[ix]
                phip = phip_arr[ix]
                a = 1./(exp(x-eta_bnu) + 1)*phip
                Fbnu_flux1 += x*a/b
                Fbnu_num1 = a/b
                Fbnu_flux1_P += x**3 * a * b * x**2
            Fbnu_flux1 *= A0*dt
            Fbnu_num1 *= A0*dt
            Fbnu_flux1_P *= A0*dt
            # for the second segment of integral: x0 < x < xmax(~50)
            Fbnu_flux2, Fbnu_num2, Fbnu_flux2_P = 0., 0., 0.
            if x0 < xmax:   # this segment exists
                for ix in range(Nx):
                    x = x_arr[ix] + x0
                    phip = phip_arr[ix]
                    a = 1./(exp(x-eta_bnu) + 1)*phip
                    d = (x-Q/the)**2 - 1/the**2
                    if d < epsilon_small:
                        d = 2*x_arr[ix]  # linear expansion for small (x-x0)
                    y = 0.51099*the*x   # energy of bnu in MeV
                    gfit = y**(-0.07056 + 0.02018*log(y) - 0.001953*(log(y))**3)
                    temp = (x*the - Q + mu)/the
                    if temp > log(1. + 1./epsilon_small):
                        f_block = 1.    # no blocking effect
                    elif temp < log(epsilon_small):  # Taylor expansion
                        f_block = exp(temp)  # extremely strong blocking
                    else:
                        f_block = 1 - 1./(exp(temp) + 1)    # in the middle
                    c = 1.48*gfit*Xp*f_block*(x-Q/the)*sqrt(d)/x**2   # absorption
                    Fbnu_flux2 += x*a/(b + c)
                    Fbnu_num2 += a/(b + c)
                    Fbnu_flux2_P += x**3 * a * (b + c) * x**2
                Fbnu_flux2 *= A*dt
                Fbnu_num2 *= A*dt
                Fbnu_flux2_P *= A*dt
            Fbnu_flux = Fbnu_flux1 + Fbnu_flux2
            Fbnu_num = Fbnu_num1 + Fbnu_num2
            Fbnu_flux_P = Fbnu_flux1_P + Fbnu_flux2_P
            K4[i, j] = log10(the**2*10**F3_intp(eta_bnu)/Fbnu_flux)
            K5[i, j] = log10(the**2*10**F2_intp(eta_bnu)/Fbnu_num)
            K6[i, j] = log10(the**2*Fbnu_flux_P/10**F3_intp(eta_bnu))
    print('writing K4, K5, K6 into files for Ye=%.3e' % Ye)
    Kdata = [K4, K5, K6]
    for m in [4, 5, 6]:
        logYe = log10(Ye)
        savename = 'K%d_logYe_m%dp%04d'\
                   % (m, floor(-logYe), round((-logYe) % 1*1e4))
        wrt_content = Kdata[m-4]
        f = open(fdir + savename + '.txt', 'w')
        f.seek(0)  # go to the beginning
        for i in range(Nlogeta):
            if i >= 1:
                f.write('\n')
            for j in range(Nlogthe):
                if j == 0:
                    f.write('%.8e' % wrt_content[i, j])
                else:
                    f.write('\t%.8e' % wrt_content[i, j])
        f.close()

print('finish generating all data, ready for interpolation!')


# I1 = np.zeros((Nlogeta, Nlogthe))   # Ppair*V/mec2
# I2 = np.zeros((Nlogeta, Nlogthe))   # Ye*rho*V/mp
# I3 = np.zeros((Nlogeta, Nlogthe))   # spair*rho*V/k_B
# I4 = np.zeros((Nlogeta, Nlogthe))   # q_nu,thin
# I5 = np.zeros((Nlogeta, Nlogthe))   # q_bnu,thin
# I6 = np.zeros((Nlogeta, Nlogthe))   # lam_nu,thin
# I7 = np.zeros((Nlogeta, Nlogthe))   # lam_bnu,thin
# I8 = np.zeros((Nlogeta, Nlogthe))   # U_pair_th*V/mec2  (not including rest mass)

# K1 = np.zeros((Nlogeta, Nlogthe))  # kap_nu for each Ye (Rosseland-mean)
# K2 = np.zeros((Nlogeta, Nlogthe))  # tkap_nu
# K3 = np.zeros((Nlogeta, Nlogthe))  # kap_P_nu   (Planck-mean)
# K4 = np.zeros((Nlogeta, Nlogthe))  # kap_bnu
# K5 = np.zeros((Nlogeta, Nlogthe))  # tkap_bnu
# K6 = np.zeros((Nlogeta, Nlogthe))  # kap_P_bnu
from plt_funcs import *
from collections import OrderedDict
from math import floor, ceil
from read_intp_data import *
from read_intp_funcs import read_data_K123456

plt_name = 'I'  # 'I', 'K'
n = 6   # 1--7 for 'I', 1--6 for 'K'
Ye = 0.3    # pick a Ye for 'K'

# ---- explanations
# I1: Ppair*V/mec2
# I2: Ye*rho*V/mp
# I3: spair*rho*V/k_B
# I4: q_nu,thin/Xp * mp/(Gam*mec2)
# I5: q_bnu,thin/Xn * mp/(Gam*mec2)
# I6: lam_nu,thin /Gam
# I7: lam_bnu,thin /Gam

# K1: kap_nu * mp/sigma0
# K2: kap_t_nu * mp/sigma0
# K3: kap_nu_P * mp/sigma0
# K4: kap_bnu * mp/sigma0
# K5: kap_t_bnu * mp/sigma0
# K6: kap_bnu_P * mp/sigma0

unit_conver, label = 1, ''
if plt_name == 'I':
    plt_arr = Idata[n-1][:]
    if n == 1:  # pressure
        unit_conver = cst.mec2/Vunit
        label = r'$\log(P_{e^-} + P_{e^+}) [\rm erg\, cm^{-3}]$'
        savename = 'pressure'
    elif n == 2:  # density
        unit_conver = 1./Vunit
        label = r'$\log(n_{e^-} - n_{e^+}) [\rm cm^{-3}]$'
        savename = 'density'
    elif n == 3:  # entropy
        unit_conver = cst.kB/Vunit
        label = r'$\log(s_{e^-} + s_{e^+}) [\rm entropy\, per\, volume]$'
        savename = 'entropy'
    elif n == 4:
        unit_conver = Gam*cst.mec2/cst.mp
        label = r'$\log(q_\nu/X_p) [\rm erg\, g^{-1}\, s^{-1}]$'
        savename = 'q_nu'
    elif n == 5:
        unit_conver = Gam*cst.mec2/cst.mp
        label = r'$\log(q_{\bar{\nu}}/X_n) [\rm erg\,g^{-1}\,s^{-1}]$'
        savename = 'q_bnu'
    elif n == 6:
        unit_conver = Gam
        label = r'$\log(\lambda_{e^-p}) [\rm s^{-1}]$'
        savename = 'lam_ep'
    else:   # n == 7
        unit_conver = Gam
        label = r'$\log(\lambda_{e^+n}) [\rm s^{-1}]$'
        savename = 'lam_en'
else:   # for plt_name = 'K'
    read_data_K123456(Ye, logYe_arr, logeta_arr, logthe_arr, Kdata, fdir)
    plt_arr = Kdata[n-1][:]
    if n == 1:  # kap_nu
        unit_conver = sigma0/cst.mp
        label = r'$\log(\kappa_{\rm R, \nu}) [\rm cm^{2}\,g^{-1}]$'
        savename = 'kap_nu'
    elif n == 2:  # kap_t_nu
        unit_conver = sigma0/cst.mp
        label = r'$\log(\tilde{\kappa}_{\rm R, \nu}) [\rm cm^{2}\,g^{-1}]$'
        savename = 'kap_t_nu'
    elif n == 3:  # kap_nu_P
        unit_conver = sigma0/cst.mp
        label = r'$\log(\kappa_{\rm P, \bar{\nu}}) [\rm cm^{2}\,g^{-1}]$'
        savename = 'kap_nu_P'
    elif n == 4:  # kap_bnu
        unit_conver = sigma0/cst.mp
        label = r'$\log(\kappa_{\rm R, \bar{\nu}}) [\rm cm^{2}\,g^{-1}]$'
        savename = 'kap_bnu'
    elif n == 5:  # kap_t_bnu
        unit_conver = sigma0/cst.mp
        label = r'$\log(\widetilde{\kappa}_{\rm R, \bar{\nu}}) [\rm cm^{2}\,g^{-1}]$'
        savename = 'kap_t_bnu'
    else:   # n == 6, kap_bnu_P
        unit_conver = sigma0/cst.mp
        label = r'$\log(\kappa_{\rm P, \bar{\nu}}) [\rm cm^{2}\,g^{-1}]$'
        savename = 'kap_bnu_P'


plt_content = plt_arr + np.log10(unit_conver)
xarr = logeta_arr
yarr = logthe_arr

fig, ax = pl.subplots(1, 1, figsize=(12, 9), sharex='none')

xlabel = r'$\log \eta$'
ylabel = r'$\log \Theta$'
zlabel = label
cmap = 'coolwarm'

min_value, max_value = plt_content.min(), plt_content.max()
mid_point = 0.5*(min_value + max_value)
print('max, min =', max_value, min_value)
# levels = [ceil(min_value + 0.05*(max_value-min_value)),
#           round(min_value + 0.3*(max_value-min_value)),
#           round(mid_point),
#           round(max_value - 0.3*(max_value-min_value)),
#           floor(max_value - 0.05*(max_value-min_value))]
# levels = list(OrderedDict.fromkeys(levels))
levels = np.arange(ceil(min_value), floor(max_value)+1, 2)
CB_levels = levels
CB_ticklabels = ['%.1f' %l for l in levels]

flag_contour = True
pltimg(ax, xarr, yarr, plt_content, xlabel, ylabel, zlabel, cmap,
       CB_levels, CB_ticklabels, flag_contour)

fig.savefig(figdir + savename + '.png')

exit()

# --- 2D plot (old script)
fig = pl.figure(figsize=(13, 13))
ax = fig.add_axes([0.16, 0.27, 0.73, 0.7])
ax.set_xlabel(r'log($\eta$)')
ax.set_ylabel(r'log($\Theta$)')
im = pl.imshow(plt_content.transpose(),
               interpolation='nearest', origin='lower',
               cmap='bwr', aspect='auto',
               extent=(min(logeta_arr), max(logeta_arr),
                       min(logthe_arr), max(logthe_arr)))
im.set_clim(vmin=plt_content.min(), vmax=plt_content.max())
# --- color bar
cbar_ax = fig.add_axes([0.16, 0.11, 0.73, 0.04])
min_value, max_value = plt_content.min(), plt_content.max()
mid_point = 0.5*(min_value + max_value)
print('max, min =', max_value, min_value)
levels = [ceil(min_value + 0.05*(max_value-min_value)),
          round(min_value + 0.3*(max_value-min_value)),
          round(mid_point),
          round(max_value - 0.3*(max_value-min_value)),
          floor(max_value - 0.05*(max_value-min_value))]
levels = list(OrderedDict.fromkeys(levels))
CB = fig.colorbar(im, cax=cbar_ax, ticks=levels,
                  orientation='horizontal')
pl.xlabel(label)
# --- contours
strs = [str(num) for num in levels]
X, Y = np.meshgrid(logthe_arr, logeta_arr)
CS = ax.contour(Y, X, plt_content,
                levels, colors='k', linewidths=3, alpha=0.5)
fmt = {}
for l, s in zip(CS.levels, strs):
    fmt[l] = s
pl.clabel(CS, CS.levels[::], inline=True, fmt=fmt,
          fontsize=20, colors='k')
fig.savefig(fdir + 'atest/' + savename + '.png')

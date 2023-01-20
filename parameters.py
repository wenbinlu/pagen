import cst
from math import pi

# directory to store/read data
fdir = './data_tables/'
figdir = './figs/'


# other parameters
Q = 2.531   # neutron-proton mass ratio in units of mec^2
Gam = 6.93e-4  # beta-interaction rate (related to free neutron decay time)
sigma0 = 1.76e-44  # typical weak interaction cross-section (cm^2)
chi = 5.8*(cst.me/cst.mp)  # correction to sigma_abs_bnup due to nucleon recoil
Vunit = cst.h**3/(8*pi*(cst.me*cst.c)**3)   # volume unit

epsilon_small = 1.3162e-15   # a very small number

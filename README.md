# pnega
Weak reaction rates and neutrino opacities in a plasma with protons (p), neutrons (n), electron-position pairs (e), gamma-rays (g), and alpha particles (a)

The underlying assumptions are: (1) all 'pnega' species are in LTE at the same temperature T; (2) neutrinoes are non-degenerate and the plasma is optically thin to neutrinos; (3) p, n, a are described by (non-degenerate, non-relativistic) Maxwell-Boltzmann distributions; (4) electron-positrons are described by Fermi-Dirac distributions (allowing for arbitrary degeneracy and relativistic temperatures); (5) gamma-rays are described by the Planck distribution; (6) the only nuclides in the system are p, n, a (other nuclides are photo-disintegrated); (7) neutrino opacities only include scattering and absorption by p, n, a (ignoring inelastic collisions with pairs).

All quantities only depends on three independent thermodynamic variables: Ye (electron fraction), Theta (dimensionless temperature), eta = (mu-1)/Theta (mu is the electron chemical potential). Some quantities only depend on Theta and eta.

The output options include: pair pressure, total density, pair entropy, cooling rates due to neutrino/anti-neutrino emission, electron/positron capture rates, Rosseland-mean/Planck-mean opacities.

Procedure:

(1) Change 'fdir' and 'figdir' to be the directories where you would like to store numerical tables and figures, respectively

(2) Set the resolutions in 'data_init' by changing Nlogeta, Nlogthe, NlogYe --- these are the number of grid points along the three dimensions. You can also change lower and upper bounds of the grid. The resolution in eta_nu (energy spectrum of neutrinoes) is by default set to be Neta_nu = 10000, and it is suggested to keep it this way since it is not computational expensive.

(3) Run 'Python gen_F_data.py' to create tables named 'FN.txt' for N=2,3,4,5,6 (this only takes a few seconds)

(4) Run 'Python gen_IK_data.py' to create tables named 'IN.txt' for N=1,2,3,4,5,6,7,8 (this takes about 20 seconds) and 'KN_logYe_*.txt' for N=1,2,3,4,5,6 (this takes about 20 minutes using 16 cpus)

Now all the relevant numerical tables are stored. You can in principle calculate whatever you want about the plasma at arbitrary density, temperature and Ye. You can visualize the numerical tables using 'plt_numeric_data.py'.

For instance, to plot electron-capture rate = neutrino emission rate per proton as a function of eta and Theta, you choose

plt_name = 'I'
n = 6
Ye = 0.3 (in fact, lam_nu does not depend on Ye)

By running 'Python plt_numeric_data.py', you will obtain a figure named 'lam_ep.png' in the directory specified by 'figdir'.

To plot the positron capture rate = anti-neutrino emission rate per neutron as a function of eta and Theta, you choose

plt_name = 'I'
n = 7
Ye = 0.3 (in fact, lam_bnu does not depend on Ye)

Then you will obtain a figure named 'lam_en.png'.

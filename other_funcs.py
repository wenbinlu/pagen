import numpy as np
import cst
from math import log, exp, sqrt
from parameters import *


def phi_p(t):   # without the normalization
    return exp(-1./t - 1./(1-t))


# mass fractions for given (Ye, rho, the) -- bisection method
def f_mass(Ye, rho, the, Xa_guess):
    rho10 = rho/1e10    # note: rho is in units of g/cc
    if the < .5:  # to avoid extremely negative exponent
        Xa = 2*min(Ye, 1.-Ye)
        Xp = max(Ye-0.5*Xa, epsilon_small)
        Xn = max(1-Ye-0.5*Xa, epsilon_small)
        return Xp, Xn, Xa
    RHS = 38.6*rho10**-1.5*the**2.25*exp(-27.65/the)
    # print('in f_mass function, rho=%.3e, the=%.3e' % (rho, the))  # for test
    # print('RHS=%.4e' % RHS)   # for test
    if RHS < epsilon_small:  # RHS is extremely close to zero
        Xa = 2*min(Ye, 1.-Ye)
        Xp = max(Ye-0.5*Xa, epsilon_small)
        Xn = max(1-Ye-0.5*Xa, epsilon_small)
        return Xp, Xn, Xa
    else:
        Xa_max = 2*min(Ye, 1.-Ye)   # maximum possible Xa
        if Xa_guess > Xa_max:
            print('incorrect Xa_guess=%.4e causing Xp < 0 or Xn < 0 at Ye=%.4e'
                  % (Xa_guess, Ye))
            exit()
        Xa = Xa_guess   # use the input guess value
        # use bisection method to find the root for: f=0
        f = (Ye-0.5*Xa)*(1-Ye-0.5*Xa)*Xa**-0.5 - RHS   # monotonic decreasing
        # print('f(Xa_guess) = %.5e' % f)   # for test
        if f > 0:
            Xa_low = Xa_guess
            Xa_up = Xa_max
        else:
            Xa_low = 0.
            Xa_up = Xa_guess
        Xa = 0.5*(Xa_low + Xa_up)
        tol = 1e-6
        # while abs(f) > tol:   # this criterion doesn't work well
        while (Xa_up-Xa_low)/Xa > tol and Xa_up-Xa_low > epsilon_small:
            f = (Ye-0.5*Xa)*(1-Ye-0.5*Xa)*Xa**-0.5 - RHS
            if f > 0:
                Xa_low = Xa
            else:
                Xa_up = Xa
            Xa = 0.5*(Xa_low + Xa_up)
            # print('Xa_low=%.5e, Xa = %.5e, Xa_up=%.5e, f(X_a)=%.5e'
            #       % (Xa_low, Xa, Xa_up, f))   # for test
            # print('c1=%.5e, c2=%.5e' % ((Xa_up-Xa_low)/Xa, Xa_up-Xa_low)) # test
        # the correct Xa has been found
        Xp = max(Ye-0.5*Xa, epsilon_small)
        Xn = max(1-Ye-0.5*Xa, epsilon_small)
        return Xp, Xn, Xa


# baryonic contribution to specific entropy, for given composition and temperature
def sb_entropy(Xp, Xn, Xa, rho, the):
    np = rho*Xp/cst.mp
    nn = rho*Xn/cst.mp
    na = 0.25*rho*Xa/cst.mp
    temp = sqrt(pi/2.)*(the*cst.mp/cst.me)**1.5/Vunit
    s_temp = 2.5*(Xp + Xn + 0.25*Xa)\
        + Xp*log(temp/np) + Xn*log(temp/nn)\
        + 0.25*Xa*log(4**1.5*temp/na)
    return s_temp*cst.kB/cst.mp


def fermi_block(x):  # Fermi blocking factor (to avoid numerical error)
    if x > log(1. + 1./epsilon_small):
        f_block = 1.    # no blocking effect
    elif x < log(epsilon_small):  # Taylor expansion
        f_block = exp(x)  # extremely strong blocking
    else:
        f_block = 1 - 1./(exp(x) + 1)    # in the middle
    return f_block


def sigmoid(x):  # smoothed version of a step function
    if x > 40:
        return 1.
    if x < -40:
        return 0.
    return (1. + exp(-x))**-1

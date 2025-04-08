# The expressions used here are from the following links:
#https://scholarlypublications.universiteitleiden.nl/access/item%3A3485847/view
#https://ar5iv.labs.arxiv.org/html/0803.4507
from setup import *
import scipy.integrate as integral

import SPIworkflow.SPIutils as spi
import numpy as np

def get_gaunt(T,nu):
    """
    Computes the Gaunt factor in the radio regime, using the expressions in references
    19 and 21 in Cox (2000). 

    OUTPUT - g:  float - Gaunt factor
    INPUT  - T:  float - Temperature of the wind,  in Kelvin
           - nu: float - Observing frequency, in Hz
    It also takes into account the ionization state (Z - read from __init__.py)
    """
    g = 10.6 + 1.90 * np.log10(T) - 1.26 * np.log10(Z * nu)
    return g

def ff_absorption(M_star, nu_ecm, T_wind, m_av, X_p, M_dot, R_ff_in, R_ff_out, NSTEPS_FF):
    '''
    Computes the optical depth, taunu, of free-free absorption at a given frequency, nu_ecm
    It uses the expression for the effective linear absorption coefficient, as given in
    Cox (2000), with a slight modification of the nomenclature. 

    OUTPUT - tau_nu   - Optical depth of free-free absoption
             kappa_nu - effective linear absorption coefficient 
             alpha_nu - kappa_nu / (ne * np)
    INPUT: 

    M_star - Star mass, in grams
    nu_ecm - Electron cyclotron frequency, in Hz
    T_wind - Temperature of the wind, in K (equal to T_corona for an isoth wind)
    m_av   - average particle mass of the stellar wind, in gr    
    X_p    - Fraction of protons
    M_dot   - stellar mass-loss rate, in units of the Sun mass-loss rate
    R_ff_in- Distance (from stellar center) where SPI emission takes place, in cm (>= 1.0)
    R_ff_out- Distance from stellar center where free-free absorption becomes negligible, in cm
    NSTEPS_FF - Number of steps for the computation of free-free absorption 

    Other parameters: 
    Z - ionisation state
    h - Planck constant
    k_B - Boltzmann constant
    '''
    #dist_absorption = np.linspace(R_ff_in, R_ff_out, NSTEPS_FF)
    # Logarithmic spacing is better if R_ff_out >> R_ff_in
    dist_absorption = np.logspace(np.log10(R_ff_in), np.log10(R_ff_out), NSTEPS_FF)

    v_sound, r_sonic, v_sw = spi.v_stellar_wind(dist_absorption, M_star, T_wind, m_av)

    # Computing densities n_p and n_e to estimate kappa_nu
    n_sw = spi.n_wind(M_dot, dist_absorption, v_sw, m_av) 
    n_p  = n_sw * X_p
    n_e  = n_sw * (1 - X_p)

    # Gaunt factor
    g = get_gaunt(T_wind, nu_ecm)

    alpha_nu = 3.692 * 10**8 * (1 - np.exp((-h * nu_ecm)/(k_B * T_wind))) * Z**2 * g * T_wind**(-1/2) * nu_ecm**(-3)
    
    # Effective linear absorption coefficient
    kappa_nu = alpha_nu * n_e * n_p

    # Optical depth
    tau_nu = integral.trapezoid(kappa_nu, dist_absorption)

    return tau_nu, kappa_nu, alpha_nu

from setup import *
import numpy as np

starname='GJ486'
d      = 8.0827 * pc               # Distance to stellar system , in  cm
R_star = 0.339 * R_sun    # Stellar radius in cm
M_star = 0.333 * M_sun      # Stellar mass in g,
P_rot_star = 49.9 * day  # Rotation period  of star, in sec
B_star =    240         # Stellar surface magnetic field
    
Exoplanet='GJ486 B'
Mp = 3.00*M_earth # Planetary mass, in grams
#Rp = spi.Rp_Zeng(1.15)
#Rp *= R_earth # Planetary radius, in cm
Rp = 1.343 *R_earth # R prop to M^0.28 for Terran worlds https://arxiv.org/pdf/2311.12593
P_orb = 1.4671205 # orbital period of planet, in days
r_orb  = 0.01713 * au   # orbital distance, in cm

#Additional parameters (type np.nan if unknown)
M_star_dot = np.nan # Stellar mass loss rate in units of the mass loss rate of the Sun (NOT the mass of the Sun). The mass loss rate of the Sun is 2e-14.  
T_corona = np.nan

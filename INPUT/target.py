from setup import *
import numpy as np

#Values of both objects are extracted from Fitzmaurice et al. 2024
starname='test'
d      =  130 * pc       # Distance to stellar system , in  cm (extracted from Table 1, page 3)
R_star =  1.5 * R_sun    # Stellar radius in cm (extracted from Table 1, page 3)
M_star =  0.5 * M_sun    # Stellar mass in g, (extracted from Table 1, page 3)
P_rot_star = 4   * day   # Rotation period  of star, in sec  (extracted from Table 1, page 3)
B_star =  2000               # Stellar surface magnetic field (in page 12 they say they use B=16-214G for this type of slow rotating 
                            #M-drawf, but as we detected the system at L-band, 1-2GHz, we assume at least B_star=360G bc 2.8x360~1GHz)
    
Exoplanet='test b' #in this case, we have a brown dwarf, not an exoplanet
M_Jupiter= 317.8*M_earth
R_Jupiter=11.209*R_earth
Mp = 1 *M_Jupiter # Planetary mass, in grams (extracted from Table 2, page 4)
Rp =1*R_Jupiter # Planetary radius, in cm (value for 2Gyr in Table 3 page 13 bc in page 11 they say they estimated 2.5Gyr for the BD)
P_orb =  30.004294 # orbital period of planet, in days  (extracted from Table 2, page 4)
r_orb  = 0.15 * au   # orbital distance, in cm  (extracted from Table 2, page 4 ,a0 esta en unidades de mas, hay que pasarlo a au)

#Additional parameters (type np.nan if unknown)
M_star_dot =  5  # Stellar mass loss rate in units of the mass loss rate of the Sun (NOT the mass of the Sun). For the Sun is 2e-14.  
T_corona = 1.5e6 #in K (extracted from Eq. 24, page 12) 
B_pl = 8.4 #interacting planet/object magnetic field, in Gauss

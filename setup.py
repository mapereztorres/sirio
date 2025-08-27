import numpy as np
from SPIworkflow.constants import *

#python spirou.py 2>&1 | tee logfile.txt

#####################################################
# ARRAY OF PLANETS TO COMPUTE RADIO EMISSION FROM SPI
#####################################################
# Choose whether to carry out the computations for a 
# table containing multiple targets (INPUT_TABLE = True), or
# for a single target (INPUT_TABLE = False). 

# Uncomment the line that applies
#INPUT_TABLE = True
INPUT_TABLE = False

#######################################################################
#  STUDY CASES
#  STUDY == 'D_ORB' - Predicted flux as a function of orbital distance
#  STUDY == 'M_DOT' - Predicted flux as a function of star mass-loss rate
#  STUDY == 'B_PL'  - Predicted flux as a function of planetary magnetic field
#  
#######################################################################

# Uncomment the line that applies
STUDY = "D_ORB"
#STUDY = "M_DOT"
#STUDY = "B_PL"

# STUDY = "D_ORB" SETUP
D_ORB_LIM = np.nan
#D_ORB_LIM = 3000

#  STUDY = 'M_DOT' SETUP
#
# M_DOT_MIN, M_DOT_MAX: Minimum and maximum mass-loss rates to carry out the study of M_dot
# In units of M_dot_sun 
# M_DOT_STRETCH: Number of points per dex in the STUDY of M_DOT
M_DOT_STRETCH = 50
M_DOT_MIN = 1e-1
M_DOT_MAX = 1e+2

#  STUDY = 'B_PL' SETUP
#
# B_PL_MIN, B_PL_MAX: Minimum and maximum planetary magnetic field to carry out 
# the study of "B_PL". Values in Gauss
#
STEP = 0.05
B_PL_MIN = 0
B_PL_MAX = 4 

##
# Distance (from the centre of the star) where SPI emission takes place (in units of R_star)
R_SPI = 1.0

####################################################
# Stellar wind FREE ABSORPTION of SPI radio emission
####################################################

### Consider free-free absorption (True => Yes; False => No)
#freefree = True
freefree = False

# Ionization state (Z = 1 - fully ionized hydrogen)
Z = 1 

# Distance, measured from R_SPI, in units of R_star, where free-free absorption becomes negligible. 
R_ff_OBSERVER = 10000

### NSTEPS_FF: number of points for the distance array
#NSTEPS_FF = 1000000
NSTEPS_FF = 10000
#####################################################
### SETTING UP VALUES TO PREDICT SPI RADIO EMISSION
#####################################################

#####################################
# STELLAR MAGNETIC FIELD SETUP
#####################################
# Setting the stellar magnetic field geometry and the value of the 
# intensity of the planetary magnetic field
# 
# Stellar magnetic field geometry
Bfield_geom_arr=['open_parker_spiral','closed_dipole','pfss']
#Bfield_geom_arr=['pfss']

# INDEX of the DIPOLE
# Q_DIPOLE = 3.0 => DIPOLE
Q_DIPOLE = 3.0

# DIPOLE TILT - Tilt of the dipolar magnetic moment of the star wrt the rotation axis
# (in radians).
# POLAR_ANGLE - Angle measured between the orbital plane of the planet and the (star's)
# polar axis, measured from the polar axis. In radians.
# NOTE! Use always a positive value, never zero. For safety reasons, the code uses a
# TOLERANCE parameter to prevent hravoc.
# Set POLAR_ANGLE to np.pi/2 for a planet in the equatorial plane of the star.
# 
DIPOLE_TILT = 0.0
POLAR_ANGLE = np.pi/2
AZIMUTH     = 0.0

TOLERANCE = 1e-2

# R_T and DELTA_R are used for the hybrid model, i.e. the transition from a pure dipole
# to an open Parker spiral.
# DELTA_R - Width of the transition region, 
# 
#R_T = 10.0 # * np.sin(POLAR_ANGLE)**2 
DELTA_R = 0.05 * R_sun

# Potential source surface radius (PFSS), in units of R_star
# R__SS is defined as below, due to the boundary conditions imposed on the components of
# the closed magnetic field geometry (see, e.g., Eqns. 5 and 6 in Jardine+2002, MNRAS)
R_SS = 4.5
#R_SS = 5.5
#R_SS = 8
#R_SS = 22
#R_SS = 10
#R_SS = 21
#R_SS = 34
#R_SS = 18
# R_ALFVEN_GUESS - Initial guess for the Alfvén radius (see get_R_alfven in SPIutils.py)
#                  in units of R_star (stellar radii)
R_ALFVEN_GUESS = 20.0

#####################################
# PLANET MAGNETIC FIELD SETUP
#####################################

# magnetized_pl_arr is a [False,True] array
# False: Unmagnetized planet 
# True : Magnetized planet
# 
#magnetized_pl_arr = [False, True]
magnetized_pl_arr = [True]

# Default planetary magnetic field, in Tesla
# bfield_earth is defined in constants.py
B_PLANET_DEFAULT = bfield_earth

# Setting the stellar magnetic field geometry and the value of the 
# Computation of planetary magnetic field 
# B_pl_law = 'Sano' => Uses Sano's scaling law (Sano 1993)
# B_pl_law = 'None' => Doesn't use any scaling law. Uses B_PLANET_DEFAULT instead.
B_planet_law = 'Sano'
#B_planet_law = 'None'

# K_MAGNETOPAUSE - factor by which the magnetopause currents enhance
# the magnetospheric magnetic field at the magnetopause, which is a value
# between 2 and 3. ADD REFERENCE
K_MAGNETOPAUSE = 2.0

###################################
# OBSERVED PARAMETERS
###################################

# Observing frequency, in  Hz
#freq_obs = 400e6
# Assumed rms noise figure, in mJy
RMS = 0.015
# Representative bandwidth of the ECMI emission, in Hz
#Delta_nu_obs = freq_obs/2 

# Observed flux density at some given frequency, in mJy  (This should be modified in the
# future)
#flux = 1.0 

###################################
# STELLAR PARAMETERS
###################################

# Temperature of the corona,  
# base density at the corona 
# magnetic field at the surface of the star
# "isothermality" of stellar plasma 

# Default value: 2.0e6 K. Solar value
#T_CORONA_DEF = 2.0e6 
T_CORONA_DEF=10**6.67



#####################################
# STELLAR WIND
#####################################

# Is the stellar wind plasma assumed to be isothermal?
ISOTHERMAL = True 

# Fully ionized, purely hydrogen stellar wind plasma (=> 50% protons, 50% electrons)
X_p = 0.5 # fraction of protons

#####################################
# SUB-ALFVENIC INTERACTION parameters
#####################################
 
# ALPHA_SPI - relative strength of the sub-Alfvénic interaction. 
# We assume that the planet has a highly conductive atmosphere, i.e., ALPHA_SPI = 1
ALPHA_SPI = 1

# Efficiency factor in converting Poynting flux into ECM radio emission.
#
# We use the naming convention in Zarka (2018) and Zarka (2024), in the 
# Handbook of Exoplanets, i.e., "beta" (his Eq. 1 in Zarka 2024), 
# which we rename as "beta_eff". 
# Zarka (2024) shows that beta_eff is in the range from 1e-4 (BETA_EFF_MIN) to 1e-2
# (BETA_EFF_MAX)
# 
#BETA_EFF_MIN = 1e-4; BETA_EFF_MAX = 1e-2 
BETA_EFF_MIN = 1e-3; BETA_EFF_MAX = 1e-3 

# Fraction of Poynting flux going to dissipated power (Eq. 2 in Zarka 2024)
# Zarka (2024) quotes a value of 0.2 +/- 0.1
EPSILON = 0.2 

# theta_M - Angle between planetary field and stellar field (rad) in the planet rest frame
THETA_M = np.pi/2

# Minimum and maximum kinetic energy of electrons emitting via ECM, in keV
# Electrons are mildly relativistic, so 20 to 200 keV energies should be good values, 
# but a choice from 10 to 511 keV (=m_e * c^2) should be also  OK. 
EKIN_MIN = 1 
EKIN_MAX = 20

## Beam solid angle of the ECM radio emission
## 
## Choose whether to compute the beam solid angle (COMPUTE_BSA = True), 
## or use fixed values (COMPUTE_BSA = False)
COMPUTE_BSA = False

## The standard avlue for OMEGA_MIN for SPI is
## taken from the Io-Jupiter interaction.
## OMEGA_JUPITER_IO = 0.16 sterradians (DAM emission from a single flux tube)
OMEGA_MIN = OMEGA_JUPITER_IO

# OMEGA_MAX due to star-planet interaction is at most a few times OMEGA_MIN
# Note: Many papers have wrongly assumed OMEGA_MAX = 1.6 (sterr), which is the 
# beam solid angle of the whole Jupiter auroral oval. 
OMEGA_MAX = 3*OMEGA_JUPITER_IO

#####################################
# PLOTTING AND WRITING SETUP
#####################################
#Plot graph for M_A
PLOT_M_A = False
## LINEWIDTH 
LW = 3
# PLOTOUT = True => plot graphs to external files
# PLOTOUT = False => plot in terminal
PLOTOUT = True 
# DRAW_RMS?
DRAW_RMS = False
# DRAW Little Earth?
DRAW_EARTH = True
LIMS_MA = True
LIM_MA_LOW = 1e-2
LIM_MA_HIGH = 1e0
FLUX_LIMS = True
FLUX_LOW  = 3*RMS * 1e-1
FLUX_HIGH = 3*RMS * 1e2
YLIMLOW   = 1e-3
YLIMHIGH  = 1e2
#PLOT_R_eff = False

### WELCOME TO SIRIO
print('###########################################################')
print('###                                                     ###')
print('### SIRIO code @Pérez-Torres & Peña-Moñino 2025         ###')
print('###                                                     ###')
print('### If you use fully or partial parts of this code,     ###')
print('### please acknowledge our work by citing the paper     ###') 
print('### XXXX and indicating the url address for the code:   ###') 
print('### https://github.com/mapereztorres/sirio              ###') 
print('###                                                     ###')
print('###########################################################\n')

###
print('You have chosen the following setup parameters for SIRIO:\n')

###
if INPUT_TABLE == False:
    print(f'RUNNING THE CODE FOR A SINGLE OBJECT: INPUT_TABLE = {INPUT_TABLE}\n')
else:
    print(f'RUNNING THE CODE FOR MULTIPLE OBJECTS: INPUT_TABLE = {INPUT_TABLE}\n')

###
if STUDY == "D_ORB":
    print('CARRYING OUT A STUDY OF RADIO EMISSION VS ORBITAL SEPARATION: STUDY == D_ORB\n')
elif STUDY == "M_DOT":
    print('CARRYING OUT A STUDY OF RADIO EMISSION vs STELLAR MASS LOSS RATE: STUDY == M_DOT\n ')
elif STUDY == "B_PL":
    print('CARRYING OUT A STUDY OF RADIO EMISSION VS PLANETARY MAGNETIC FIELD: STUDY == B_PL\n')

###
'''
#['open_parker_spiral','closed_dipole','closed_pfss']
if 'open_parker_spiral' in Bfield_geom_arr:
    print('RUNNING FOR AN OPEN PARKER SPIRAL MAGNETIC FIELD GEOMETRY\n')    
if 'closed_dipole' in Bfield_geom_arr:
    print('RUNNING FOR A CLOSED DIPOLAR MAGNETIC FIELD GEOMETRY\n')    
if 'pfss' in Bfield_geom_arr:   
    print('RUNNING FOR A CLOSED PFSS MAGNETIC FIELD GEOMETRY\n')   
'''   
ind=0
printing_str='RUNNING FOR THE FOLLOWING GEOMETRIES:'
while ind < len(Bfield_geom_arr):
   if ind == len(Bfield_geom_arr)-1 and ind>0:
       printing_str += ' AND'
   elif ind>0:
       printing_str += ','
   printing_str += ' ' + Bfield_geom_arr[ind].replace('_', ' ').upper()
   ind=ind+1
print(printing_str+'\n') 
    
###    
if len(magnetized_pl_arr)==2:    
    print('RUNNING FOR BOTH A MAGNETIZED AND A NON-MAGNETIZED PLANET\n')
else:
    if magnetized_pl_arr[0] == False:
        print('RUNNING FOR A NON-MAGNETIZED PLANET\n')    
    else:
        print('RUNNING FOR A MAGNETIZED PLANET\n')    
        
if (B_planet_law) == 'Sano':
    print('THE PLANETARY MAGNETIC FIELD IS ESTIMATED FOLLOWING SANO\'s SCALING LAW\n')
else:
    print('THE PLANETARY MAGNETIC FIELD IS ASSUMED TO BE THAT OF THE EARTH\n')


###
if freefree == True:
    print('FREE FREE ABSORPTION EFFECTS ARE CONSIDERED\n')
else:
    print('FREE FREE ABSORPTION EFFECTS ARE IGNORED\n')    
    
###
print(f'THE EFFICIENCY FACTOR IN CONVERTING POYNTING FLUX \nTO RADIO EMISSION RANGES FROM {BETA_EFF_MIN:.1e} UP TO {BETA_EFF_MAX:.1e}\n')

###
if COMPUTE_BSA == False: 
    print(f'THE BEAM SOLID ANGLE OF THE AURORAL EMISSION RANGES \nFROM {OMEGA_MIN:.2f} STERRADIANS UP TO {OMEGA_MAX:.2f} STERRADIANS\n')
else:
    print(f'THE BEAM SOLID ANGLE OF THE AURORAL EMISSION RANGES WILL BE COMPUTED BASED\nON THE VALUES OF THE MINIMUM AND MAXIMUM KINETIC ENERGY OF THE ELECTRONS\nEMITTING VIA THE ECM MECHANISM: FROM {EKIN_MIN:.0f} KEV UP TO {EKIN_MAX:.0f} KEV\n')

### 
print(f'You are assuming an observing rms value of {RMS*1e3:.1f} μJy\n')    

print('###########################################################\n')

'''
import os
import subprocess

# Save originals
_original_os_system = os.system
_original_subprocess_run = subprocess.run
_original_subprocess_call = subprocess.call
_original_subprocess_popen = subprocess.Popen

# Patch os.system
def debug_os_system(cmd):
	print(f"[DEBUG os.system] Command: {cmd!r}")
	return _original_os_system(cmd)

# Patch subprocess.run
def debug_subprocess_run(*args, **kwargs):
	print(f"[DEBUG subprocess.run] args: {args}, kwargs: {kwargs}")
	return _original_subprocess_run(*args, **kwargs)

# Patch subprocess.call
def debug_subprocess_call(*args, **kwargs):
	print(f"[DEBUG subprocess.call] args: {args}, kwargs: {kwargs}")
	return _original_subprocess_call(*args, **kwargs)

# Patch subprocess.Popen
def debug_subprocess_popen(*args, **kwargs):
	print(f"[DEBUG subprocess.Popen] args: {args}, kwargs: {kwargs}")
	return _original_subprocess_popen(*args, **kwargs)

# Apply patches
os.system = debug_os_system
subprocess.run = debug_subprocess_run
subprocess.call = debug_subprocess_call
subprocess.Popen = debug_subprocess_popen

'''
## Star-planet Interaction (model)
## Sub-Alfvenic flux for both a Parker spiral magnetic field configuration and a closed dipolar field configuration


# import statements
import os
import shutil
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from scipy.special import lambertw
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
#matplotlib.rc_file_defaults()
plt.style.use(['bmh','SPIworkflow/spi.mplstyle'])

from output import OutputWriter

## Getting parameters to predict SPI radio emission
#
# Observing parameters: Observin frequency, assumed rms noise, Delta_nu
# Stellar parameters: T_corona, n_base_corona, B_field, Isothermality of plasma
# Geometry of the sub-Alfvénic interaction: ALPHA_SPI, THETA_M
#
from setup import *

# Import useful constants and functions to be used in the code
from SPIworkflow.constants import *
import SPIworkflow.SPIutils as spi
import SPIworkflow.freefree as ff
from SPIworkflow.load_data import get_spi_data, create_data_tables, load_target, table_lists

import importlib

# Create output directory for the results 
# Return df_planets and df_no_planets
# Create CARMENES tables for targets 
# with planets only, and with no planets, unless those already exist
# 
# outdir, df_planets, df_no_noplanets = create_data_tables()

# Call empty lists to be used later in out_table
# all_lists = table_lists()
# LPM  - DOCUMENT BETTER
#dipole_mag_pl_lists   = table_lists()
#dipole_unmag_pl_lists = table_lists()
#spiral_mag_pl_lists   = table_lists()
#spiral_unmag_pl_lists = table_lists()

if not os.path.exists('./OUTPUT'):
    os.makedirs('./OUTPUT')
#################################################################
################## GENERAL EMISSION PARAMETERS  #################
#################################################################

# Compute min and max speed of electrons emitting via ECM, in units of the speed of light 
beta_min = spi.beta_keV(EKIN_MIN);  beta_max = spi.beta_keV(EKIN_MAX)
            
# Beam solid angle covered by the ECM emission, in sterradians
Omega_min, Omega_max = spi.beam_solid_angle(COMPUTE_BSA, beta_min, beta_max)

# Read in the input data to estimate radio emission from SPI
# INPUT_TABLE is defined in setup.py
if INPUT_TABLE == True:
    # Read in the input data to estimate radio emission from SPI
    data = get_spi_data(infile_data = './INPUT/table.csv')

    ############## CHECK THAT THE DATA TABLE IS CORRECT
    print('Reading table: ')
    pd.set_option('display.max_columns', None)
    print(data)

    ############## TABLE INITIALIZATION 
    # 
    # Create column for M_star_dot to fill it with values
    data['M_star_dot(M_sun_dot)']=''
    # If bfield_star(gauss) is missing, set it to np.nan
    #data['bfield_star(gauss)'].replace('', np.nan, inplace=True)
    data.replace({'bfield_star(gauss)': ''}, {'bfield_star(gauss)': np.nan}, inplace=True)

    # If p_rot is missing, set it to np.nan
    #data['p_rot(days)'].replace('', np.nan, inplace=True)
    # Remove targets without p_rot
    data.dropna(subset=['p_rot(days)'], inplace=True)

    # Do not use stars with P_rot smaller than 10 days
    #data = data[data['p_rot(days)'] > 10.0]
    #data['radius_planet(r_earth)'].replace('', np.nan, inplace=True)
    data.replace({'radius_planet(r_earth)': ''}, {'radius_planet(r_earth)': np.nan}, inplace=True)
    data.reset_index(inplace = True) # to prevent funny jumps in the indices

    ############## PRINT INDEX OF EACH PLANET AFTER RESETTING INDICES IN data
    print('SIRIO WILL BE RUN FOR THE FOLLOWING TARGETS\n')
    print(data['planet_name'])
    print('\n')
    planet_array = range(len(data))
else:
    # else we deal with one single target
    planet_array = [0] 

    file_name = 'INPUT.target' 
    
    print(f'READING INFO FROM {file_name.replace(".","/")}'+'.py\n')

    # import all parameters into a single variable
    imported_parameters = importlib.import_module(file_name) 
        
    # convert imported_parameters into global variables
    globals().update({k: v for k, v in imported_parameters.__dict__.items() if not k.startswith('__')})

#all_targets = []
#flux_r_S_list = []
#flux_reconnect_list = []
#flux_sb_list = []
    
    
flux_data = {
    'open_parker_spiral': {
        'x': [], 'Flux_r_S': [], 'Flux_reconnect': [], 'Flux_sb': [], 'obs_freq': [], 'name': []
    },
    'closed_dipole': {
        'x': [], 'Flux_r_S': [], 'Flux_reconnect': [], 'Flux_sb': [], 'obs_freq': [], 'name': []
    },
    'pfss': {
        'x': [], 'Flux_r_S': [], 'Flux_reconnect': [], 'Flux_sb': [], 'obs_freq': [], 'name': []
    }
}


    
for indi in planet_array:
# Read parameters from a table containing multiple targets or from a single target file
    if INPUT_TABLE == True: 
        starname,d, R_star, M_star, P_rot_star, B_star, Exoplanet, Mp, Rp, r_orb, P_orb,eccentricity, q, Q = load_target(data, indi)
        T_corona = data['T_corona(MK)'][indi]*1e6
        M_star_dot = data['mdot(mdot_sun)'][indi]
        print(f'RUNNING SIRIO FOR EXOPLANET {planet_array[indi]}: {Exoplanet}\n')
    else: 
        print(f'RUNNING SIRIO FOR EXOPLANET {Exoplanet}\n')
    # Fill B_star column if empty. Uses original units from table
    # B_star - stellar magnetic field at the equator, in Gauss
    if pd.isna(B_star):
         B_star= spi.B_starmass(star_mass=data['mass_star(m_sun)'][indi],Prot=data['p_rot(days)'][indi])
         data['bfield_star(gauss)'][indi]=B_star
         # Eventually include uncertainties in B _star
         data['e_bfield_star'][indi]='TBD'
    Exoplanet=Exoplanet.replace("'s","")
    starname=starname.replace("'s","")
    
    B_spi = B_star * (R_SPI)**-3                        # Magnetic field where the SPI emission takes place (R_SPI)  

    #nu_ecm = 2.8e6 * B_star # cyclotron freq, in Hz

    # cyclotron (= ECM) frequency, and ECM bandwith, in Hz
    nu_ecm = e * B_spi/(2*np.pi * m_e * c) 
    Delta_nu_cycl = nu_ecm 

    #  Check whether M_star_dot is read from input table/file
    if np.isnan(M_star_dot):
        if INPUT_TABLE == True:       
            #data['M_star_dot(M_sun_dot)'][indi] = spi.Mdot_star(R_star=data['radius_star(r_sun)'][indi], M_star=data['mass_star(m_sun)'][indi], Prot_star=data['p_rot(days)'][indi])/M_sun_dot
            #M_star_dot = data['M_star_dot(M_sun_dot)'][indi]     # Stellar mass loss rate in solar units 
            M_star_dot = spi.Mdot_star(R_star=R_star/R_sun, M_star=M_star/M_sun, Prot_star=P_rot_star/day)/M_sun_dot
        else:
            M_star_dot = spi.Mdot_star(R_star=R_star/R_sun, M_star=M_star/M_sun, Prot_star=P_rot_star/day)/M_sun_dot
        print(f'Estimated value of M_star_dot: {M_star_dot:.2e} M_dot_sun')
    #  Check whether M_star_dot is read from input table/file
    if np.isnan(T_corona):
        T_corona = T_CORONA_DEF
     


    ###############################################

    # Electron gyrofrequency and ECM bandwidth 
    #gyrofreq = e*B_spi/(2*np.pi * m_e * c) # in cgs units
    #Delta_nu_cycl = gyrofreq # Hz - width of ECMI emission  assumed to be  (0.5 * gyrofreq), 
    
    # Max. orbital distance, in units of R_star
    d_orb_max = max(2*r_orb/R_star, D_ORB_LIM) 

    # The type of STUDY (D_ORB, M_DOT or B_PL) is set up in setup.py 
    # and tells us whether the computaion is done 
    # as a function of the planetary orbital distance: STUDY = "D_ORB"
    # as a function of the stellar mass loss rate    : STUDY = "M_DOT"
    # as a function of the planetary magnetic field  : STUDY = "B_PL"
    # 
    # Nsteps defines the size of the array
    #
    if STUDY == "D_ORB":
        Nsteps = int(2*d_orb_max)
        d_orb  = np.linspace(1.02, d_orb_max, Nsteps) * R_star # Array of (orbital) distances to the star, in cm 
        #d_orb  = np.linspace(1.02, r_orb/R_star, 2)
        M_star_dot_arr = np.array([M_star_dot]) # Convert to a numpy array of 1 element for safety reasons
    elif STUDY == "M_DOT":
        d_orb  = np.array([r_orb])
        Nsteps = int(M_DOT_STRETCH * np.log10(M_DOT_MAX/M_DOT_MIN))
        M_star_dot_arr = np.logspace(np.log10(M_DOT_MIN), np.log10(M_DOT_MAX), Nsteps)
    elif STUDY == "B_PL":
        Nsteps = round( (B_PL_MAX - B_PL_MIN) / STEP)
        d_orb = np.array([r_orb])
        M_star_dot_arr = np.array([M_star_dot]) # Convert to a numpy array of 1 element for safety reasons

    # get array of orbital and corotation speeds (v_orb and v_corot) and Omega_star (float)
    v_orb, v_corot, Omega_star = spi.get_velocity_comps(M_star, d_orb, P_rot_star) 

    # Angular speed of the planet, in s^(-1). Note that it's an array
    Omega_planet =  v_orb / d_orb 

    # Get wind composition, from the fraction of protons
    X_e, mu, m_av = spi.wind_composition(X_p)

    # Compute stellar wind velocity at each value of d_orb
    v_sound, r_sonic, v_sw = spi.v_stellar_wind(d_orb, M_star, T_corona, m_av)

    # Stellar wind velocity at the closest distance to the star, in cm/s
    v_sw_base = v_sw[0]    
     
    # Terminal speed of the isothermal stellar wind, in cm/s
    v_sw_terminal = spi.get_v_sw_terminal(R_star, M_star, T_corona, m_av)
    
    # Magnetic confinement parameter at the stellar equator
    eta_star = spi.get_eta_star(B_star, R_star, M_star_dot_arr, v_sw_terminal)

    #Alfven radius at the specified POLAR_ANGLE, in stellar radii
    if STUDY == "D_ORB" or STUDY == "B_PL":
        R_alfven = spi.get_R_alfven(eta_star, colatitude=POLAR_ANGLE)
        R_alfven =np.array(R_alfven)
        R_alfven_pole = spi.get_R_alfven(eta_star, colatitude = 0)
        R_alfven_pole = np.array(R_alfven_pole)
    if STUDY == "M_DOT":
        R_alfven = spi.get_R_alfven_alt(eta_star, colatitude=POLAR_ANGLE)
        R_alfven =np.array(R_alfven)
        R_alfven_pole = spi.get_R_alfven_alt(eta_star, colatitude = 0)
        R_alfven_pole = np.array(R_alfven_pole)

    #Alfven radius at the pole, in stellar radii
    #R_alfven_pole = spi.get_R_alfven(eta_star, colatitude = 0)

    # Latitude above which all magnetic field lines are open
    theta_A_deg, latitude = spi.get_theta_A(R_alfven_pole)

    # Plasma number density at base of the corona, in cm^(-3)
    n_base_corona = spi.n_wind(M_star_dot_arr, R_star, v_sw_base, m_av) 

    # Maximum plasma frequency at the base of the corona, in Hz. If the ECM
    # freq is less than the plasma frequency, the emission is completely absorbed 
    nu_plasma_corona = spi.plasma_freq(n_base_corona * X_e) 

    # Relative speed between stellar wind and obstacle
    # Angle between radial vector and relative velocity
    # From Eq. 23 of Turnpenney+18 - Second term of RHS 
    # The vector v_rel = v_sw - v_orb (Eq. 5 in Saur+13, and see also Fig. 1 in Turnpenney+18)
    v_rel = np.sqrt(v_orb**2 + v_sw**2)  # in cm/s
    #angle_v_rel = np.arctan(v_orb/v_sw)  # in radians
    angle_v_rel = np.arctan2(v_orb,v_sw)  # in radians
    # Compute n_sw_planet, the number density of the wind at the orbital distance to the planet. 
    # If the stellar plasma is assumed to be isothermal, then 
    # the density falls down as ~ R^(-2) * v_sw^(-1).
    # Alternatively, we fix the density at the distance of the planet from the host star.
    if ISOTHERMAL:
        #n_sw_planet = n_sw_base / (d_orb/R_star)**2 / (v_sw/v_sw_base) # Plasma density at distance (R/R_star)
        n_sw_planet = spi.n_wind(M_star_dot_arr, d_orb, v_sw, m_av) # Plasma number density at distance (R/R_star)
        
    else:
        # WARNING: This (arbitrary) value of 1e4 for n_sw_planet to be set up in setup.py
        #n_sw_planet = np.ones(len(d_orb)) * 1e4  
        n_sw_planet = np.ones(Nsteps) * 1e4  

    rho_sw_planet = m_av * n_sw_planet #wind density at the distance to the planet, in g * cm^(-3)
    
    for ind in range(len(Bfield_geom_arr)):
        for ind1 in range(len(magnetized_pl_arr)):
            # Bfield_geom_arr defines the geometry of the magnetic field (closed / open)
            #if Bfield_geom_arr[ind]:
            #    selected_geometry="OPEN PARKER MAGNETIC FIELD GEOMETRY"
            #else:
            #    selected_geometry="CLOSED DIPOLAR MAGNETIC FIELD GEOMETRY"
            #['open_parker_spiral','closed_dipole','pfss']
            if Bfield_geom_arr[ind] == 'open_parker_spiral':
                selected_geometry="OPEN PARKER MAGNETIC FIELD GEOMETRY"
            elif Bfield_geom_arr[ind] == 'closed_dipole':
                selected_geometry="CLOSED DIPOLAR MAGNETIC FIELD GEOMETRY"
            elif Bfield_geom_arr[ind] == 'pfss':   
                selected_geometry="HYBRID PFSS - PARKER SPIRAL MAGNETIC FIELD GEOMETRY"
            # get magnetic field components
            B_r, B_phi, B_sw, angle_B, theta, geom_f = spi.get_bfield_comps(Bfield_geom_arr[ind], B_spi, d_orb, R_star, v_corot, v_sw, angle_v_rel)
            
            # Compute Alfvén parameters in the stellar wind at a distance d_orb 
            v_alf, M_A, v_alf_r, M_A_radial = spi.get_alfven(rho_sw_planet, B_sw, B_r, v_rel, v_sw)
            # defines whether planet is unmagnetized (magnetized_pl_arr[ind1] = 0), or magnetized (magnetized_pl_arr[ind1] = 1)
            if magnetized_pl_arr[ind1]: # magnetized planet
                planet_magnetized='MAGNETIZED PLANET'
                if B_planet_law == 'Sano':
                    # Planetary magnetic field, using Sano's (1993) scaling law, in units of B_earth # Assumes a tidally locked planet, i.e., the rotation period of the
                    # planet equals its orbital one. 
                    # WARNING: For small rotation periods, the inferred magnetic field
                    # is too large to be reliable at all.
                    r_core, rho_core, magn_moment_planet, B_planet_arr = spi.bfield_sano(M_planet = Mp / M_earth, 
                                               R_planet = Rp / R_earth, 
                                               Omega_rot_planet = Omega_planet / Omega_earth)  
                    B_planet_arr *= bfield_earth  # B_planet_arr, in Tesla
                else: 
                    B_planet_arr = np.ones(len(Omega_planet)) * B_PLANET_DEFAULT  # B_planet_arr, in Tesla
                
                B_planet_arr    *=  Tesla2Gauss #  B_planet_arr, in Gauss 

            else:  # unmagnetized planet
                planet_magnetized='UNMAGNETIZED PLANET'
                B_planet_arr  = np.zeros(len(d_orb)) # unmagnetized planet
            print(f'\nRUNNING FOR {selected_geometry} \nAND {planet_magnetized}')
            
            if STUDY == "B_PL":
                B_planet_Sano = B_planet_arr # Planet magnetic field at r_orb. 1-element array, in Gauss. 
                B_planet_arr  = np.linspace(B_PL_MIN, B_PL_MAX, Nsteps)
                

                    
            # Effective radius of the obstacle, in cm
            # Case 1. À la Saur+2013. (NOT currently used)
            #R_obs_Saur = spi.get_Rmp_Saur(Rp, THETA_M, B_planet_arr, B_sw)


            # Case 2. À la Zarka (2007), Turnpenney+2018, etc.
            #
            # Compute radius of magnetopause, Rmp as balance of wind and planet's
            # pressures
            
            # Planet pressure, in erg/cm3 - only the magnetic component is considered
            P_B_planet  = spi.get_P_B_planet(B_planet_arr) 

            # Stellar wind pressure, in erg/cm3
            P_sw, P_dyn_sw, P_th_sw, P_B_sw = spi.get_P_sw(n_sw_planet, v_rel, T_corona, B_sw, mu)

            #if Bfield_geom_arr[ind] == 'closed_dipole':
            #    r_ss_estimate = spi.get_rss(P_th_sw, P_dyn_sw,P_B_sw, d_orb)

            eta = spi.get_confinement(P_dyn_sw, P_B_sw)
            #alfven_alt = spi.get_alfven_alt(eta, POLAR_ANGLE)
            
            
            #magnetic field of the planet:
            v_orb_pl = (G * M_star / r_orb)**0.5
            Omega_pl =  v_orb_pl / r_orb 
            _,_,_,Bplanet_field  = spi.bfield_sano(M_planet = Mp / M_earth, R_planet = Rp / R_earth, Omega_rot_planet = Omega_pl / Omega_earth)
            Bplanet_field  *=  bfield_earth
            Bplanet_field  *=  Tesla2Gauss
            
            Reff_lanza=spi.get_Reff_lanza(B_sw,B_pl=Bplanet_field)
            
            
            # Radius of magnetopause, in cm
            Rmp_saur= spi.get_Rmp(P_B_planet, P_dyn_sw, P_th_sw, P_B_sw)
            Rmp=Rmp_saur* Rp
            #Rmp=Reff_lanza* Rp*np.ones(len(Rmp_saur))
    
            # The effective radius (in cm) is the radius of the magnetopause
            # If R_pl_eff < R_p, force R_pl_eff = R_p (R_obs cannot be smaller than Rp
            R_obs = np.copy(Rmp)
            R_obs[R_obs < Rp] = Rp 
            R_obs_normalized = R_obs/Rp 

            # Find value of Bp where R_obs where is larger than Rp
            #indices_R_obs_larger_Rp = np.argwhere(R_obs > Rp)
            #if indices_R_obs_larger_Rp.size > 0:
                #Bp value where R_magnetosphere is larger than Rp
            #    B_planet_eff_rad = B_planet_arr[indices_R_obs_larger_Rp[0]]          

            
            ## Get Poynting fluxes for the Alfvén wing model (Zarka 2007, Saur 2013,
            # Turnpenney+2018)
            
            # Get Poynting flux using Eq. 55 in Saur+2013 (S_poynt) and Eq. 8 in Zarka
            # 2007 (S_poyn_Z), in cgs units. They coincide, except for a factor 2. 
            # In mks units
            S_poynt, S_poynt_Z = spi.get_S_poynt(R_obs, B_sw, v_alf, v_rel, M_A, ALPHA_SPI, geom_f)

            # Get fluxes at Earth, in cgs units for both Saur+ (Flux_r_S...) and
            # Zarka/Lanza (Flux_r_S_Z...),  in erg/s/Hz/cm2
            Flux_r_S_min, Flux_r_S_max =  spi.get_Flux(Omega_min, Omega_max, Delta_nu_cycl, d, S_poynt)
            Flux_r_S_Z_min, Flux_r_S_Z_max =  spi.get_Flux(Omega_min, Omega_max, Delta_nu_cycl, d, S_poynt_Z)
            # Compute flux density for an intermediate value of eps (in log scale)
            Flux_r_S_inter = 10**((np.log10(Flux_r_S_max) + np.log10(Flux_r_S_min))/2)
            
            
            # Get flux for the reconnection model (Lanza 2009, A&A)
            R_obs_reconnect=np.copy(R_obs)
            R_obs_reconnect[np.isclose(R_obs_reconnect, Rp, atol=1e-2)] = np.nan
            
            #S_reconnect, P_d, P_d_mks = spi.get_S_reconnect(geom_f,R_obs_reconnect, B_sw, v_rel, gamma = 0.5)
            S_reconnect, P_d, P_d_mks = spi.get_S_reconnect(geom_f,R_obs_reconnect, B_sw, v_rel, gamma = 1)
            
            
            
            #20250124-TBC: if R_obs < Rp:
            Flux_reconnect_min, Flux_reconnect_max = spi.get_Flux(Omega_min, Omega_max, Delta_nu_cycl, d, S_reconnect)
            Flux_reconnect_inter = 10**((np.log10(Flux_reconnect_max) + np.log10(Flux_reconnect_min))/2)




            # Get flux for the stretch and break model (Lanza 2013 and Strugarek et al.2022)            
            S_sb = spi.get_S_stretch_and_break(R_obs=R_obs_reconnect, B_sw=B_sw, v_rel=v_rel, B_planet_arr=B_planet_arr,geom_f=geom_f)
            
            
            Flux_sb_min, Flux_sb_max = spi.get_Flux(Omega_min, Omega_max, Delta_nu_cycl, d, S_sb)
            Flux_sb_inter = 10**((np.log10(Flux_sb_max) + np.log10(Flux_sb_min))/2)           
            
            
            ###
            ### COMPUTATION OF FREE-FREE Absorption by the stellar wind 
            ###
            alphamatrix=[]
            
            #altitude over stellar surface where SPI takes place, in cm
            R_ff_in  = R_star * R_SPI 

            # Limit for integration of free-free absorption, in cm
            R_ff_out = R_star * R_ff_OBSERVER 
            pdn=pd.DataFrame(columns=np.linspace(R_ff_in, R_ff_out, NSTEPS_FF))

            # Unabsorbed flux density
            Flux_r_S_min_no_abs=np.copy(Flux_r_S_min)
            Flux_r_S_max_no_abs=np.copy(Flux_r_S_max)
            Flux_r_S_inter_no_abs=np.copy(Flux_r_S_inter)
            
            Flux_reconnect_min_no_abs = np.copy(Flux_reconnect_min)
            Flux_reconnect_max_no_abs = np.copy(Flux_reconnect_max)
            Flux_reconnect_inter_no_abs = np.copy(Flux_reconnect_inter)
            
            
            Flux_sb_min_no_abs = np.copy(Flux_sb_min)
            Flux_sb_max_no_abs = np.copy(Flux_sb_max)
            Flux_sb_inter_no_abs = np.copy(Flux_sb_inter)           
            # Compute flux density, taking into account free-free absorption
            if freefree == True:
                print('Applying ff-absorption')
                # Generate empty list to be filled with absorption values
                absorption = []
                for elem in M_star_dot_arr:
                    # mdot - one-element array from each value of M_star_dot_arr
                    M_dot = np.array(elem) 

                    #Compute optical depth for f-f absorption (tau_nu)
                    tau_nu, kappa_nu, alpha_nu = ff.ff_absorption(M_star, nu_ecm, T_corona, m_av, X_p, M_dot, R_ff_in, R_ff_out, NSTEPS_FF)
                    absorption.append(np.exp(-tau_nu))

                # Convert list into numpy array
                absorption_factor = np.array(absorption)

                Flux_r_S_min    *= absorption_factor
                Flux_r_S_max    *= absorption_factor
                Flux_r_S_Z_min *= absorption_factor
                Flux_r_S_Z_max *= absorption_factor    
                Flux_r_S_inter  *= absorption_factor
                Flux_reconnect_min *= absorption_factor
                Flux_reconnect_max *= absorption_factor
                Flux_reconnect_inter *= absorption_factor
                Flux_sb_min *= absorption_factor
                Flux_sb_max *= absorption_factor   
                Flux_sb_inter *= absorption_factor         
                
                

            # Find out the position of the planet in the distance array
            d_diff = np.abs((d_orb-r_orb)/R_star)
            loc_pl = np.where(d_diff == d_diff.min())
            
            #Bplanet_field = B_planet_arr[loc_pl][0]

            M_star_dot_diff = np.abs(M_star_dot_arr - M_star_dot)
            M_star_dot_loc  = np.where(M_star_dot_diff == M_star_dot_diff.min())
            
            if STUDY == 'D_ORB':
                xnom = r_orb/R_star
                closest_index = np.argmin(np.abs(d_orb/R_star - xnom))
                
            elif STUDY == 'M_DOT':
                xnom = M_star_dot
                closest_index = np.argmin(np.abs(M_star_dot_arr - M_star_dot))
                
            elif STUDY == 'B_PL':
                xnom = Bplanet_field
                closest_index = np.argmin(np.abs(B_planet_arr - Bplanet_field))
            

            

            v_alf_planet_closest=v_alf[closest_index]
            
            
            Sigma_P, Sigma_A, alpha_interaction_strength=spi.get_interaction_strength(r_orb,B_star,Bplanet_field,v_alf_planet_closest,M_A,geom_f)

            loc_max_alpha = np.argmax(alpha_interaction_strength)
            loc_min_alpha = np.argmin(alpha_interaction_strength)            
            
            Flux_r_S_inter_planet = Flux_r_S_inter[closest_index]
            Flux_reconnect_inter_planet = Flux_reconnect_inter[closest_index]
            Flux_sb_inter_planet = Flux_sb_inter[closest_index]
            
            #all_targets.append(Exoplanet)
            #flux_r_S_list.append(Flux_r_S_inter_planet)
            #flux_reconnect_list.append(Flux_reconnect_inter_planet)
            #flux_sb_list.append(Flux_sb_inter_planet)
            
            # Plotting is different, depending on the "STUDY" case
            if STUDY == 'D_ORB':
                x = d_orb / R_star # (distance array, in units of R_star)
            elif STUDY == 'M_DOT':
                x = M_star_dot_arr # (M_star_dot_arr array, in units of M_dot_sun)
            elif STUDY == 'B_PL':
                x = B_planet_arr # (B_planet_arr array, in Gauss )
            plt.figure(figsize=(10,6))
            
            geom = Bfield_geom_arr[ind] 
            flux_data[geom]['x'].append(x[closest_index])
            flux_data[geom]['Flux_r_S'].append(Flux_r_S_inter_planet)
            flux_data[geom]['Flux_reconnect'].append(Flux_reconnect_inter_planet)
            flux_data[geom]['Flux_sb'].append(Flux_sb_inter_planet)
            flux_data[geom]['name'].append(Exoplanet)
            flux_data[geom]['obs_freq'].append(2.8*B_star)

            '''      
            #M_star,d_orb_planet,P_rot_star,M_star,T_corona,m_av,Bfield_geom_arr[ind], B_spi, R_star
            d_orb_planet=r_orb           
            v_orb_planet, v_corot_planet, Omega_star_planet = spi.get_velocity_comps(M_star, d_orb_planet, P_rot_star) 
            v_sound_planet, r_sonic_planet, v_sw_planet = spi.v_stellar_wind(d_orb_planet, M_star, T_corona, m_av)
            v_rel_planet = np.sqrt(v_orb_planet**2 + v_sw_planet**2)  # in cm/s
            angle_v_rel_planet = np.arctan2(v_orb_planet,v_sw_planet)  # in radians
            B_r_planet, B_phi_planet, B_sw_planet_planet, angle_B_planet, theta_planet, geom_f = spi.get_bfield_comps(Bfield_geom_arr[ind], B_spi, d_orb_planet, R_star, v_corot_planet, v_sw_planet, angle_v_rel_planet)
            v_alf_planet, M_A_planet, v_alf_r_planet, M_A_radial_planet = spi.get_alfven(rho_sw_planet, B_sw, B_r, v_rel, v_sw)
            #flux at the planet position:
            S_poynt_planet, S_poynt_Z_planet = spi.get_S_poynt(R_obs, B_sw_planet, v_alf_planet, v_rel_planet, M_A_planet, ALPHA_SPI, geom_f)
            Flux_r_S_min_planet, Flux_r_S_max_planet =  spi.get_Flux(Omega_min, Omega_max, Delta_nu_cycl, d_planet, S_poynt_planet)
            Flux_r_S_Z_min_planet, Flux_r_S_Z_max_planet =  spi.get_Flux(Omega_min, Omega_max, Delta_nu_cycl, d_planet, S_poynt_Z_planet)
            Flux_r_S_inter_planet = 10**((np.log10(Flux_r_S_max_planet) + np.log10(Flux_r_S_min_planet))/2)
            
            
            R_obs_reconnect=np.copy(R_obs)
            R_obs_reconnect[np.isclose(R_obs_reconnect, Rp, atol=1e-2)] = np.nan
            S_reconnect_planet, P_d_planet, P_d_mks_planet = spi.get_S_reconnect(geom_f,R_obs_reconnect, B_sw_planet, v_rel_planet, gamma = 1)
            Flux_reconnect_min_planet, Flux_reconnect_max_planet = spi.get_Flux(Omega_min, Omega_max, Delta_nu_cycl, d_planet, S_reconnect_planet)
            Flux_reconnect_inter_planet = 10**((np.log10(Flux_reconnect_max_planet) + np.log10(Flux_reconnect_min_planet))/2)   
            
            
            S_sb_planet = spi.get_S_stretch_and_break(R_obs=R_obs_reconnect, B_sw=B_sw_planet, v_rel=v_rel_planet, B_planet_arr=B_planet_arr_planet,geom_f=geom_f)         
            Flux_sb_min_planet, Flux_sb_max_planet = spi.get_Flux(Omega_min, Omega_max, Delta_nu_cycl, d_planet, S_sb_planet)
            Flux_sb_inter_planet = 10**((np.log10(Flux_sb_max_planet) + np.log10(Flux_sb_min_planet))/2)   
            '''

            ###########################################################################
            ####                  PLOTTING                                         ####
            ###########################################################################
            # Create OUTPUT folder if it doesn't exist
            FOLDER = 'OUTPUT/' + str(Exoplanet.replace(" ", "_"))
            if not(os.path.isdir(FOLDER)):
                os.system('mkdir OUTPUT/' + str(Exoplanet.replace(" ", "_")))
            if not(os.path.isdir(FOLDER+'/FLUX_PDF')):
                os.system('mkdir OUTPUT/' + str(Exoplanet.replace(" ", "_"))+'/FLUX_PDF')
            if not(os.path.isdir(FOLDER+'/FLUX_PNG')):
                os.system('mkdir OUTPUT/' + str(Exoplanet.replace(" ", "_"))+'/FLUX_PNG')
            if not(os.path.isdir(FOLDER+'/CSV')):
                os.system('mkdir OUTPUT/' + str(Exoplanet.replace(" ", "_"))+'/CSV')
            if not(os.path.isdir(FOLDER+'/TXT')):
                os.system('mkdir OUTPUT/' + str(Exoplanet.replace(" ", "_"))+'/TXT')
            if not(os.path.isdir(FOLDER+'/DIAG_PDF')):
                os.system('mkdir OUTPUT/' + str(Exoplanet.replace(" ", "_"))+'/DIAG_PDF')
            if not(os.path.isdir(FOLDER+'/R_EFF_PDF')):
                os.system('mkdir OUTPUT/' + str(Exoplanet.replace(" ", "_"))+'/R_EFF_PDF')      
            if not(os.path.isdir(FOLDER+'/COMPARISON_PDF')):
                os.system('mkdir OUTPUT/' + str(Exoplanet.replace(" ", "_"))+'/COMPARISON_PDF')    
            if not(os.path.isdir(FOLDER+'/S_poynting')):
                os.system('mkdir OUTPUT/' + str(Exoplanet.replace(" ", "_"))+'/S_poynting')    
            # generate output files for the different magnetic field geometries
            geometry = "-" + Bfield_geom_arr[ind].replace('_','-') + '-Bstar'

            ### Plot received flux density as a function of distance from the star
            
            filename = 'plotting/plot_flux_density.py'
                #filename = 'plotting/plot_flux_and_radius.py'            
            with open(filename) as file:
                exec(file.read())            
            filename = 'plotting/plot_poynting.py'
            with open(filename) as file:
                exec(file.read())   
            
            if freefree == True and STUDY == 'M_DOT': #################
                #plt.figure(figsize=(8,11))
                plt.figure(figsize=(8,8))
                ax = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
                ax.plot(M_star_dot_arr, absorption_factor, color='k')
                ax.set_xscale('log')
                ax.set_xlabel(r"Mass Loss rate [$\dot{M}_\odot$]",fontsize=20)
                ax.set_ylabel(r"Fraction of transmitted flux")
                #ax.text(1e-1¡, 0, r'T$_{c} = $'+"{:.1f}".format(T_corona/1e6)+' MK', fontsize = 22)
                ax.text(1e1, 0.9, r'T$_{c} = $'+"{:.1f}".format(T_corona/1e6)+' MK', fontsize = 22)
                ax.set_facecolor("white")	
                secax = ax.secondary_yaxis('right', functions=(spi.identity,spi.identity))
                #Specific for YZ Cet from Pineda2023
                    
                plt.savefig(FOLDER + '/' + str(Exoplanet.replace(" ", "_"))
                        +'-'+'absorption_vs_mdot'+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'-'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'.pdf', bbox_inches='tight')
                plt.close()
                
            #### Plot effective radius variation
            
            
            
            filename = 'plotting/plot_effective_radius.py'
            with open(filename) as file:
                exec(file.read())
    
            
           
            li = [x,y_min.tolist(), y_max.tolist()]   
            
              
            #df = pd.DataFrame(data=li)
            #df = df.assign(column_name=column_series)
            #df.index = [STUDY, 'flux_min'+str(T_corona/1e6)+'MK', 'flux_max'+str(T_corona/1e6)+'MK']
            
            if freefree == True: 
                df = pd.DataFrame(zip(x,y_min, y_max), columns=[STUDY, 'flux_min'+str(T_corona/1e6)+'MK', 'flux_max'+str(T_corona/1e6)+'MK'])
                df.to_csv(os.path.join(FOLDER + '/CSV/' +outfile+ ".csv"))            
                df2= pd.DataFrame(zip(x,absorption_factor),columns=[STUDY,'abs_factor_'+str(T_corona/1e6)+'MK'])
                #df2.to_csv(FOLDER + '/' + str(Exoplanet.replace(" ", "_"))+'-'+'absorption_vs_mdot'+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'-'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'.csv')   
                df2.to_csv(FOLDER + '/CSV/' + str(Exoplanet.replace(" ", "_"))
                         +'-'+'absorption_vs_mdot'+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'-'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'.csv')   
                         
            
            ################################
            # DIAGNOSTIC PLOTS
            ################################
            
            #PLOTTING MAGNETIC FIELD LINES
            if STUDY=='D_ORB':
                filename = 'plotting/plot_magnetic_fields.py'
                with open(filename) as file:
                    exec(file.read()) 
                    
            
            filename = 'plotting/plot_diagnostic_plots.py'
            #Specific plots to benchmark against Turnpenney 2
                    
            with open(filename) as file:
                exec(file.read())
            
            #Plot for just the variation of M_A
            filename = 'plotting/plot_MA.py'    
            with open(filename) as file:
                exec(file.read())    
    
            ###########################################################
            ################### Send OUTPUT to external text file/s
            ###########################################################

            outfileTXT = os.path.join(FOLDER + '/TXT/' +outfile+'.txt')
            out_to_file = OutputWriter(outfileTXT)

            out_to_file.write_parameters(T_corona, M_star_dot, mu, d, R_star, M_star, P_rot_star, B_star, Exoplanet, Rp, Mp, r_orb, P_orb, loc_pl, M_star_dot_loc, n_base_corona,
                nu_plasma_corona, nu_ecm, Flux_r_S_min, Flux_r_S_max, rho_sw_planet, n_sw_planet, v_sw_base, Flux_r_S_Z_min,
                Flux_r_S_Z_max, v_sw, v_rel, v_alf, M_A, B_sw,Bplanet_field, Rmp, R_obs,x_larger_rms,x_smaller_rms,STUDY,Omega_min, Omega_max,R_obs_normalized,x_superalfv,
                Flux_r_S_inter_planet,Flux_reconnect_inter_planet,Flux_sb_inter_planet)

            print(f'\nSAVING USEFUL VALUES TO {outfileTXT}')
            

    filename = 'plotting/plot_model_comparison.py'
    with open(filename) as file:
        exec(file.read())          

    print(f'\nDONE WITH PLANET {Exoplanet}!!\n')
    print('###########################################################\n')






'''

colors = {
    'open_parker_spiral': 'orange',
    'closed_dipole': 'blue',
    'pfss': 'green'
}

for geom, data in flux_data.items():
    plt.scatter(data['x'], data['Flux_r_S'], color=colors[geom], label=f"{geom} - Flux_r_S")
    plt.scatter(data['x'], data['Flux_reconnect'], marker='x', color=colors[geom], label=f"{geom} - Flux_reconnect")
    plt.scatter(data['x'], data['Flux_sb'], marker='^', color=colors[geom], label=f"{geom} - Flux_sb")

plt.xlabel('x')
plt.ylabel('Flux')
plt.legend()
plt.tight_layout()
plt.show()
'''





plt.close('all') 

geom_list = ['open_parker_spiral', 'closed_dipole', 'pfss']
colors = {
    'Flux_r_S': 'orange',
    'Flux_reconnect': 'blue',
    'Flux_sb': 'green'
}




for geom in geom_list:
    #plot_data = flux_data[geom]
    #plot_data['Flux_r_S'] = plot_data['Flux_r_S'].apply(lambda x: x[0] if isinstance(x, list) else x)
    #plot_data['Flux_r_S'] = pd.to_numeric(plot_data['Flux_r_S'], errors="coerce")
    
    
    plot_data = pd.DataFrame(flux_data[geom])   # ensure DataFrame
    plot_data['Flux_r_S'] = plot_data['Flux_r_S'].apply(lambda x: x[0] if isinstance(x, list) else x)
    plot_data['Flux_r_S'] = pd.to_numeric(plot_data['Flux_r_S'], errors="coerce")
    #plot_data = plot_data[(plot_data['Flux_r_S'] > 1e-3) & (plot_data['Flux_r_S'] < 1e2)]
    
    
    
    
    
    
    
    
    #plot_data = plot_data[(plot_data['Flux_r_S'] > 1e-3) & (plot_data['Flux_r_S'] < 1e2)]
        
        
        
        
        
        
        
        
    #### Plotting flux vs orbital separation for all targets
    
        
    fig, ax = plt.subplots(figsize=(8, 12))

    #ax.set_ylim([1.001e-3, 1e2])
    #ax.set_xlim([10, 50])
    ax.set_yscale('log')
    
    ax.set_ylim([1.001e-6, 1e-1])
    ax.set_xlim([1, 150])

    # Scatter plots
    ax.scatter(plot_data['x'], plot_data['Flux_r_S'], color=colors['Flux_r_S'], s=60)
    #ax.scatter(plot_data['x'], plot_data['Flux_reconnect'], color=colors['Flux_reconnect'], marker='x', s=60)
    #ax.scatter(plot_data['x'], plot_data['Flux_sb'], color=colors['Flux_sb'], marker='^', s=60)

    # Assign numbers to names
    number_map = {name: i+1 for i, name in enumerate(plot_data['name'])}

    # Annotate points with numbers
    for xi, frs, name in zip(plot_data['x'], plot_data['Flux_r_S'], plot_data['name']):
        if pd.notna(xi) and pd.notna(frs) and np.isfinite(xi) and np.isfinite(frs):
            if 1e-6 < frs < 1e-1:
                if 1 < xi < 150:
                    ax.text(xi, frs, str(number_map[name]), fontsize=6, ha='right', va='bottom', rotation=0)

    # Create a legend mapping numbers to names
    mapping_handles = [Line2D([0], [0], color='none', label=f"{num} = {name}") 
                       for name, num in number_map.items()]

    ax.legend(handles=mapping_handles, loc='upper right', fontsize=6, title='Targets')

    # Axis labels & title
    ax.set_xlabel('Orbital separation / Stellar radius')
    ax.set_ylabel('Flux (mJy)')
    ax.set_title(geom.replace('_', ' ').title())
    #ax.set_xlabel(xlabel,fontsize=20)
    
    plt.tight_layout()
    plt.savefig('OUTPUT/alltargets_flux_separation_' + geom + '.pdf', bbox_inches='tight')
    plt.close()   
    
    
    
    
    
    
    
    
    #### Plotting flux vs observing frequency for all targets    
    
    fig, ax = plt.subplots(figsize=(8, 12))

    #ax.set_ylim([1.001e-3, 1e2])
    #ax.set_xlim([10, 50])
    ax.set_yscale('log')
    
    ax.set_ylim([1.001e-6, 1e-1])
    ax.set_xlim([1, 4000])

    # Scatter plots
    ax.scatter(plot_data['obs_freq'], plot_data['Flux_r_S'], color=colors['Flux_r_S'], s=60)
    #ax.scatter(plot_data['x'], plot_data['Flux_reconnect'], color=colors['Flux_reconnect'], marker='x', s=60)
    #ax.scatter(plot_data['x'], plot_data['Flux_sb'], color=colors['Flux_sb'], marker='^', s=60)

    # Assign numbers to names
    number_map = {name: i+1 for i, name in enumerate(plot_data['name'])}

    # Annotate points with numbers
    for xi, frs, name in zip(plot_data['obs_freq'], plot_data['Flux_r_S'], plot_data['name']):
        if pd.notna(xi) and pd.notna(frs) and np.isfinite(xi) and np.isfinite(frs):
            if 1e-6 < frs < 1e-1:
                if 1 < xi < 4000:
                    ax.text(xi, frs, str(number_map[name]), fontsize=6, ha='right', va='bottom', rotation=0)

    # Create a legend mapping numbers to names
    mapping_handles = [Line2D([0], [0], color='none', label=f"{num} = {name}") 
                       for name, num in number_map.items()]

    ax.legend(handles=mapping_handles, loc='upper right', fontsize=6, title='Targets')

    # Axis labels & title
    ax.set_xlabel('Observing frequency (MHz)')
    ax.set_ylabel('Flux (mJy)')
    ax.set_title(geom.replace('_', ' ').title())
    #ax.set_xlabel(xlabel,fontsize=20)
    
    plt.tight_layout()
    plt.savefig('OUTPUT/alltargets_flux_freq_' + geom + '.pdf', bbox_inches='tight')
    plt.close()   
        
    
    
    
    
    
    
    
    
    
    
    
    
print('###########################################################')
print(f'SIRIO HAS FINISHED SUCCESSFULLY!!\n')
print('###########################################################')

    
'''
for geom in geom_list:
    plot_data = flux_data[geom]

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.set_ylim([1.001e-3, 1e2])
    ax.set_xlim([10, 50])
    ax.set_yscale('log')

    # Scatter plots
    s1 = ax.scatter(plot_data['x'], plot_data['Flux_r_S'], color=colors['Flux_r_S'], label='Flux_r_S', s=60)
    s2 = ax.scatter(plot_data['x'], plot_data['Flux_reconnect'], color=colors['Flux_reconnect'], marker='x', label='Flux_reconnect', s=60)
    s3 = ax.scatter(plot_data['x'], plot_data['Flux_sb'], color=colors['Flux_sb'], marker='^', label='Flux_sb', s=60)

    # Assign numbers to names
    number_map = {name: i+1 for i, name in enumerate(plot_data['name'])}

    # Annotate points with numbers
    for xi, frs, name in zip(plot_data['x'], plot_data['Flux_r_S'], plot_data['name']):
        if pd.notna(xi) and pd.notna(frs) and np.isfinite(xi) and np.isfinite(frs):
            ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=90)

    # Axis labels & title
    ax.set_xlabel('x')
    ax.set_ylabel('Flux')
    ax.set_title(geom.replace('_', ' ').title())

    # Legend for scatter markers (top-right)
    ax.legend(loc='upper left', title='Flux Types')

    # Create inset axes for number → name mapping
    ax2 = fig.add_axes([0.75, 0.55, 0.2, 0.4])  # [left, bottom, width, height] in figure coords
    ax2.axis('off')
    for i, (name, num) in enumerate(number_map.items()):
        ax2.text(0, 1 - i*0.05, f"{num} : {name}", fontsize=8, va='top')

    plt.tight_layout()
    plt.savefig('OUTPUT/alltargetsflux_' + geom + '.pdf', bbox_inches='tight')
    plt.close()

'''





'''
for geom in geom_list:
    plot_data = flux_data[geom]
    
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    ax.set_ylim([1.001e-3, 1e2])
    ax.set_xlim([10, 50])
    #ax.set_xlim([np.min(plot_data['x'])*0.9, np.max(plot_data['x'])*1.5])
    #ax.set_xscale('log')
    ax.set_yscale('log')
    
    ax.scatter(plot_data['x'], plot_data['Flux_r_S'], color=colors['Flux_r_S'], label='Flux_r_S', s=60)
    ax.scatter(plot_data['x'], plot_data['Flux_reconnect'], color=colors['Flux_reconnect'], marker='x', label='Flux_reconnect', s=60)
    ax.scatter(plot_data['x'], plot_data['Flux_sb'], color=colors['Flux_sb'], marker='^', label='Flux_sb', s=60)
    

    for xi, frs, name in zip(plot_data['x'], plot_data['Flux_r_S'], plot_data['name']):
        if pd.notna(xi) and pd.notna(frs) and np.isfinite(xi) and np.isfinite(frs):
            ax.text(xi, frs, name, fontsize=8, ha='right', va='bottom', rotation=90)
    
    ax.set_xlabel('x')
    ax.set_ylabel('Flux')
    ax.set_title(geom.replace('_', ' ').title())
    ax.legend()
    
    plt.tight_layout()
    #plt.show()
    plt.savefig('OUTPUT/alltargetsflux_'+geom+'.pdf')
'''







































'''
# Loop over geometries for plotting

# Create a new blank figure
plt.close('all')
#figure, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
figure, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
plt.subplots_adjust(hspace=0.4)  # increase vertical spacing between subplots
plt.tight_layout()
geom_list = ['open_parker_spiral', 'closed_dipole', 'pfss']
colors = {
    'Flux_r_S': 'orange',
    'Flux_reconnect': 'blue',
    'Flux_sb': 'green'
}

for ax, geom in zip(axes, geom_list):
    data = flux_data[geom]
    ax.set_ylim([1.001e-3,1e2])
    ax.set_xscale('log') 
    ax.set_yscale('log') 
    ax.scatter(data['x'], data['Flux_r_S'], color=colors['Flux_r_S'], label='Flux_r_S', s=60)
    ax.scatter(data['x'], data['Flux_reconnect'], color=colors['Flux_reconnect'], marker='x', label='Flux_reconnect', s=60)
    ax.scatter(data['x'], data['Flux_sb'], color=colors['Flux_sb'], marker='^', label='Flux_sb', s=60)
    
    # Add text labels
    for xi, frs, name in zip(data['x'], data['Flux_r_S'], data['name']):
        ax.text(xi, frs, name, fontsize=8, ha='right', va='bottom',rotation=90)
    
    ax.set_ylabel('Flux')
    ax.set_title(geom.replace('_', ' ').title())
    ax.legend()
    
ax.set_xlim([1,100])
axes[-1].set_xlabel('x')
plt.tight_layout()
plt.show()
'''







'''
figure = plt.figure(figsize=(10,6))  # new blank figure

colors = {
    'open_parker_spiral': 'orange',
    'closed_dipole': 'blue',
    'pfss': 'green'
}

for geom, data in flux_data.items():
    plt.scatter(data['x'], data['Flux_r_S'], color=colors[geom], label=f"{geom} - Flux_r_S")
    plt.scatter(data['x'], data['Flux_reconnect'], marker='x', color=colors[geom], label=f"{geom} - Flux_reconnect")
    plt.scatter(data['x'], data['Flux_sb'], marker='^', color=colors[geom], label=f"{geom} - Flux_sb")

plt.xlabel('x')
plt.ylabel('Flux')
plt.legend()
plt.tight_layout()
plt.show()
'''








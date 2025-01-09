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
dipole_mag_pl_lists   = table_lists()
dipole_unmag_pl_lists = table_lists()
spiral_mag_pl_lists   = table_lists()
spiral_unmag_pl_lists = table_lists()

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
    print('All table planets')
    print(data['planet_name'])
    planet_array = range(len(data))
else:
    # else we deal with one single target 
    planet_array = [0] 

#print(data)
#print(Data2)

for indi in planet_array:
# Read parameters from a table containing multiple targets or from a single target file
    if INPUT_TABLE == True: 
        #starname,d, R_star, M_star, P_rot_star, B_star, Exoplanet, Mp, Rp, r_orb, P_orb,eccentricity, q, Q = load_target(data, indi)
        starname,d, R_star, M_star, P_rot_star, B_star, Exoplanet, Mp, Rp, r_orb, P_orb,eccentricity, q, Q = load_target(data, indi)
        #print('starname,d, R_star, M_star, P_rot_star, B_star, Exoplanet, Mp, Rp, r_orb, P_orb,eccentricity, q, Q');print(starname,d, R_star, M_star, P_rot_star, B_star, Exoplanet, Mp, Rp, r_orb, P_orb,eccentricity, q, Q)
        T_corona = data['T_corona(MK)'][indi]*1e6
        M_star_dot = data['mdot(mdot_sun)'][indi]
    else:   # Read from target.py (for individual targets)
        file_name = 'INPUT.target' 
        # import all parameters into a single variable
        imported_parameters = importlib.import_module(file_name) 

        # convert imported_parameters into global variables
        globals().update({k: v for k, v in imported_parameters.__dict__.items() if not k.startswith('__')})
        
        # LPM  -  MODIFY PRINTING
        # print(starname, d, R_star, M_star, P_rot_star, B_star, Exoplanet, Mp, Rp, r_orb, P_orb)

    # Fill B_star column if empty. Uses original units from table
    if pd.isna(B_star):
         B_star= spi.B_starmass(star_mass=data['mass_star(m_sun)'][indi],Prot=data['p_rot(days)'][indi])
         data['bfield_star(gauss)'][indi]=B_star
         # Eventually include uncertainties in B _star
         data['e_bfield_star'][indi]='TBD'

    #B_star = data['bfield_star(gauss)'][indi]/R_SPI**3    # Stellar surface magnetic field
    B_star = B_star * (R_SPI)**-3                        # Magnetic field where the SPI emission takes place (R_SPI)  
    #nu_ecm = 2.8e6 * B_star # cyclotron freq, in Hz

    # cyclotron (= ECM) frequency, and ECM bandwith, in Hz
    nu_ecm = e*B_star/(2*np.pi * m_e * c) 
    Delta_nu_cycl = nu_ecm 

    #  Check whether M_star_dot is read from input table/file
    if np.isnan(M_star_dot):
        if INPUT_TABLE == True:       
            data['M_star_dot(M_sun_dot)'][indi] = spi.Mdot_star(R_star=data['radius_star(r_sun)'][indi], M_star=data['mass_star(m_sun)'][indi], Prot_star=data['p_rot(days)'][indi])/M_sun_dot
            #M_star_dot = data['M_star_dot(M_sun_dot)'][indi]     # Stellar mass loss rate in solar units 
            M_star_dot = spi.Mdot_star(R_star=R_star/R_sun, M_star=M_star/M_sun, Prot_star=P_rot_star/day)/M_sun_dot
        else:
            M_star_dot = spi.Mdot_star(R_star=R_star/R_sun, M_star=M_star/M_sun, Prot_star=P_rot_star/day)/M_sun_dot
        print('M_star_dot: ',M_star_dot)
        print(type(M_star_dot))    
    #  Check whether M_star_dot is read from input table/file
    if np.isnan(T_corona):
        T_corona = T_CORONA_DEF
     

    print('Exoplanet name: ', Exoplanet)

    ###############################################

    # Electron gyrofrequency and ECM bandwidth 
    #gyrofreq = e*B_star/(2*np.pi * m_e * c) # in cgs units
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
        print('You asked to carry out a study of radio emission vs orbital separation: STUDY == "D_ORB" ')
        Nsteps = int(2*d_orb_max)
        d_orb  = np.linspace(1.02, d_orb_max, Nsteps) * R_star # Array of (orbital) distances to the star, in cm 
        M_star_dot_arr = np.array([M_star_dot]) # Convert to a numpy array of 1 element for safety reasons
    elif STUDY == "M_DOT":
        print('You asked to carry out a study of radio emission vs Mass loss rate: STUDY == "M_DOT" ')
        d_orb  = np.array([r_orb])
        Nsteps = int(M_DOT_STRETCH * np.log10(M_DOT_MAX/M_DOT_MIN))
        M_star_dot_arr = np.logspace(np.log10(M_DOT_MIN), np.log10(M_DOT_MAX), Nsteps)
    elif STUDY == "B_PL":
        print('You asked to carry out a study of radio emission vs B_planet: STUDY == "B_PL" ')
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
    angle_v_rel = np.arctan(v_orb/v_sw)  # in radians
    
    # Compute n_sw_planet, the number density of the wind at the orbital distance to the planet. 
    # If the stellar plasma is assumed to be isothermal, then 
    # the density falls down as ~ R^(-2) * v_sw^(-1).
    # Alternatively, we fix the density at the distance of the planet from the host star.
    if isothermal:
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
            if Bfield_geom_arr[ind]:
                print("\nOpen Parker magnetic field geometry")
            else:
                print("\nClosed dipolar magnetic field geometry")
           
            # get magnetic field components
            B_r, B_phi, B_sw, angle_B, theta, geom_f = spi.get_bfield_comps(Bfield_geom_arr[ind], B_star, d_orb, R_star, v_corot,
                    v_sw, angle_v_rel)
            
            # Compute Alfvén parameters in the stellar wind at a distance d_orb 
            v_alf, M_A, v_alf_r, M_A_radial = spi.get_alfven(rho_sw_planet, B_sw, B_r, v_rel, v_sw)

            # defines whether planet is unmagnetized (magnetized_pl_arr[ind1] = 0), or magnetized (magnetized_pl_arr[ind1] = 1)
            if magnetized_pl_arr[ind1]: # magnetized planet
                print('Magnetized planet\n')
                if B_planet_law == 'Sano':
                    # Planetary magnetic field, using Sano's (1993) scaling law, in units of B_earth 
                    # Assumes a tidally locked planet, i.e., the rotation period of the
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
                print('Unmagnetized planet\n')
                B_planet_arr  = np.zeros(len(d_orb)) # unmagnetized planet

            if STUDY == "B_PL":
                B_planet_Sano = B_planet_arr # Planet magnetic field at r_orb. 1-element array, in Gauss. 
                B_planet_arr  = np.linspace(B_PL_MIN, B_PL_MAX, Nsteps)
            
            # Effective radius of the obstacle, in cm
            # Case 1. À la Saur+2013. 
            R_planet_eff_Saur = spi.get_Rmp_Saur(Rp, THETA_M, B_planet_arr, B_sw)


            # Case 2. À la Zarka (2007), Turnpenney+2018, etc.
            #
            # Compute radius of magnetopause, Rmp as balance of wind and planet's
            # pressures
            
            # Planet pressure, in erg/cm3 - only the magnetic component is considered
            P_B_planet  = spi.get_P_B_planet(B_planet_arr) 

            # Stellar wind pressure, in erg/cm3
            P_sw, P_dyn_sw, P_th_sw, P_B_sw = spi.get_P_sw(n_sw_planet, v_rel, T_corona, B_sw, mu)

            # Radius of magnetopause, in cm
            Rmp = spi.get_Rmp(P_B_planet, P_dyn_sw, P_th_sw, P_B_sw) * Rp
    
            # The effective radius (in cm) is the radius of the magnetopause
            R_planet_eff = Rmp
            R_planet_eff_normalized = R_planet_eff/Rp 

            # Find value of Bp where R_planet_eff where is larger than Rp
            indices_R_planet_eff_larger_Rp = np.argwhere(R_planet_eff > Rp)
            if indices_R_planet_eff_larger_Rp.size > 0:
                B_planet_eff_rad = B_planet_arr[indices_R_planet_eff_larger_Rp[0]]          
                print('value of Bp where magnesphere is larger than Rp: ',B_planet_eff_rad)

            # If R_pl_eff < R_p, force R_pl_eff = R_p (R_planet_eff cannot be smaller
            # than Rp)
            R_planet_eff[R_planet_eff < Rp] = Rp 
            
            ## Get Poynting fluxes for the Alfvén wing model (Zarka 2007, Saur 2013,
            # Turnpenney+2018)
            
            # Get Poynting flux using Eq. 55 in Saur+2013 (S_poynt) and Eq. 8 in Zarka
            # 2007 (S_poyn_Z), in cgs units. They coincide, except for a factor 2. 
            # In mks units
            S_poynt, S_poynt_Z = spi.get_S_poynt(R_planet_eff, B_sw, v_alf, v_rel, M_A, ALPHA_SPI, geom_f)

            # Get fluxes at Earth, in cgs units for both Saur+ (Flux_r_S...) and
            # Zarka/Lanza (Flux_r_S_Z...),  in erg/s/Hz/cm2
            Flux_r_S_min, Flux_r_S_max =  spi.get_Flux(Omega_min, Omega_max, Delta_nu_cycl, d, S_poynt)
            Flux_r_S_Z_min, Flux_r_S_Z_max =  spi.get_Flux(Omega_min, Omega_max, Delta_nu_cycl, d, S_poynt_Z)
            # Compute flux density for an intermediate value of eps (in log scale)
            Flux_r_S_inter = 10**((np.log10(Flux_r_S_max) + np.log10(Flux_r_S_min))/2)
            
            
            # Get flux for the reconnection model (Lanza 2009, A&A)
            S_reconnect, P_d, P_d_mks = spi.get_S_reconnect(R_planet_eff, B_sw, v_rel, gamma = 0.5)
            
            Flux_reconnect_min, Flux_reconnect_max = spi.get_Flux(Omega_min, Omega_max, Delta_nu_cycl, d, S_reconnect)
            Flux_reconnect_inter = 10**((np.log10(Flux_reconnect_max) + np.log10(Flux_reconnect_min))/2)


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
            Flux_r_S_min_no_abs=Flux_r_S_min
            Flux_r_S_max_no_abs=Flux_r_S_max
            Flux_r_S_inter_no_abs=Flux_r_S_inter
            
            Flux_reconnect_min_no_abs = Flux_reconnect_min
            Flux_reconnect_max_no_abs = Flux_reconnect_max
            Flux_reconnect_inter_no_abs = Flux_reconnect_inter

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
             
            """
            Moving parts of plotting outside the loop
            """
            # Find out the position of the planet in the distance array
            d_diff = np.abs((d_orb-r_orb)/R_star)
            loc_pl = np.where(d_diff == d_diff.min())

            M_star_dot_diff = np.abs(M_star_dot_arr - M_star_dot)
            M_star_dot_loc  = np.where(M_star_dot_diff == M_star_dot_diff.min())

            ###########################################################################
            ####                  PLOTTING                                         ####
            ###########################################################################
            # Create OUTPUT folder if it doesn't exist
            FOLDER = 'OUTPUT/' + str(Exoplanet.replace(" ", "_"))
            if not(os.path.isdir(FOLDER)):
                os.system('mkdir OUTPUT/' + str(Exoplanet.replace(" ", "_")))
            else:
                print(FOLDER + ' already exists.')
                
                
            ### Plot received flux density as a function of distance from the star
            filename = 'plotting/plot_flux_density.py'
            with open(filename) as file:
                exec(file.read())            

            
            if freefree == True and STUDY == 'M_DOT': #################
                plt.figure(figsize=(8,11))
                ax = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
                ax.plot(M_star_dot_arr, absorption_factor, color='k')
                ax.set_xscale('log')
                ax.set_xlabel(r"Mass Loss rate [$\dot{M}_\odot$]",fontsize=20)
                ax.set_ylabel(r"Fraction of transmitted flux")
                ax.text(1e-1, 0, r'T$_{c} = $'+"{:.1f}".format(T_corona/1e6)+' MK', fontsize = 22)
                ax.set_facecolor("white")	
                plt.savefig(FOLDER + '/' + str(Exoplanet.replace(" ", "_"))
                        +'-'+'absorption_vs_mdot'+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'-'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'.pdf')
                plt.close()
                
            #### Plot effective radius variation
            filename = 'plotting/plot_effective_radius.py'
            with open(filename) as file:
                exec(file.read())
    
            
            # LPM  - REVISE
            li = [x,y_min.tolist(), y_max.tolist()]   
            
              
            #df = pd.DataFrame(data=li)
            #df = df.assign(column_name=column_series)
            #df.index = [STUDY, 'flux_min'+str(T_corona/1e6)+'MK', 'flux_max'+str(T_corona/1e6)+'MK']
            
            if freefree == True: 
                df = pd.DataFrame(zip(x,y_min, y_max), columns=[STUDY, 'flux_min'+str(T_corona/1e6)+'MK', 'flux_max'+str(T_corona/1e6)+'MK'])
                df.to_csv(os.path.join(outfile + ".csv"))            
                df2= pd.DataFrame(zip(x,absorption_factor),columns=[STUDY,'abs_factor_'+str(T_corona/1e6)+'MK'])
                df2.to_csv(FOLDER + '/' + str(Exoplanet.replace(" ", "_"))
                         +'-'+'absorption_vs_mdot'+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'-'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'.csv')   
            
            ################################
            # DIAGNOSTIC PLOTS
            ################################


            filename = 'plotting/plot_diagnostic_plots.py'
            with open(filename) as file:
                exec(file.read())
    
            ###########################################################
            ################### Send OUTPUT to external text file/s
            ###########################################################

            outfileTXT = os.path.join(outfile+'.txt')
            out_to_file = OutputWriter(outfileTXT)
            ### NOTE: Add B_planet_arr[loc_pl] in the output table
            print(" ")
            #print('n_base_corona = ', n_base_corona)
            #print('type(n_base_corona) = ', type(n_base_corona))
            print(" ")
            print('M_star_dot_loc = ', M_star_dot_loc)
            print('Type of M_star_dot_loc : ', type(M_star_dot_loc))
            print('n_base_corona[M_star_dot_loc] = ', n_base_corona[M_star_dot_loc])
            print('############################')
            #print(x_larger_rms)
            if x_larger_rms is np.nan or 'nan':
                print('NO value of '+STUDY+' where there is clear detection for the Alfvén Wing model')
            else:
                print('value of '+STUDY+' where there is clear detection for the Alfvén Wing model: ',x_larger_rms)
                
            out_to_file.write_parameters(T_corona, M_star_dot, mu, d, R_star, M_star, P_rot_star, B_star, 
                Exoplanet, Rp, Mp, r_orb, P_orb, loc_pl, M_star_dot_loc, n_base_corona,
                nu_plasma_corona, nu_ecm, Flux_r_S_min, Flux_r_S_max, rho_sw_planet, n_sw_planet, v_sw_base, Flux_r_S_Z_min,
                Flux_r_S_Z_max, v_sw, v_rel, v_alf, M_A, B_sw, Rmp, R_planet_eff,x_larger_rms,x_smaller_rms,STUDY,Omega_min, Omega_max,R_planet_eff_normalized,x_superalfv)

            # Print out the expected flux received at Earth from the SPI at the position of the planet

            print("\nPrint out minimum and maximum values of flux density at the planet location")
            print('B_planet_ref = {0:.3f} G'.format(B_planet_ref * bfield_earth*Tesla2Gauss))
            print("Saur/Turnpenney (mJy): ", Flux_r_S_min[loc_pl], Flux_r_S_max[loc_pl])
            print("Zarka: (mJy)", Flux_r_S_Z_min[loc_pl], Flux_r_S_Z_max[loc_pl])
            print("Reconnection: (mJy)", Flux_reconnect_min[loc_pl], Flux_reconnect_max[loc_pl])
            print(f"Done with planet {Exoplanet}")
            print(spi.Kepler_r(M_star/M_sun,P_orb))
            print(spi.Kepler_r(M_star/M_sun,np.array([2.34,4.12,6.74])))
            print(spi.Kepler_r(M_star/M_sun,np.array([2.34,4.12,6.74]))*au/R_star)
            print(str(np.round(spi.Kepler_r(M_star/M_sun,np.array([2.34,4.12,6.74]))*au/R_star, 2)))
            #print(B_planet_Sano)
            #### TEMPORARY TABLE
            ####################

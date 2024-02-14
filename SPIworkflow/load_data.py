import os
import pandas as pd
import numpy as np
import SPIworkflow.SPIutils as spi
from SPIworkflow.__init__ import *
from SPIworkflow.constants import *


def get_spi_data(infile_data='./INPUT/SPI-targets.csv',
        out_latex=True, outfile='latex_table.tex', 
        distance_min = 0.1, distance_max=1000.0, 
        p_orb_min=0.1, p_orb_max=1000.0, 
        bfield_min=10., bfield_max=10000., 
        dec_min = -90.0, dec_max=90.): 
    """
    Read in the input data to estimate radio emission from SPI

    Returns
    -------
    data : pandas.Dataframe
        The SPI useful data, after having applied different masks

    latex_table.tex: Latex table 
        Latex table of the data pandas.Dataframe (optional)

    Parameters
    ----------
    infile_data : string (optional)
        location to read in the data

    out_latex : Boolean (optional)
        If True, create table

    outfile   : string (optional)
        Name of the output latex table

    distance_min, distance_max : float (optional)
        Distance cuts (min and max), in pc

    p_orb_min, p_orb_max : float (optional)
        Orbital period cuts, in days

    bfield_min, bfield_max : float (optional)
        Magnetic field cuts, in Gauss

    dec_min, dec_max : float (optional)
        Declination cuts, in degrees
    """
    #df=df.rename(columns={"star_radius(R_Sun)": "radius_star(r_sun)", "Planet_radius(R_Earth)": "radius_planet(r_earth)"
    #                  , "P_orb(days)": "p_orb(days)", "M_star": "mass_star(m_sun)"
    #                   , "semi_major_axis(AU)": "a(au)", "RA(deg)": "ra(deg)"
    #                   , "DEC(deg)": "dec(deg)", "Planet_mass(M_Earth)": "mass_planet(m_earth)"
    #                   , "starname": "star_name", "P_rot(days)": "p_rot(days)"
    #                   , "<B>": "bfield_star(gauss)", "sptype": "SP_TYPE"
    #                   , "distance(pc)": "d_star(pc)"
    #                  })

    # Read in data file using pandas
    df = pd.read_csv(infile_data)

    # Create new pandas.Dataframe with a subset of rows (easy to display in,
    # e.g., jupyter notebook.
    if infile_data =='./INPUT/SPI-sources_planets_MASTER.csv':
        df2 = df[["star_name", "SP_TYPE", "d_star(pc)", "mass_star(m_sun)",
            "radius_star(r_sun)", "p_rot(days)", "bfield_star(gauss)",
            "planet_name", "p_orb(days)", "a(au)", "radius_planet(r_earth)",
            "mass_planet(m_earth)","ra(hms)","dec(dms)", "ra(deg)", "dec(deg)"]]

    elif infile_data =='./INPUT/my-SPI-sources.csv':
        df2 = df[["planet_name", "star_name", 
            "ra(hms)","dec(dms)", "ra(deg)", "dec(deg)", 
            "d_star(pc)", "radius_star(r_sun)", "mass_star(m_sun)", "p_rot(days)", "bfield_star(gauss)",
            "radius_planet(r_earth)", "mass_planet(m_earth)", "p_orb(days)", "a(au)" ]]
    else: 
        df2 = df[["planet_name", "star_name",  
            "d_star(pc)", "radius_star(r_sun)", "mass_star(m_sun)", "p_rot(days)", "bfield_star(gauss)", 
            "e_bfield_star", "radius_planet(r_earth)", "mass_planet(m_earth)",
            "mass_sini", "p_orb(days)", "a(au)" ,"eccentricity"]]

    #
    #df2.to_csv(r'./INPUT/SPI-sources.csv', index=True, header=True)

    # Define masks
    #
    #mask_Rpl = df2['radius_planet(r_earth)'] < 1.5
    #p_rot_mask = df2['p_rot(days)'] < 200.0
    distance_mask_min =  df2['d_star(pc)'] > distance_min
    distance_mask_max =  df2['d_star(pc)'] < distance_max 
    p_orb_mask_min = df2['p_orb(days)'] > p_orb_min
    p_orb_mask_max = df2['p_orb(days)'] < p_orb_max
    bfield_mask_min = df2['bfield_star(gauss)'] > bfield_min
    bfield_mask_max = df2['bfield_star(gauss)'] < bfield_max
    #
    #MPT@20230116 -  The masks in declination give a TypeError when reading
    # './INPUT/SPI-sources.csv' 
    # (but not when reading the file './INPUT/SPI-sources_planets_MASTER.csv'    
    # TypeError: '>' not supported between instances of 'str' and 'float'
    # Needs to be understood and fixed.
    # It may have to do with "NaN" values, which may conflict
    #
    #declination_mask_min = df2['dec(deg)'] > dec_min
    #declination_mask_max = df2['dec(deg)'] < dec_max

    # Apply masks
    #data = df2[distance_mask_min & distance_mask_max & p_orb_mask_min &
    #        p_orb_mask_max & bfield_mask_min & bfield_mask_max &
    #        declination_mask_min & declination_mask_max]
    data = df2[distance_mask_min & distance_mask_max & p_orb_mask_min &
            p_orb_mask_max & bfield_mask_min & bfield_mask_max]


    # Create column with cyclotron freq, in GHz
    # It gives a SetWarningCopy message, but it's seems to work fine.
    #
    #data['freq_cycl(ghz)'] = data.loc[:,('bfield_star(gauss)')]*2.8e-3
    data['freq_cycl(ghz)'] = data['bfield_star(gauss)']*2.8e-3

    data.reset_index(inplace=True)


    # Generate also latex_table? out_latex==True => Return. It will truncate 
    if out_latex==True:
        with open('./OUTPUT/latex_table.tex', 'w') as tf:
            tf.write(data.to_latex(index=False))
        #print(data.to_latex(index=False))

    return data

def create_data_tables(infile_data='./INPUT/SPI-targets.csv'):
    """
    Read in an input table with stars and, optional, planets, and
    creates two tables: one including targets hosting planets, and another
    one with targets that do not host planets.

    Returns
    -------
    df_planets : pandas.Dataframe
         It contains the targets with planets 

    df_no_planets : pandas.Dataframe
         It contains the targets without planets 

    outfile_planets: File, for the moment, it's fixed:  'StarTable-CARMENES_only_planets.csv')

    outfile_no_planets: File, for the moment, it's fixed:  'StarTable-CARMENES_no_planets.csv')

    Parameters
    ----------
    infile_data : string (optional)
        location to read in the data

    """

    # Create OUTPUT directory for keeping output from code.
    outdir = 'OUTPUT'
    try:
        os.mkdir(outdir)
        print('Directory', outdir, 'created.\n')
    except FileExistsError:
        print('Directory', outdir, 'already exists.\n')

    df = pd.read_csv(infile_data)
    # Copy Dataframes to generate two tables: one for sources with planets,
    # and a second one for sources not hosting planets
    df_planets=df.copy()
    df_planets=df_planets.dropna(subset=['planet_name'], inplace=False)

    df_no_planets = df.copy()
    df_no_planets = df_no_planets[df_no_planets['planet_name'].isnull()]

    # Create CARMENES tables for targets with planets only, and with no planets,
    # unless they already exist
    # 
    outfile_planets = os.path.join(outdir, 'StarTable-CARMENES_only_planets.csv')
    outfile_no_planets = os.path.join(outdir, 'StarTable-CARMENES_no_planets.csv')

    if os.path.exists(outfile_planets): 
        print(f'File {outfile_planets} already exists. I will not overwrite it.\n')
    else:
        print(f'Creating table file: {outfile_planets}.\n')
        df_planets.to_csv(outfile_planets)

    if os.path.exists(outfile_no_planets): 
        print(f'File {outfile_no_planets} already exists. I will not overwrite it.\n')
    else:
        print(f'Creating table file: {outfile_no_planets}.\n')
        df_no_planets.to_csv(outfile_no_planets)

    return outdir, df_planets, df_no_planets 

def load_target(data, indi):
    """Loads data for each target in input data frame
    OUTPUT:
    INPUT: 
    """
    starname=data['star_name'][indi]
    d      = data['d_star(pc)'][indi] * pc               # Distance to stellar system , in  cm
    R_star = data['radius_star(r_sun)'][indi] * R_sun    # Stellar radius in cm
    M_star = data['mass_star(m_sun)'][indi] * M_sun      # Stellar mass in g,
    P_rot_star = float(data['p_rot(days)'][indi]) * day  # Rotation period  of star, in sec
    B_star = data['bfield_star(gauss)'][indi]            # Stellar surface magnetic field

    # Read info for planets in table
    if source_data == './INPUT/SPI-targets.csv':
        # Planet - 
        Exoplanet = data['planet_name'][indi]
        Mp = float(data['mass_planet(m_earth)'][indi])*M_earth # Planetary mass, in grams
        if pd.isna(Mp): # If there is no mass value, use mass * sin(i)
            Mp = data['mass_sini'][indi] * M_earth 
        Rp = data['radius_planet(r_earth)'][indi]
        if pd.isna(Rp): 
            Rp = spi.Rp_Zeng(data['mass_planet(m_earth)'][indi])
        Rp *= R_earth # Planetary radius, in cm
        P_orb = data['p_orb(days)'][indi] # orbital period of planet, in days
        r_orb  = data['a(au)'][indi] * au   # orbital distance, in cm
        if pd.isna(r_orb): # If there is no mass value, use mass * sin(i)
            r_orb = spi.Kepler_r(M_star/M_sun, P_orb) * au
        eccentricity=data['eccentricity'][indi]
    else:
        # If no planet, set exoplanet as Earth, and semi-major axis = 0.2 * au 
        Exoplanet = 'Earth'
        Mp = M_earth # Planetary mass, in grams
        Rp = R_earth # Planetary radius, in cm
        r_orb  = 0.2 * au    # orbital distance, in cm
        eccentricity=0
        P_orb = spi.Kepler_P(data['mass_star(m_sun)'][indi], 0.2)   #orbital period of planet, in days

    # compute periastron (q) and apoastron (Q), in cm
    q = (1 - eccentricity) * r_orb 
    Q = (1 + eccentricity) * r_orb

    return starname,d, R_star, M_star, P_rot_star, B_star, Exoplanet, Mp, Rp, r_orb, P_orb,eccentricity, q, Q


def define_lists():
    planet_name_list=[]
    star_name_list=[]
    d_star_list=[]
    mass_star_list=[]
    radius_star_list=[]
    p_rot_list=[]
    bfield_star_list=[]
    a_list=[]
    p_orb_list=[]
    eccentricity_list=[]
    q_list=[]
    Q_list=[]
    mass_planet_list=[]
    radius_planet_list=[]
    T_cor_list=[]
    m_dot_list=[]
    nu_pl_list=[]
    nu_cycl_list=[]
    rho_pl_list=[]
    B_pl_list=[]
    B_sw_list=[]
    v_alf_list=[]
    M_A_list=[]
    Flux_r_S_ZL_min_list=[]
    Flux_r_S_ZL_max_list=[]
    P_Bpl_list=[]
    P_dyn_list=[]
    P_th_list=[]
    P_Bsw_list=[]
    Rmp_list=[]
    
    return locals()

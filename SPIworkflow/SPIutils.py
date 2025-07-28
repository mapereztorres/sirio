import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

from setup import *
from SPIworkflow.constants import *
from scipy.special import lambertw
from scipy.optimize import fsolve
import math as math

# FUNCTION DEFINITIONS 

def get_velocity_comps(M_star, d_orb, P_rot_star):
    """ 
    OUTPUT: v_orb - Orbital/Keplerian speed of planet (array), in cm/s
            v_corot -  Corotation speeds (array), in cm/s
            Omega_star - Angular rotation velocity of the star (float), in s^(-1)
    INPUT:  M_star - Star mass (float), in gr
            d_orb  - orbital distances (array), in cm
            P_rot_star - Rotation period of star (float), in sec
    """
    v_orb = (G * M_star / d_orb)**0.5 
    Omega_star = 2.0*np.pi / P_rot_star 
    v_corot = d_orb * Omega_star 

    return v_orb, v_corot, Omega_star
 
def get_bfield_comps(Bfield_geom, B_star, d_orb, R_star, v_corot, v_sw, angle_v_rel):
    """
    Computes the radial and azimuthal components of the stellar wind magnetic field for
    three different geometries: 
    1.- An open Parker spiral (open_parker_spiral)
    2.- A pure (closed) dipole (closed_dipole)
    3.- A modified dipole (potential field source surface, PFSS). 
    4.- A hybrid model, using a sigmoid to model the transition from a pure closed
    dipole to an open Parker spiral magnetic field geometry. 

    OUTPUT:
    INPUT: 
        Bfield_geom (string): 
        B_star      (float) : surface magnetic field at the stellar equator, in Gauss
        d_orb       (array) : array of orbital distances to the planet, in cm
        R_star      (float) : Radius of star, in cm
        v_corot     (array) : array of corotation speeds, in cm/s
        v_sw        (array) : stellar wind speed, in cm/s
        angle_v_rel (array) : Angle of the relative velocity of the planet wrt the
                              stellar wind velocity (radians)
        R_alfven    (array) : Alfvén radius, in units of stellar radii
    """

    
    #R_SS = R_T = R_alfven
    R_T = R_SS
    Delta_R_norm = DELTA_R / R_star

    r_dipole_geom     = np.cos(DIPOLE_TILT) * np.cos(POLAR_ANGLE) + np.sin(DIPOLE_TILT) * np.sin(POLAR_ANGLE) * np.cos(AZIMUTH)
    theta_dipole_geom = np.cos(DIPOLE_TILT) * np.sin(POLAR_ANGLE) - np.sin(DIPOLE_TILT) * np.cos(POLAR_ANGLE) * np.cos(AZIMUTH)
    phi_dipole_geom   = np.sin(DIPOLE_TILT) * np.sin(AZIMUTH)

    # Magnetic moment of the star
    Magnetic_m_dipole =  R_star**3 * B_star
    
    # Pure dipole equations
    B_r_dipole     =  Magnetic_m_dipole / d_orb**3  * 2 * r_dipole_geom
    B_theta_dipole =  Magnetic_m_dipole / d_orb**3  * theta_dipole_geom
    B_phi_dipole   =  Magnetic_m_dipole / d_orb**3  * phi_dipole_geom
    
    
    ### Pure Parker spiral
    B_r_parker =  B_star * (d_orb/R_star)**(-2)
    B_theta_parker = 0.0
    B_phi_parker = B_r_parker * v_corot/v_sw * np.sin(POLAR_ANGLE)
    
    
    # PFSS model
    # Dipole magnetic field at the PFSS
    B_r_dipole_pfss     =  Magnetic_m_dipole / (R_SS*R_star)**3  * 2 * r_dipole_geom
    B_theta_dipole_pfss =  Magnetic_m_dipole / (R_SS*R_star)**3  * theta_dipole_geom
    B_phi_dipole_pfss   =  Magnetic_m_dipole / (R_SS*R_star)**3  * phi_dipole_geom
    B_sw_pfss  = np.sqrt(B_r_dipole_pfss**2 + B_theta_dipole_pfss**2 + B_phi_dipole_pfss**2)
        
    #Parker spiral equations imposing continuity of B_sw at the PFSS    
    B_continuity = B_sw_pfss * R_SS**2
    B_r_parker_pfss =  B_continuity * (d_orb/R_star)**(-2)
    B_theta_parker_pfss = 0.0
    B_phi_parker_pfss = B_r_parker_pfss * v_corot/v_sw * np.sin(POLAR_ANGLE)
    #print('B_r_parker_pfss)

    # pffs_r_factor and pfss_theta_factor taken as in Jardine+2002 (MNRAS). See Eqs. 7 and 8 
    #pfss_r_factor=  ((d_orb/R_star)**3 + 2*R_SS**3) / (1 + 2*R_SS**3) 
    #pfss_theta_factor =  (-2*(d_orb/R_star)**3 + 2*R_SS**3)/(1+2*R_SS**3) 
    #pfss_theta_factor tends to zero at the PFSS. 
    #on other hand the radial field tends to zero at the equator
    #so that the stellar wind magnetic field at the PFSS is not null
    
    #We set these two factors to 1 since we consider a classical dipole up to the PFSS.    
    pfss_r_factor = pfss_theta_factor = 1
 
    #The heaviside function is 0 if d_orb/R_star <= R_T, 1 if d_orb/R_star > R_T 
    heaviside = np.heaviside((d_orb/R_star)- R_T, 0)
    
    # The planet cannot be in the same axis of the dipolar magnetic moment of the star  
    # If that's the case, a Parker spiral geometry is enforced.
    #Currently not used, replaced by the heaviside function
    if (np.abs(AZIMUTH) < TOLERANCE and np.abs(POLAR_ANGLE - DIPOLE_TILT) < TOLERANCE) or (np.abs(AZIMUTH-np.pi) < TOLERANCE and np.abs(POLAR_ANGLE - (np.pi - DIPOLE_TILT)) < TOLERANCE):
        sigmoid = 1.0
        print('NOTE: The combination of selected values for AZIMUTH, POLAR_ANGLE and DIPOLE TILT are such that \nthe planet is located along the axis of the dipolar magnetic moment of the star.')
        print('Therefore, there are no closed (dipolar) lines, so the PFSS_Parker model will be equivalent to \nrunning a pure Parker spiral.')
        #else:
        #    print('planet and magnetic dipole are both in the XZ plane, but neither alligned nor anti-alligned')
    else:  
        sigmoid = 1 / (1 + np.exp(- ((d_orb/R_star) - R_T) / Delta_R_norm))
        
        

    
    if Bfield_geom == 'open_parker_spiral': 
        # Open Parker Spiral - Falls down with distances as R^(-2) rather than R^(-3) as in the dipole case
        # B_r   - Radial component of B_sw (e.g., Eqn 20 Turnpenney 2018).
        # B_phi - Azimuthal component of B_sw
        # B_theta - Polar component. Here we assume symmetry in the radial direction, so
        # B_theta = 0
        # B_sw - Total stellar wind B-field at the orbital distance of the planet
        B_r     = B_r_parker  
        B_theta = B_theta_parker 
        B_phi   = B_phi_parker

    elif Bfield_geom == 'closed_dipole':
        # Closed, dipolar configuration - It falls with distance as R^(-3)
        # B_star - magnetic field at the magnetic equator on the stellar surface
        # B_r - Radial component of the dipole magnetic field of the stellar wind as f(distance to star) 
        # B_theta_dipole - polar component of the dipole magnetic field 
        # B_phi_dipole   - azimuthal component of the dipole magnetic field (0 for a
        # dipole with spherical symmetry)
        # B_sw - Total stellar wind B-field at planet orbital distance
        B_r     = B_r_dipole
        B_theta = B_theta_dipole
        B_phi   = B_phi_dipole

    elif Bfield_geom == 'pfss': 
        # PFSS with a transition to an open Parker spiral
        #
        # It is a classical dipole up to the PFSS and a modified Parker spiral beyond
        
        B_r     = B_r_dipole * pfss_r_factor * (1 - heaviside) + B_r_parker_pfss * heaviside 
        B_theta = B_theta_dipole *pfss_theta_factor *(1 - heaviside) + B_theta_parker_pfss * heaviside
        B_phi   = B_phi_dipole* (1 - heaviside) + B_phi_parker_pfss * heaviside

    # Total stellar wind magnetic field
    B_sw  = np.sqrt(B_r**2 + B_theta**2 + B_phi**2) 

    # Eq. 23 of Turnpenney 2018 -  First term of RHS
    #angle_B = np.arctan(B_phi/B_r) # Angle the B-field makes with the radial direction
    print(type(B_r),type(B_phi))
    angle_B = np.arctan2(B_phi, B_r) # Angle the B-field makes with the radial direction
    # Angle between the stellar wind magnetic field and the impinging plasma velocity
    # Eqn 23 in Turnpenney 2018. It's also Eq. 13 in Zarka 2007
    theta = np.absolute(angle_B - angle_v_rel) # theta is the angle between the B_sw (the insterstellar magnetic field), and the
    # incident stellar wind velocity.  See Fig. 1 in Turnpenney+2018
    #
    
    if Bfield_geom == 'closed_dipole':
        theta = np.pi/2 # np.pi/2 for closed dipole configuration, different for the open field configuration
        #geom_f = (np.sin(theta))**2 # Geometric factor is 1 for the dipolar configuration
        #print('geom_f for dipolar configuration: ',geom_f)
 
    #if Bfield_geom == 'open_parker_spiral': 
        # theta is the angle between the B_sw (the insterstellar magnetic field), and the
        # incident stellar wind velocity.  See Fig. 1 in Turnpenney+2018
        #
        #geom_f = (np.sin(theta))**2 # Geometric factor in efficiency 
       
    geom_f = (np.sin(theta))**2
    
         
    if Bfield_geom == 'pfss': 
        geom_f = (1 - heaviside)+(np.sin(theta))**2 * heaviside # Geometric factor in efficiency 
    #geom_f = 1.0
    return B_r, B_phi, B_sw, angle_B, theta, geom_f

def getImage(path):
    return OffsetImage(plt.imread(path, format="jpg"), zoom=.02)

def Kepler_r(M_star, P_orb):
    """Computes the orbital radius of a planet (=semimajor axis, in units of au)
       orbiting a star of mass M_star (in solar masses) and orbital period, P_orb (in days)
    """
    r_orb = (G * M_star*M_sun)**(1/3) * (P_orb*day/2/np.pi)**(2/3)/au
    return r_orb

def Kepler_P(M_star, r_orb):
    """Computes the orbital period (P_orb) of a planet (in days) of a planet
       orbiting a star of mass M_star (in solar masses) and semi-major axis (r_orb), in au.
    """
    
    P_orb = 2*np.pi * (G * M_star*M_sun)**(-1/2) * (r_orb*au)**(3/2)/day
    return P_orb

def beta_keV(E_kin=511):
    """ Computes the velocity of an electron with kinetic energy E_k, in units of the speed of light
        OUTPUT: beta (= particle speed / speed of light)
        INPUT: E_kin - kinetic energy of the particle, in keV
    """ 
    E_0 = 511.0 # Rest mass of the electron, in keV
    beta = np.sqrt(1 - (1 + E_kin/E_0)**(-2)) 
    return beta

def Rp_Mueller(Mp = 1.0):
    """ Computes the planet radius using the mass-radius relationship in 
        Mueller+2019 (A&A, 2024): https://www.aanda.org/10.1051/0004-6361/202348690
        
        OUTPUT: Rp - planet radius, in units of Earth radius 
        INPUT:  Mp - planet mass, in units of Earth mass
    """
    if Mp <= 0:
        raise ValueError("Planet mass (Mp) must be a positive value.")
    
    if Mp < 4.37:
        Rp = 1.02 * Mp**0.27
    elif 4.37 <= Mp < 127:
        Rp = 0.56 * Mp**0.67
    else:
        Rp = 18.6 * Mp**(-0.06)

    return Rp

def Rp_Zeng(Mp = 1.0):
    """ Computes the planet radius using the mass-radius relationship in 
        Zeng+2016 (ApJ). https://iopscience.iop.org/article/10.3847/0004-637X/819/2/127
        See also Zeng et al. (PNAS, 2019): https://www.pnas.org/content/116/20/9723
        R [R_earth] = f * (M/M_earth)^(1/3.7)  (for rocky planets, f=1.)
        
        OUTPUT: Rp - planet radius, in units of Earth radius 
        INPUT:  Mp - planet mass, in units of Earth mass
        We assume f = 1 (Earth-like rocky planet)
    """
    f = 1.0
    Rp = f * Mp**(1/3.7)
    return Rp

def bfield_sano(M_planet = 1.0, R_planet = 1.0, Omega_rot_planet = 1.0):
    """ Computes the surface magnetic field strength of a planet using the Sano scaling law.
        (From Sano, J. Geomag. Geolectr., 45, 65-77, 1993), 
        OUTPUT: B_planet: Surface magnetic field of the planet, in units of B_Earth 
        INPUT : M_planet: Planet radius, in units of the Earth mass
                R_planet: in units of the Earth radius
                Omega_rot_planet: Rotational speed of the planet. It is assumed that the 
                        planets are tidally locked, hence P_rot_planet = P_orb_planet. 
                        In units of the rotational speed of Earth (radians/sec) 
                rho_core - Planet outer core density, in g/cm^3. 
    """
    # Scaling law for the planet's core radius, from Curtis & Ness (1986)
    r_core = M_planet**0.44 # in units of r_core_Earth
    rho_core = M_planet/R_planet**3
    magn_moment_planet = Omega_rot_planet * rho_core**(1/2) * r_core**(7/2) # Magnetic moment, in units of Earth magn. moment
    B_planet  = magn_moment_planet / R_planet**3 # in units of B_earth
 
    return r_core, rho_core, magn_moment_planet, B_planet

def Lrad_leto(B_star=1.0, R_star=1.0, P_rot=1.0):
    """ Returns the radio luminosity of an early-type magnetic star, 
        using the scaling law of Leto et al. (2021, MNRAS)
        OUTPUT: Lnu_rad, in erg/s/Hz
        INPUT:  B_star, in Gauss
                R_star, in R_sun
                P_rot , in days
    """
    alpha = 13.6; beta = 1.94 # from Leto+2021
    nu_ecm = B_star * 2.8e6 # cyclotron freq, in Hz 
    Lnu_rad = 10**alpha * (B_star/1e3)**beta * R_star**(2*beta) * P_rot**(-beta) # erg/s/Hz
    Ltot_rad = Lnu_rad * nu_ecm
    return Ltot_rad
    
def v_stellar_wind(d_orb, M_star, T_corona, m_av):
    """ This function computes the sound speed (v_sound, in cm/s) and the sonic point (r_sonic) of
        an isothermal stellar wind. v_sound only depends on the temperature of the corona,
        T_corona. The sonic point, r_sonic, also depends on the mass of the star, M_star. 
        It also computes the stellar wind speed, in cm/s, at every distance (d_orb) from the star. 

        The stellar wind speed is computed as in Turnpenney+18, using the Lambert W
        function,  with D_r obtained using Eq. 18. In turn, Eq. 18 is taken from Eq. 15 in Cranmer 2004. 

       OUTPUT:  v_sound (cm/s) - Float: Sound speed of the isothermal stellar wind, in cm/s.
                v_sw  (cm/s)   - Array: Stellar wind speed at every distance in d_orb
       INPUT:   d_orb (cm)     - Array: Distances, measured from the center of the star
                M_star (g)     - Float: Mass of star, in grams
                T_corona (K)   - Float: Temperature at the base of the corona, in Kelvin
                m_av     (g)   - Float: Average particle mass of the stellar wind, in grams
    """
    v_sound = np.sqrt(k_B * T_corona / m_av) 
    r_sonic =  G * M_star / (2 * v_sound**2) # Radius of sonic point, in cgs units

    D_r = (d_orb/r_sonic)**(-4) * np.exp(4*(1 - r_sonic/d_orb) - 1)
    v_sw2 = np.zeros(len(d_orb), dtype=complex)
    v_sw  = np.zeros(len(d_orb))

    for i in range(len(d_orb)):
        if (d_orb[i]/r_sonic) >= 1.0:
            v_sw2[i] = -v_sound**2 * lambertw(-D_r[i], k=-1)
        else: 
            v_sw2[i] = -v_sound**2 * lambertw(-D_r[i], k=0)
        # The actual speeed is the real part of v_sw2[i]
        v_sw[i]  = np.sqrt(v_sw2[i].real)

    return v_sound, r_sonic, v_sw

def get_v_sw_terminal(R_star, M_star, T_corona, m_av):
    """ This function computes the terminal speed of an isothermal stellar wind
        The stellar wind speed is computed as in Turnpenney+18, using the Lambert W
        function,  with D_r obtained using Eq. 18. In turn, Eq. 18 is taken from Eq. 15 in Cranmer 2004. 

       OUTPUT:  v_sw_terminal (cm/s) - Float: terminal speed of the isothermal stellar wind, in cm/s.
       INPUT:   R_star (cm)          - Float Radius of the star
                M_star (g)           - Float: Mass of star, in grams
                T_corona (K)         - Float: Temperature at the base of the corona, in Kelvin
                m_av     (g)         - Float: Average particle mass of the stellar wind, in grams
    """
    v_sound = np.sqrt(k_B * T_corona / m_av) 
    r_sonic =  G * M_star / (2 * v_sound**2) # Radius of sonic point, in cgs units

    # Terminal speed -- computed at 200 stellar radii 
    r_terminal = 200 * R_star

    D_r = (r_terminal / r_sonic)**(-4) * np.exp(4*(1 - r_sonic/r_terminal) - 1)
    #v_sw2 = np.zeros(len(r_terminal), dtype=complex)
    #v_sw  = np.zeros(len(r_terminal))

    v_sw_terminal2 = -v_sound**2 * lambertw(-D_r, k=-1)
    # The actual speeed is the real part 
    v_sw_terminal  = np.sqrt(v_sw_terminal2.real)

    return v_sw_terminal


def get_eta_star(B_star, R_star, M_star_dot_arr, v_sw_terminal):
    """
    Compute the magnetic confinement parameter at the stellar equator.
    OUTPUT: 
        eta_star (array):  Magnetic confinement parameter at the equator of the stellar
                           surface. Adimensional. As defined in ud-Doula & Owocki (2002, ApJ).
                           It returns an array of the same length as M_star_dot_arr.

    INPUT: 
        B_star (float)          - stellar magnetic field in the equator (Gauss)
        R_star (float)          - Stellar radius, in cm
        M_star_dot_arr (array)  - Stellar mass loss rate,  in units of M_sun_dot
        v_sw_terminal  (float)  - terminal wind speed, in cm/s
    """

    #Convert M_star_dot_arr to cgs units
    M_star_dot_cgs = M_star_dot_arr * M_sun_dot * (M_sun / yr2sec) 

    # magnetic confinement parameter
    # B_star = B_0 / 2 (as defined in ud-Doula & Owocki 2002)
    eta_star = B_star**2 * R_star**2 / (M_star_dot_cgs * v_sw_terminal)

    return eta_star

def get_R_alfven(eta_star, colatitude):
    """
    Computes the Alfvén radius for a given colatitude value 
    (colatitude = 90 deg implies at the magnetic equator)
    We follow the formalism in ud-Doula & Owocki (2002, ApJ)

    OUTPUT: 
        R_alfvén (array) - Alfvén radius at a certain value of the POLAR_ANGLE (theta ==
                           colatitude). 
                           It returns an array of the same length as M_star_dot_arr.
                           In units of stellar radii
    INPUT: 
        eta_star (array):  Magnetic confinement parameter at the equator of the stellar
                           surface. Adimensional. As defined in ud-Doula & Owocki (2002, ApJ).
        colatitude (float)  - colatitude at which to determine R_alfven, in radians
    """
    factor = 4 - 3 * np.sin(colatitude)**2  ## factor in RHS of Eq. 8 

    def equation(R_ratio, eta_star):
        return (R_ratio**(2*Q_DIPOLE-2) - R_ratio**(2*Q_DIPOLE-3)) - eta_star * factor


    R_alfven = fsolve(equation, R_ALFVEN_GUESS, args=(eta_star))[0]

    return R_alfven
    
def get_R_alfven_alt(eta_star, colatitude):
    """
    Computes the Alfvén radius for a given colatitude value 
    (colatitude = 90 deg implies at the magnetic equator)
    We follow the formalism in ud-Doula & Owocki (2002, ApJ)

    OUTPUT: 
        R_alfvén (array) - Alfvén radius at a certain value of the POLAR_ANGLE (theta ==
                           colatitude). 
                           It returns an array of the same length as M_star_dot_arr.
                           In units of stellar radii
    INPUT: 
        eta_star (array):  Magnetic confinement parameter at the equator of the stellar
                           surface. Adimensional. As defined in ud-Doula & Owocki (2002, ApJ).
        colatitude (float)  - colatitude at which to determine R_alfven, in radians
    """
    factor = 4 - 3 * np.sin(colatitude)**2  ## factor in RHS of Eq. 8 

    def equation(R_ratio, eta_star):
        return (R_ratio**(2*Q_DIPOLE-2) - R_ratio**(2*Q_DIPOLE-3)) - eta_star * factor
    R_alfven=[]
    for ind in range(len(eta_star)):
        R_alfven.append(fsolve(equation, R_ALFVEN_GUESS, args=(eta_star[ind]))[0])

    return R_alfven
    
def get_theta_A(R_alfven_pole):
    """
    OUTPUT: theta_A (array): colatitude at which the last closed magnetic line of the
            dipole closes the stellar surface. 
            See Eq. 9 in ud-Doula & Owocki (2002). Note that here R_alfven_pole is
            already normalized to R_star.
    INPUT: 
        R_alfvén_pole (array):  Alfvén radius at the pole of the star (theta = 0)
                                It returns an array.  In units of stellar radii
    """

    theta_A = np.arcsin(np.sqrt(1./R_alfven_pole))
    theta_A_deg = np.degrees(theta_A)
    latitude = np.abs(90. - theta_A_deg)
    
    return theta_A_deg, latitude


def wind_composition(X_p = 0.5):
    """ Computes fraction of electrons, mean "molecular weight" and average particle
    mass of the stellar wind, given the fraction of protons.
    OUTPUT: X_e - fraction of electrons of the stellar wind
            mu  - mean "molecular" weight of the stellar wind
            m_av - average particle mass of the stellar wind (in gr)
    INPUT : X_p - fraction of protons
    """
    X_e = 1 - X_p
    mu = (X_p * m_p + X_e * m_e)/(m_p + m_e) 
    m_av = mu * m_p  

    return X_e, mu, m_av

def n_wind(M_star_dot=1.0, d=7e10, v_sw=25.6e5, m_av = 0.5 * m_p):
    """ Computes the particle density of the stellar wind at some distance d from the
        center of the star.
        OUTPUT: n_sw - particle density of the stellar wind at distance d, in #/cm^3
        INPUT:  M_star_dot - stellar mass-loss rate, in units of the Sun mass-loss rate
                d          - distance from the center of the star, in cm 
                v_sw       - Speed of stellar wind at distance d, in cm/s
                m_av       - average particle mass of the stellar wind (= mu * m_p) in grams
    """
    
    ### WARNING: when using the " *= " format below, the original value of M_star_dot is
    # overwritten outside the function, even if not present in return!!
    #
    #M_star_dot *= M_sun_dot * M_sun/yr2sec  # mass-loss rate, in grams/sec
    #rho  = M_star_dot/(4*np.pi * d**2 * v_sw) # Density of the stellar wind, in gr/cm3
    M_star_dot_grams = M_star_dot * M_sun_dot * M_sun/yr2sec  # mass-loss rate, in grams/sec
    rho  = M_star_dot_grams / (4*np.pi * d**2 * v_sw) # Density of the stellar wind, in gr/cm3
    n_sw = rho / m_av

    return n_sw

def get_alfven(rho_sw_planet, B_sw, B_r, v_rel, v_sw):
    """Computes Alfvén parameters in the stellar wind at a distance d_orb (though
    n_sw_planet, B_sw, v_rel).
       OUTPUT
       INPUT
    """
    # Alfven speed and Mach Number
    #v_alf = 2.18e11 * B_sw / np.sqrt(n_sw_planet) # Alfven speed at the distance of the planet, in cm/s
    #v_alf = B_sw / np.sqrt(4.0 * np.pi * rho_sw) 
    # Relativistically corrected Alfvén speed, as v_alf must be less than the speed of light
    v_alf = B_sw / np.sqrt(4.0 * np.pi * rho_sw_planet) * 1./np.sqrt(1 + (B_sw**2/(4.*np.pi * rho_sw_planet * c**2)))
    M_A   = v_rel/v_alf # Alfven mach number
    
    #Radial Alfvén speed
    mu_0_cgs = 1.0 # magnetic permeability in vacuum, in cgs units
    v_alf_r = B_r / np.sqrt(mu_0_cgs * rho_sw_planet) # in cm/s
    M_A_radial = np.abs(v_sw / v_alf_r)

    return v_alf, M_A, v_alf_r, M_A_radial

def plasma_freq(n_e = 1.0):
    """ Returns the plasma frequency.
        Any wave emitting at a frequency less than the plasma frequency will
        not propate in this medium.

        OUTPUT: 
              nu_plasma - plasma frequency, in Hz

        INPUT : n_e - electron density, in cm^(-3)
    """
    nu_plasma = 9.0e3 * np.sqrt(n_e)
    return nu_plasma

def get_Rmp_Saur(Rp, THETA_M, B_planet_arr, B_sw):
    """It computes the effective radius, R_obs, of the Alfvén wing, in cm, using Eq. 57 in 
       Saur+2013, A&A).  It depends on the orientation, THETA_M, of the intrinsic planetary
       magnetic field (B_planet_arr) wrt the external magnetic field of the stellar wind (B_sw).
    OUTPUT: R_obs (cm) - Array: Effective planet radius, in cm
    INPUT : Rp           (cm) - Float: Planet radius, in cm
            THETA_M      (rad)- Float: Angle of the planetary magnetic field wrt stellar
                                wind magnetic field, in radians
            B_planet_arr (G)  - Array: Planetary magnetic field, in Gauss
            B_sw         (G)  - Array: Stellar wind magnetic field, in Gauss
    """
    Rmp = Rp * np.sqrt(3*np.cos(THETA_M/2)) * (B_planet_arr/B_sw)**(1./3.) # in cm
    #R_obs[ R_obs < Rp] = Rp # R_obs cannot be smaller than Rplanet    

    return Rmp

def get_P_sw(n_sw, v_rel, T_e, B_sw, mu):
    """
       Computes the total pressure of the stellar wind exerted on the planet.

       OUTPUT: P_sw, P_dyn_sw, P_th_sw, P_B_sw - all in erg/cm3
               P_sw     - total pressure stellar wind at the orbital distance of the planet
               P_dyn_sw - dynamical pressure of the stellar wind at the orbital distance of the planet
               P_th_sw  - thermal pressure of the stellar wind at the orbital distance of the planet
               P_B_sw   - magnetic pressure of the stellar wind at the orbital distance of the planet

       INPUT : n_sw  - number density of the stellar wind at the planet's orbital position, in  # cm^(-3)
               v_rel - Relative velocity between the stellar wind to the orbital velocity of the planet, in cm/s
               T_e  - Temperature of the stellar wind ath the planet's position. 
               B_sw  - Stellar wind magnetic field, in Gauss
               mu    - mean molecular weight (adimensional)
    """

    rho_sw  = mu * m_p * n_sw
    P_dyn_sw = rho_sw * v_rel**2
    P_th_sw = n_sw * k_B * T_e
    P_B_sw = B_sw**2 / (8 * np.pi)

    P_sw = P_dyn_sw + P_th_sw + P_B_sw

    return P_sw, P_dyn_sw, P_th_sw, P_B_sw

def get_P_dyn_sw(n_sw=1e7, mu=0.5, v_rel=1e7):
    """Computes the dynamical pressure of the stellar wind at the orbital distance 
               of the planet, P_dyn_sw, in erg/cm3

       OUTPUT: P_dyn_sw, in erg/cm3

       INPUT : n_sw  - number density of the stellar wind at the planet's orbital position, in  # cm^(-3)
               mu    - mean molecular weight 
               v_rel - Relative velocity between the stellar wind to the orbital velocity of the planet, in cm/s
    """
    rho_sw  = mu * m_p * n_sw
    P_dyn_sw = rho_sw * v_rel**2

    return P_dyn_sw

def get_P_th_sw(n_sw=1e7, T_e=2e6):
    """Computes the thermal pressure of the stellar wind at the orbital distance 
                of the planet, P_th_sw, in erg/cm3

       OUTPUT: P_th_sw, in erg/cm3

       INPUT : n_sw - number density of the stellar wind at the planet's orbital position, in  # cm^(-3)
               k_B  - Boltzmann's constant
               T_e  - Temperature of the stellar wind ath the planet's position. 
    """
    P_th_sw = n_sw * k_B * T_e

    return P_th_sw

def get_P_B_sw(B_sw=1.0):
    """Computes the magnetic pressure of the stellar wind at the orbital distance 
               of the planet, P_dyn_sw, in erg/cm3

       OUTPUT: 
              P_B_sw, in erg/cm3

       INPUT : 
              B_sw  - Stellar wind magnetic field, in Gauss
    """

    P_B_sw = B_sw**2 / (8 * np.pi)

    return P_B_sw


def get_P_B_planet(B_planet = 1.0):
    """Computes the magnetic pressure of planet at the equator, P_Bp, in erg/cm3

       OUTPUT: 
              P_B_planet - magnetic pressure of planet at the equator, in erg/cm3

       INPUT: 
              B_planet - planetary magnetic field (dipole) at the poles, in Gauss

    """
    P_B_planet = (B_planet / 2)**2 / (8*np.pi)

    return P_B_planet

def get_Rmp(P_Bp=1.0, P_dyn_sw=1.0, P_th_sw=1.0, P_B_sw=1.0):
    """Computes the radius of the planetary magnetopause, in units of planet radii

    It's calculated as the radius that balances the pressure between the
    stellar wind (ram, thermal, and magnetic) and the planet (magnetic)
    pressure. Note that, in the equation below, the ram and thermal components
    of the planetary pressure are neglected.

    OUTPUT: Rmp - Radius of the planetary magnetopause, in units of Rp

    INPUT : P_Bp     - pressure of the planetary magnetic field (dipole) at the equator
            P_dyn_sw - dynamic pressure of the stellar wind at the orbital distance 
                       of the planet, in erg/cm3
            P_th_sw  - thermal pressure of the stellar wind at the orbital distance 
                       of the planet, in erg/cm3
            B_sw     - Stellar wind magnetic field, in Gauss
    """
    # Total pressure of the planet (only magnetic pressure is considered)
    P_planet = P_Bp

    # Total pressure of the stellar wind 
    P_sw = P_dyn_sw + P_th_sw + P_B_sw
    
    # Expression in Donati & Vidotto (2017).
    # Similar to others in, e.g., Zarka (2007), Turnpenney+2018, etc.
    # 
    Rmp = K_MAGNETOPAUSE**(1./3.) * (P_planet / P_sw)**(1./6) 

    return Rmp
 
 
    
def get_Reff_lanza( B_sw,B_pl):

    #Reff_lanza in R_pl units
    
    Reff_lanza=(B_sw/B_pl)**(-1/3)


    return Reff_lanza
    
        
    
    
def B_starmass(star_mass,Prot):
  """Calculation of the Stellar magnetic field at the surface as a function of the stellar mass and its rotation period
  The value of the magnetic field depends on the Rossby number (Ro), as defined in https://www.aanda.org/articles/aa/pdf/2022/06/aa43251-22.pdf.
  There are two relations, depending no whether the star rotates fast (Ro<0.13) or slow (Ro>0.13)
  
  The Rossby number in turn is calculated using the stellar mass, as defined in equation 6 of https://academic.oup.com/mnras/article/479/2/2351/5045253.
  
  
  OUTPUT: B_mass     -  Stellar magnetic field at the surface (in G)
  
  
  INPUT : star_mass  -  Mass of the star (solar units)
          Prot       -  Rotation period of the star (in days)
  
  
  
  """
  #https://www.aanda.org/articles/aa/pdf/2022/06/aa43251-22.pdf (table )
  #https://academic.oup.com/mnras/article/479/2/2351/5045253  (eq 6)
  C=2.33
  D=-1.50
  E=0.31
  alpha=-1.26
  alpha_fast=-0.11
  tau_mass=10**(C+D*star_mass+E*star_mass**2)
  Ro_mass=Prot/tau_mass
  if Ro_mass>0.13:
    print('Slow rotator')
    B_mass= 199*Ro_mass**alpha
  elif Ro_mass<0.13:
    print('Fast rotator')
    B_mass= 2050*Ro_mass**alpha_fast    
  else:
    print('Uncertain regime: Ro=',Ro_mass)
  return B_mass

def B_color(starname,star_mass,Prot):
  #in development
  from astroquery.simbad import Simbad
  Simbad.add_votable_fields('flux(V)')
  Simbad.add_votable_fields('flux(K)')
  #https://www.aanda.org/articles/aa/pdf/2022/06/aa43251-22.pdf (table )
  #https://academic.oup.com/mnras/article/479/2/2351/5045253  (eq 5)
  A=0.64
  B=0.25
  alpha=-1.26
  #eachstar=data['star_name'][indi]
  result_table = Simbad.query_object(starname)
  print(result_table['MAIN_ID'][0],result_table['FLUX_V'][0],result_table['FLUX_K'][0])
  #pdn['B_color'][ind]=result_table['FLUX_B'][0]
  V_color=result_table['FLUX_V'][0]
  K_color=result_table['FLUX_K'][0]  
  VK=V_color-K_color
  logtau_color=A+B*VK
  tau_color=10**(A+B*VK)
  Ro_color=Prot/tau_color
  B_color= 199*Ro_color**alpha

  return B_color


def Mdot_star(R_star, M_star, Prot_star):
    """
    OUTPUT: Stellar mass-loss rate, in Msolar/yr
    INPUT: Radius of star (in units of Rsun)
           M_star    (in units of Msun)
           Prot_star (days)
           
           It uses Eq. 7 in Johnstone and Güdel (2015, A&A, 577, A27) = J+2015
           to estimate the mass-loss rate of a low-mass main sequence star. 
           M_dot_sun_fit and Omega_sun_fit are taken as defined in J+2015. 
           The exponents a and b below are best fits obtained by J+2015, using the 
           fixed values of M_dot_sun_fit (1.4e-14 Msolar/yr) and Omega_sun_fit (2.67e-6
           rad/sec). 
           NOTE that M_dot_sun_fit is taken to be 1.4e-14 Msolar/yr in J+2015, 
           instead of the usual value of 2E-14 Msolar/yr.

    """
    #R_sun = 7e10 #Sun radius, in cm 
    #M_sun = 2e33 #Sun mass, in g
    R_sun = 1; M_sun = 1
    a = 4./3
    b = -10./3 
    M_dot_sun_fit = 1.4e-14 #Sun mass-loss rate, in g/s
    Omega_sun_fit = 2.67e-6 # Sun angular speed, in rad/s - equals to Prot_sun = 27.2367 days
    Omega_star = 2*np.pi / (Prot_star*86400)
    Mdot = M_dot_sun_fit * (R_star/R_sun)**2 * (Omega_star/Omega_sun_fit)**a * (M_star/M_sun)**b
    return Mdot

def beam_solid_angle(COMPUTE_BSA,beta_min, beta_max):
    """ Computes beam solid angle for the ECM emission cone 
        We assume an emission cone with half-opening angle theta and angular width d_theta.
        theta and d_theta are related to the speed of the electrons as 
        theta ~ d_theta ~ v/c = beta (see Melrose and Dulk YYYY and Dulk ARA&A, YYYY

        The solid angle, Omega, of a cone with half-opening angle, theta, is 
        Omega = 2*np.pi * (1. - np.cos(theta))
        Therefore, a cone sheet with inner half-opening angle, theta, and angular width d_theta
        has a solid angle Omega = 2*np.pi* ( np.cos (theta - d_theta/2) - np.cos(theta +d_theta/2)) 

        If we don't fix the beam solid angle of the emission, we can constrain other parameters, e.g., 
        eps, the efficienty factor in converting energy into Poynting flux.

        OUTPUT: Omega_min and Omega_max
                Omega_min - Float: Minimun beam solid angle, in sterradians
                Omega_max - Float: Maximum beam solid angle, in sterradians

        INPUT:  beta_min and beta_max
                beta_min - minimum speed of electrons emitting via ECM, in units of the speed of light
                beta_max - Maximum speed of electrons emitting via ECM, in units of the speed of light
    """
    if COMPUTE_BSA == False:
        Omega_min = OMEGA_MIN
        Omega_max = OMEGA_MAX
    else:
        Omega_1 = 2*np.pi * (np.cos(np.arccos(beta_min) - beta_min/2) - np.cos(np.arccos(beta_min) + beta_min/2)) 
        Omega_2 = 2*np.pi * (np.cos(np.arccos(beta_max) - beta_max/2) - np.cos(np.arccos(beta_max) + beta_max/2)) 
        Omega_min = min(Omega_1, Omega_2)
        Omega_max = max(Omega_1, Omega_2)

    return Omega_min, Omega_max

def get_S_poynt(R_obs, B_sw, v_alf, v_rel, M_A, ALPHA_SPI, geom_f):
    """
    # Returns Total Poynting flux 
    #
    # S_poynt_mks - Computed as in Saur+2013 - Eq. 55 (page 7 of 20)
    # Applies if  M_A is small (<< 1)
    # Note that for the geometric factor, we follow Turnpenney's definition, so 
    # the factor is sin^2(theta), not cos^2(theta)
    # Saur says that the power is "per hemisphere", as Zarka below
    #
    # Total Poynting flux (S_mks), in mks units [kg * m * s^(-2) * A^(-2)]
    # Poynting flux, in mks units

    # S_poynt_Z_mks - Total Poynting flux, as in Zarka 2007 (Eq. 8), but using the
    # Alfvén conductance as defined in Neubaur.
    # Eq. 8 in Zarka2007 explicityly states that it is the power "per hemisphere",
    # In this sense, this is the same as in the expresion by Saur, 
    # so there seems to be only a factor of two discrepancy, 
    #
    """
    S_poynt_mks = 2 * np.pi * (R_obs/1e2)**2 * (ALPHA_SPI*M_A)**2  \
                    * (v_alf/1e2) * (B_sw/1e4)**2 / mu_0_mks * geom_f
    S_poynt = S_poynt_mks * 1e7 # in cgs units (erg/s) 
    
    S_poynt_Z_mks = 1./ np.sqrt(1 + 1/M_A**2) *  (v_rel/1e2) \
                    * (B_sw/1e4)**2 * geom_f / mu_0_mks * np.pi*(R_obs/1e2)**2 
    S_poynt_Z     = S_poynt_Z_mks * 1e7  # in cgs units
    
    return S_poynt, S_poynt_Z

def get_S_reconnect(geom_f, R_obs, B_sw, v_rel, gamma = 1.0):
        """
        OUTPUT: 
            S_reconnect - Poynting flux (array), in cgs
            P_d_mks     - Dissipated power (array), in SI units (Watt)
            P_d         - Dissipated power (array), in cgs units (erg/s)
                          
            P_d is computed using Eq. (8) in Lanza (2009; A&A 505, 339–350).
            
        """
        P_d_mks = gamma * np.pi / mu_0_mks * (B_sw/1e4)**2 * (R_obs/1e2)**2 * (v_rel/1e2)*geom_f
        P_d     = P_d_mks * 1e7 # in cgs units 
        S_reconnect = P_d / EPSILON
        #print(type(S_reconnect))
        #print('R_obs') 
        #print(R_obs)

        return S_reconnect, P_d, P_d_mks
        
        
def get_S_stretch_and_break(R_obs, B_sw, v_rel, B_planet_arr, geom_f): 
    """
        OUTPUT: 
            S_sb - Poynting flux (array), in cgs
                          
            S_sb is computed using Eq. (5) in Strugarek 2022: MNRAS 512, 4556–4572 (2022)
            
    """  
    xi = B_sw / B_planet_arr
    f_ap = 1 - (1-(3 * xi**(1/3)/(2+xi)))**(1/2)
    S_sb = 2 * np.pi * (R_obs/1e2)**2 * (xi**(-2) * f_ap) * geom_f * (B_sw/1e4)**2 * (v_rel/1e2)/mu_0_mks
    S_sb = S_sb * 1e7 # in cgs units 
    return S_sb     
                

def get_Flux(Omega_min, Omega_max, Delta_nu_cycl, d, S_poynt):
    """ Computes the minimum and maximum expected flux densities to be received at
        Earth, in erg/s/Hz/cm2
    """
    dilution_factor_min = BETA_EFF_MIN / (Omega_max * d**2 * Delta_nu_cycl) 
    dilution_factor_max = BETA_EFF_MAX / (Omega_min * d**2 * Delta_nu_cycl)
    Flux_min = S_poynt * dilution_factor_min
    Flux_min *= 1e26 # Flux density, in mJy
    Flux_max = S_poynt * dilution_factor_max
    Flux_max *= 1e26 # Flux density, in mJy
    
    return Flux_min, Flux_max
    
    
    
def get_confinement(P_dyn_sw, P_B_sw): 
    """
     Computes the "Magnetic confinement parameter" as in ud-Doula & Owocki (2002)
    """
    P_kin_sw = P_dyn_sw/2
    eta = P_B_sw / P_kin_sw
    return eta

def identity(x):
    return x


'''
#def break_up_ratio(B_star,M_star,r_star)   

  

def get_interaction_strength(r_orb,B_star,Bplanet_field,v_alf):    
  
  kappa=15.475
  lamba=-2.082
  mu_parameter=0.5
  B_jupiter=4.170
  L_XUV=100 #times the solar value
  
  mu_0=4*np.pi*10**(-7)
  v_alf_SI=v_alf*1e-2
  
  Sigma_A=1/v_alf
  Sigma_A_SI=1/(mu_0*v_alf_SI)
  #print('Sigma_A/Sigma_A_SI :',Sigma_A/Sigma_A_SI)
  #Sigma_A=1/(mu_0*v_alf_SI)
  #Sigma_A=1/(mu_0*v_alf)
  
  Sigma_P=kappa*(r_orb/au)**lamba*(B_jupiter/Bplanet_field)*(L_XUV)**mu_parameter

  alpha_interaction_strength=Sigma_P/(Sigma_P+2*Sigma_A)
  alpha_interaction_strength_SI=Sigma_P/(Sigma_P+2*Sigma_A_SI)
  #print('alpha_interaction_strength :',alpha_interaction_strength)
  #print('alpha_interaction_strength_SI :',alpha_interaction_strength_SI)

  return Sigma_P, Sigma_A, alpha_interaction_strength
'''


def get_rss(Omega_star,M_star, R_star,v_alf_equator):
    # https://iopscience.iop.org/article/10.1088/0004-637X/814/2/99/pdf
    # convert to cgs:
    #B_star = B_star * 1e-4
    M_star = M_star * 1e-3
    R_star = R_star * 1e-2

    #Rho_star = M_star / (4 / 3 * np.pi * R_star ** 3)
    #v_A_v_esc = B_star / np.sqrt(4 * np.pi * Rho_star) / (np.sqrt(2 * G * M_star / R_star))
    G_SI=6.67*10**(-11)
    v_esc=np.sqrt(2*G_SI*M_star/R_star)
    v_esc = v_esc * 1e+2
    print('v_alf_equator: ',v_alf_equator)
    print('v_esc: ', v_esc)
    v_A_v_esc=v_alf_equator/v_esc



    
    f = Omega_star * R_star ** (3 / 2) * (np.sqrt(G * M_star)) ** (-1)

    return v_A_v_esc, f


def get_interaction_strength(r_orb, B_star, Bplanet_field, v_alf):
    kappa = 15.475
    lamba = -2.082
    mu_parameter = 0.5
    B_jupiter = 4.170
    L_XUV = 1  # times the solar value

    mu_0 = 4 * np.pi * 10 ** (-7)
    v_alf_SI = v_alf * 1e-2


    Sigma_A_SI = 1 / (mu_0 * v_alf_SI)
    Sigma_A=1/(mu_0*v_alf_SI)


    Sigma_P = kappa * (r_orb / au) ** lamba * (B_jupiter / Bplanet_field) * (L_XUV) ** mu_parameter

    alpha_interaction_strength = Sigma_P / (Sigma_P + 2 * Sigma_A)
    alpha_interaction_strength_SI = Sigma_P / (Sigma_P + 2 * Sigma_A_SI)
    # print('alpha_interaction_strength :',alpha_interaction_strength)
    # print('alpha_interaction_strength_SI :',alpha_interaction_strength_SI)

    return Sigma_P, Sigma_A, alpha_interaction_strength
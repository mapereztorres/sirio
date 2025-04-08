from SPIworkflow.constants import *

# class to handle output to file
class OutputWriter:
    def __init__(self, outfileTXT):
        self.outfileTXT = outfileTXT
    
    def write_parameters(self, T_corona, M_star_dot, mu, d, R_star, M_star, P_rot_star, B_star,
        Exoplanet, Rp, Mp, r_orb, P_orb, 
        loc_pl, M_star_dot_loc, n_base_corona, nu_plasma_corona, gyrofreq,
        Flux_r_S_min, Flux_r_S_max, rho_sw, n_sw_planet, v_sw_base, Flux_r_S_Z_min,
        Flux_r_S_Z_max, v_sw, v_rel, v_alf, M_A, B_sw,Bplanet_field, Rmp, R_planet_eff,x_larger_rms,x_smaller_rms,STUDY,
        Omega_min, Omega_max,R_planet_eff_normalized,x_superalfv):

        with open(self.outfileTXT, 'w') as f:
             f.write('# INPUT PARAMETERS:               ########\n') 
             f.write('#                                 ########\n') 
             f.write('# GENERIC (fixed for all targets) ########\n') 
             f.write('#                                 ########\n') 
             f.write('T_corona                 = {0:.0e} K\n'.format(T_corona))
             f.write('Stellar wind mass loss rate = {0:.3f} Sun mass loss rate\n'.format(M_star_dot))
             f.write('mean molecular weight    = {0:.2f}\n'.format(mu))
             #f.write('mu_0 [mks units] = {0:.0e}\n \n'.format(mu_0))
             f.write('#                                 ########\n') 
             f.write('# STAR:                           ########\n') 
             f.write('#                                 ########\n') 
             #f.write('B_star = {0:.0f} G; B_planet = {1:.0f} G\n'.format(B_star, Bp[loc_pl][0]))
             f.write('Distance to star         = {0:.2f} pc\n'.format(d/pc))
             f.write('Stellar radius           = {0:.3f} Rsun\n'.format(R_star/R_sun))
             f.write('Stellar mass             = {0:.3f} Msun\n'.format(M_star/M_sun))
             f.write('Stellar rotation period  = {0:.3f} days\n'.format(P_rot_star/day))
             f.write('Stellar magnetic field   = {0:.2f} Gauss\n'.format(B_star))
             f.write('#                                 ########\n') 
             f.write('# PLANET:                         ########\n') 
             f.write('#                                 ########\n') 
             f.write('Exoplanet name           = {0} \n'.format(Exoplanet)) 
             f.write('Planetary radius         = {0:.3f} R_Earth\n'.format(Rp/R_earth))
             f.write('Planetary mass           = {0:.3f} M_Earth\n'.format(Mp/M_earth))
             f.write('Semimajor axis of orbit  = {0:.3f} au\n'.format(r_orb/au))
             f.write('Orbital period           = {0:.3f} days\n'.format(P_orb))
             #f.write('Stellar wind particle density at the base = {0:.0e} cm-3\n'.format(n_sw_base))
             f.write('#                                 ########\n') 
             f.write('#                                 ########\n') 
             f.write('# OUTPUT PARAMETERS:              ########\n')
             f.write('#                                 ########\n') 
             f.write('n_base_corona = {0:.3e} cm^-3\n'.format(n_base_corona[M_star_dot_loc][0]))
             f.write('nu_plasma_corona = {0:.2e} MHz\n'.format(nu_plasma_corona[M_star_dot_loc][0]/1e6))
             f.write('ECMI freq (fundamental) = {0:.0f} MHz\n'.format(gyrofreq/1e6))
             f.write('Flux_ST: ({0:.2e}, {0:.2e}) mJy\n'.format(Flux_r_S_min[loc_pl][0], Flux_r_S_max[loc_pl][0]))
             #f.write('Flux_Z: ({0:.2e}, {0:.2e}) mJy\n'.format(Flux_r_S_Z_min[loc_pl][0], Flux_r_S_Z_max[loc_pl][0]))
             f.write('rho_sw at r_orb: {0:.2e} \n'.format(rho_sw[loc_pl][0]))
             f.write('n_sw_planet at r_orb: {0:.2e} \n'.format(n_sw_planet[loc_pl][0]))
             f.write('n_sw_planet at base: {0:.2e} \n'.format(n_sw_planet[0]))
             f.write('v_sw at base of the wind: {0:.2e} \n'.format(v_sw_base))
             f.write('Flux_Z: ({0:.2e}, {0:.2e}) mJy\n'.format(Flux_r_S_Z_min[loc_pl][0], Flux_r_S_Z_max[loc_pl][0]))
             f.write('v_sw at the base of the wind: {0:.2e} \n'.format(v_sw_base))
             f.write('v_sw at r_orb: {0:.2e} \n'.format(v_sw[loc_pl][0]))
             f.write('v_sw(r_orb)/v_sw_base: {0:.2e} \n'.format(v_sw[loc_pl][0]/v_sw_base))
             f.write('v_rel at r_orb: {0:.2e} \n'.format(v_rel[loc_pl][0]))
             f.write('v_alf at r_orb: {0:.2e} \n'.format(v_alf[loc_pl][0]))
             f.write('M_A at r_orb: {0:.2e} \n'.format(M_A[loc_pl][0]))
             #f.write('v_rel/v_alf at r_orb: {0:.2e} \n'.format(v_rel[loc_pl][0]/v_alf[loc_pl][0]))
             f.write('B_sw at r_orb: {0:.2e} \n'.format(B_sw[loc_pl][0]))
             f.write('Magnetic field of the planet: {0:.2e} \n'.format(Bplanet_field))
             f.write('Rmp at r_orb: {0:.2e} \n'.format(Rmp[loc_pl][0]/Rp))
             f.write('R_planet_eff at r_orb: {0:.2e} \n'.format(R_planet_eff[loc_pl][0]/Rp))
             #f.write('value of  {0} where there is clear detection:  {1:.2e}\n'.format(STUDY,x_larger_rms))
             f.write('value of '+STUDY+' where there is clear detection: '+x_larger_rms+' \n')
             f.write('value of '+STUDY+' where there is a clear NON-detection: '+x_smaller_rms+' \n')
             f.write('minimum value of solid angle of the emission: {0:.2f} \n'.format(Omega_min))
             f.write('maximum value of solid angle of the emission: {0:.2f} \n'.format(Omega_max))
             f.write('Maximum magnetopause radius in planet units: {0:.2f} \n'.format(np.max(R_planet_eff_normalized)))
             f.write('value of '+STUDY+' where the planet enters super-Alfvénic regime: '+ str(x_superalfv)+' \n')

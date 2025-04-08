import numpy as np

from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
#matplotlib.rc_file_defaults()
plt.style.use(['bmh','SPIworkflow/spi.mplstyle'])

##comparison plots for the fluxes
def construct_filename(model_type, topology, freefree):
    #planet_field_str = '[{:.3f}]'.format(Bplanet_field)
    planet_field_str='['+"{:.3f}".format(Bplanet_field)+']'
    base_path = FOLDER + '/CSV/' + STUDY + "_" + Exoplanet.replace(" ", "_")
    base_path += f"-{topology}-Bstar{B_star:.1f}G-Bplanet{planet_field_str}G"
    base_path += f"-{BETA_EFF_MIN:.1e}-{BETA_EFF_MAX:.1e}-T_corona{T_corona/1e6}MK"
    base_path += f"SPI_at_{R_ff_in/R_star}R_star"
    if freefree:
        base_path += f"_freefree_{model_type}_model.csv"
    else:
        base_path += f"_{model_type}_model.csv"
    return base_path
    
alfven_wing_parker = pd.read_csv(construct_filename('alfven_wing', 'open-parker-spiral', freefree))
reconnection_parker = pd.read_csv(construct_filename('reconnection', 'open-parker-spiral', freefree))
alfven_wing_dipole = pd.read_csv(construct_filename('alfven_wing', 'closed-dipole', freefree))
reconnection_dipole = pd.read_csv(construct_filename('reconnection', 'closed-dipole', freefree))
alfven_wing_pfss_parker = pd.read_csv(construct_filename('alfven_wing', 'pfss', freefree))
reconnection_pfss_parker = pd.read_csv(construct_filename('reconnection', 'pfss', freefree))
    
    

plt.figure(figsize=(8, 7.5))
lw=2
ax2 = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)
ax2.set_facecolor("white")    
ax2.set_ylim([YLIMLOW * 10**-2, YLIMHIGH * 10**3])     
ax2.set_yscale('log')
ax2.plot(alfven_wing_parker[STUDY], alfven_wing_parker['FLUX'], color='orange')
ax2.plot(reconnection_parker[STUDY], reconnection_parker['FLUX'], color='blue')
ax2.plot(alfven_wing_dipole[STUDY], alfven_wing_dipole['FLUX'], color='orange', linestyle='dotted')
ax2.plot(reconnection_dipole[STUDY], reconnection_dipole['FLUX'], color='blue', linestyle='dotted')
ax2.plot(alfven_wing_pfss_parker[STUDY], alfven_wing_pfss_parker['FLUX'], color='orange', linestyle='dashed')
ax2.plot(reconnection_pfss_parker[STUDY], reconnection_pfss_parker['FLUX'], color='blue', linestyle='dashed')
legend_elements = [
Line2D([0], [0], color='orange', lw=lw, label='Alfvén wing (Open Parker Spiral)'),
Line2D([0], [0], color='blue', lw=lw, label='Reconnection (Open Parker Spiral)'),
Line2D([0], [0], color='orange', linestyle='dotted', lw=lw, label='Alfvén wing (Dipole)'),
Line2D([0], [0], color='blue', linestyle='dotted', lw=lw, label='Reconnection (Dipole)'),
Line2D([0], [0], color='orange', linestyle='dashed', lw=lw, label='Alfvén wing (PFSS)'),
Line2D([0], [0], color='blue', linestyle='dashed', lw=lw, label='Reconnection (PFSS)'),
]
ax2.set_xlabel(xlabel,fontsize=20)
ax2.set_ylabel(r"Flux density [mJy]")
ax2.legend(handles=legend_elements, loc='lower left', fontsize=12, facecolor='white', edgecolor='white', framealpha=0)
if STUDY == "D_ORB":
    ax2.set_xlim(left=2)
    ax2.set_xlim(right=d_orb_max) 
    ax2.set_xscale('log')     
ax2.margins(x=0)    
secax = ax2.secondary_yaxis('right', functions=(spi.identity,spi.identity))
ax2.axvspan(x[0], R_SS, facecolor='gray', alpha=0.7)   
#plt.savefig(FOLDER+'/'+'COMPARISON_PDF'+'/'+'Flux'+'_model_comparison-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-Bstar'+"{:.1f}".format(B_star)+'G-Bplanet' + str(B_planet_arr[loc_pl]) + 'G' + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'.pdf', bbox_inches='tight')
plt.savefig(FOLDER+'/'+'COMPARISON_PDF'+'/'+'Flux'+'_model_comparison-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-Bstar'+"{:.1f}".format(B_star)+'G-Bplanet' +'['+"{:.3f}".format(Bplanet_field)+']' + 'G' + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'.pdf', bbox_inches='tight')

##ratio between reconnection and Alfven wing fluxes
plt.figure(figsize=(8, 7.5))
lw=2
ax2 = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)
ax2.set_facecolor("white")    
ax2.set_ylim(1e-1,1e3)     
ax2.set_yscale('log')
ax2.plot(alfven_wing_parker[STUDY],  reconnection_parker['FLUX']/alfven_wing_parker['FLUX'], color='blue', linestyle='dashed')
ax2.plot(alfven_wing_dipole[STUDY],  reconnection_dipole['FLUX']/alfven_wing_dipole['FLUX'], color='orange', linestyle='dotted')
ax2.plot(alfven_wing_pfss_parker[STUDY],  reconnection_pfss_parker['FLUX']/alfven_wing_pfss_parker['FLUX'], color='black')
legend_elements = [
Line2D([0], [0], color='blue', lw=lw, label='Open Parker Spiral'),
Line2D([0], [0], color='orange', linestyle='dotted', lw=lw, label='Dipole'),
Line2D([0], [0], color='black', linestyle='dashed', lw=lw, label='PFSS'),
]
ax2.set_xlabel(xlabel,fontsize=20)
ax2.set_ylabel(r"Reconnection/Alfven wing flux")
ax2.legend(handles=legend_elements, loc='lower left', fontsize=12, facecolor='white', edgecolor='white', framealpha=0)
if STUDY == "D_ORB":
    ax2.set_xlim(left=2)
    ax2.set_xlim(right=d_orb_max) 
    ax2.set_xscale('log') 
    ax2.axvspan(x[0], R_SS, facecolor='gray', alpha=0.7)
ax2.margins(x=0)    
secax = ax2.secondary_yaxis('right', functions=(spi.identity,spi.identity))
plt.savefig(FOLDER+'/'+'COMPARISON_PDF'+'/'+'Flux'+'_ratio_reconnect_Alfven-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-Bstar'+"{:.1f}".format(B_star)+'G-Bplanet' +'['+"{:.3f}".format(Bplanet_field)+']' + 'G' + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'.pdf', bbox_inches='tight')




##comparison plots for the stellar wind magnetic field
FOLDER + '/CSV/' +"diagnostic-" + STUDY + "_" + str(Exoplanet.replace(" ", "_")) +  geometry + diagnostic_string +'_B_sw.csv'
diagnostic_string = "{:.1f}".format(B_star) + "G" + "-Bplanet" + str(B_planet_arr[loc_pl]) + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'


def construct_diagnostic_filename(topology):
    #planet_field_str = '[{:.3f}]'.format(Bplanet_field)
    planet_field_str='['+"{:.3f}".format(Bplanet_field)+']'
    base_path = FOLDER + '/CSV/diagnostic-' + STUDY + "_" + Exoplanet.replace(" ", "_")
    base_path += f"-{topology}-Bstar{B_star:.1f}G-Bplanet{planet_field_str}G"
    base_path += f"-{BETA_EFF_MIN:.1e}-{BETA_EFF_MAX:.1e}-T_corona{T_corona/1e6}MK"
    base_path += f"SPI_at_{R_ff_in/R_star}R_star_B_sw.csv"
    return base_path	
    
bsw_parker = pd.read_csv(construct_diagnostic_filename('open-parker-spiral'))
bsw_dipole = pd.read_csv(construct_diagnostic_filename('closed-dipole'))
bsw_pffs_parker = pd.read_csv(construct_diagnostic_filename('pfss'))

plt.figure(figsize=(8, 7.5))
lw=2
ax2 = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)
ax2.set_facecolor("white")  


    
ax2.set_yscale('log')
ax2.plot(bsw_parker[STUDY], bsw_parker['Bsw'], color='blue', linestyle='dashed')
ax2.plot(bsw_dipole[STUDY], bsw_dipole['Bsw'], color='orange', linestyle='dotted')
#ax2.plot(bsw_hybrid[STUDY], bsw_hybrid['Bsw'], color='blue', linestyle='dashed')
ax2.plot(bsw_pffs_parker[STUDY], bsw_pffs_parker['Bsw'], color='black')
legend_elements = [
Line2D([0], [0], color='blue', linestyle='dashed', lw=lw, label='Open Parker Spiral'),
Line2D([0], [0], color='orange', linestyle='dotted', lw=lw, label='Dipole'),
#Line2D([0], [0], color='blue', linestyle='dashed', lw=2, label='Hybrid'),
Line2D([0], [0], color='black', lw=lw, label='PFSS'),
]

ax2.legend(handles=legend_elements, loc='lower left', fontsize=12, facecolor='white', edgecolor='white', framealpha=0)
if STUDY == "D_ORB":
    ax2.set_xlim(left=2)   
    ax2.set_xlim(right=d_orb_max)
    ax2.set_xscale('log')      
ax2.set_xlabel(xlabel,fontsize=20)
ax2.set_ylabel(r"$B_{\rm sw}$ $[G]$")
secax = ax2.secondary_yaxis('right', functions=(spi.identity,spi.identity))
ax2.margins(x=0)
if STUDY == 'MDOT':
    ax2.set_xscale('log')

ax2.axvspan(x[0], R_SS, facecolor='gray', alpha=0.7)       
#plt.savefig(FOLDER + '/' +'COMPARISON_PDF'+'/'+'B_sw_'+'model_comparison-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-hybrid-Bstar'+"{:.1f}".format(B_star) + "G" + "-Bplanet" + str(B_planet_arr[loc_pl]) + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'.pdf', bbox_inches='tight')
plt.savefig(FOLDER + '/' +'COMPARISON_PDF'+'/'+'B_sw_'+'model_comparison-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-hybrid-Bstar'+"{:.1f}".format(B_star) + "G" + "-Bplanet" + '['+"{:.3f}".format(Bplanet_field)+']'+ "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'.pdf', bbox_inches='tight')




##comparison plots for M_A

def construct_diagnostic_ma_filename(topology):
    #planet_field_str = '[{:.3f}]'.format(Bplanet_field)
    planet_field_str='['+"{:.3f}".format(Bplanet_field)+']'
    base_path = FOLDER + '/CSV/diagnostic-' + STUDY + "_" + Exoplanet.replace(" ", "_")
    base_path += f"-{topology}-Bstar{B_star:.1f}G-Bplanet{planet_field_str}G"
    base_path += f"-{BETA_EFF_MIN:.1e}-{BETA_EFF_MAX:.1e}-T_corona{T_corona/1e6}MK"
    base_path += f"SPI_at_{R_ff_in/R_star}R_star_M_A.csv"
    return base_path
    
bsw_parker = pd.read_csv(construct_diagnostic_ma_filename('open-parker-spiral'))
bsw_dipole = pd.read_csv(construct_diagnostic_ma_filename('closed-dipole'))
bsw_pffs_parker = pd.read_csv(construct_diagnostic_ma_filename('pfss'))
    
    
plt.figure(figsize=(8, 7.5))
lw=2
ax2 = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)
ax2.set_facecolor("white")  
if STUDY == 'MDOT':
    ax2.set_xscale('log')
ax2.set_yscale('log')

ax2.plot(bsw_parker[STUDY], bsw_parker['M_A'], color='blue', linestyle='dashed')
ax2.plot(bsw_dipole[STUDY],  bsw_dipole['M_A'], color='orange', linestyle='dotted')
#ax2.plot(bsw_hybrid[STUDY], bsw_hybrid['Bsw'], color='blue', linestyle='dashed')
ax2.plot(bsw_pffs_parker[STUDY], bsw_pffs_parker['M_A'], color='black')
legend_elements = [
Line2D([0], [0], color='blue', linestyle='dashed', lw=lw, label='Open Parker Spiral'),
Line2D([0], [0], color='orange', linestyle='dotted', lw=lw, label='Dipole'),
#Line2D([0], [0], color='blue', linestyle='dashed', lw=2, label='Hybrid'),
Line2D([0], [0], color='black', lw=lw, label='PFSS'),
]




ax2.legend(handles=legend_elements, loc='upper left', fontsize=12, facecolor='white', edgecolor='white', framealpha=0)

if STUDY == "D_ORB":
    ax2.set_xlim(left=2)     
    ax2.set_xlim(right=d_orb_max)  
    ax2.set_xscale('log')   
ax2.set_xlabel(xlabel,fontsize=20)
ax2.set_ylabel(r"$M_A$")
secax = ax2.secondary_yaxis('right', functions=(spi.identity,spi.identity))

if Exoplanet=='YZCet b Model A' or Exoplanet=='YZCet b Model B':
    ax2.set_ylim([1e-3, 1e2]) 
    
ax2.axvspan(x[0], R_SS, facecolor='gray', alpha=0.7)   
ax2.margins(x=0)
#fig.tight_layout()
#plt.savefig(FOLDER + '/'+'COMPARISON_PDF'+'/'+'M_A_'+'model_comparison-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-hybrid-Bstar'+"{:.1f}".format(B_star) + "G" + "-Bplanet" + str(B_planet_arr[loc_pl]) + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'.pdf', bbox_inches='tight')
plt.savefig(FOLDER + '/'+'COMPARISON_PDF'+'/'+'M_A_'+'model_comparison-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-hybrid-Bstar'+"{:.1f}".format(B_star) + "G" + "-Bplanet" + '['+"{:.3f}".format(Bplanet_field)+']'+ "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'.pdf', bbox_inches='tight')







#print(FOLDER + '/' +'diagnostic-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-pfss-Bstar'+"{:.1f}".format(B_star) + "G" + "-Bplanet" + str(B_planet_arr[loc_pl]) + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'_B_sw.csv')
'''

'''

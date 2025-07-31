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
    #ax2.axvspan(x[0], R_SS, facecolor='gray', alpha=0.7)
    ax2.axvline(R_SS, color='grey', alpha=0.9, linestyle='-', lw=3)
    ax2.text(R_SS*0.8, 1, rf'$R_{{SS}}$ = {R_SS}', fontsize=11, alpha=1,
             bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.2'))

ax2.margins(x=0)    
secax = ax2.secondary_yaxis('right', functions=(spi.identity,spi.identity))



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
Line2D([0], [0], color='blue', lw=lw, label='Parker Spiral'),
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
    #ax2.axvspan(x[0], R_SS, facecolor='gray', alpha=0.7)
    ax2.axvline(R_SS, color='grey', alpha=0.9, linestyle='-', lw=3)
    ax2.text(R_SS*0.8, 1, rf'$R_{{SS}}$ = {R_SS}', fontsize=11, alpha=1,
             bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.2'))
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
Line2D([0], [0], color='blue', linestyle='dashed', lw=lw, label='Parker Spiral'),
Line2D([0], [0], color='orange', linestyle='dotted', lw=lw, label='Dipole'),
#Line2D([0], [0], color='blue', linestyle='dashed', lw=2, label='Hybrid'),
Line2D([0], [0], color='black', lw=lw, label='PFSS'),
]

ax2.legend(handles=legend_elements, loc='lower left', fontsize=12, facecolor='white', edgecolor='white', framealpha=0)
if STUDY == "D_ORB":
    ax2.set_xlim(left=2)   
    ax2.set_xlim(right=d_orb_max)
    ax2.set_xscale('log')
    #ax2.axvspan(x[0], R_SS, facecolor='gray', alpha=0.7)
    ax2.axvline(R_SS, color='grey', alpha=0.9, linestyle='-', lw=3)
    ax2.text(R_SS*0.8, 1, rf'$R_{{SS}}$ = {R_SS}', fontsize=11, alpha=1,
             bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.2'))

ax2.set_xlabel(xlabel,fontsize=20)
ax2.set_ylabel(r"$B_{\rm sw}$ $[G]$")
secax = ax2.secondary_yaxis('right', functions=(spi.identity,spi.identity))
ax2.margins(x=0)
if STUDY == 'MDOT':
    ax2.set_xscale('log')

#ax2.axvspan(x[0], R_SS, facecolor='gray', alpha=0.7)
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
Line2D([0], [0], color='blue', linestyle='dashed', lw=lw, label='Parker Spiral'),
Line2D([0], [0], color='orange', linestyle='dotted', lw=lw, label='Dipole'),
#Line2D([0], [0], color='blue', linestyle='dashed', lw=2, label='Hybrid'),
Line2D([0], [0], color='black', lw=lw, label='PFSS'),
Line2D([], [], color='none', label=f'$B_{{*}}$ = {B_star} G')
]




ax2.legend(handles=legend_elements, loc='lower right', fontsize=16, facecolor='white', edgecolor='white', framealpha=0)

if STUDY == "D_ORB":
    ax2.set_xlim(left=2)     
    ax2.set_xlim(right=d_orb_max)  
    ax2.set_xscale('log')   
    
if STUDY == "M_DOT":
    ax2.set_xscale('log')       
ax2.set_xlabel(xlabel,fontsize=20)
ax2.set_ylabel(r"$M_A$")
secax = ax2.secondary_yaxis('right', functions=(spi.identity,spi.identity))

if Exoplanet=='YZCet b Model A' or Exoplanet=='YZCet b Model B':
    ax2.set_ylim([1e-3, 1e2]) 

'''
#if Exoplanet == 'Proxima b Turn-regular' or Exoplanet == 'Proxima b':
ax2.axvline(x = 0.04856*au/R_star, ls='-.', color='k', lw=1.5)  
ax2.axvline(x = 0.02885*au/R_star, ls='-.', color='k', lw=1.5)  
ax2.text(0.04856*au/R_star-1,2e4,'b',ha='center',fontsize=11) 
ax2.text(0.02885*au/R_star-1,2e4,'d',ha='center',fontsize=11)
ax2.set_xlim([2,3e3])
'''

if Exoplanet == 'Proxima b Turn-regular' or Exoplanet == 'Proxima b Turnpenney' or Exoplanet == 'Proxima b kavanagh' or Exoplanet=='Proxima b Reville':
    if STUDY == 'D_ORB':
        #ax2.axvline(x = 0.04856*au/R_star, ls='-.', color='k', lw=1.5)  
        ax2.axvline(x = 0.02885*au/R_star, ls='-.', color='k', lw=1.5)  
        #ax2.text(0.04856*au/R_star-10,1e4,'b',ha='center',fontsize=15) 
        #ax2.text(0.02885*au/R_star-7,1e4,'d',ha='center',fontsize=15)
        #ax2.text(0.04856*au/R_star-10,9e1,'b',ha='center',fontsize=15,bbox=dict(facecolor='white', edgecolor='none', boxstyle='square,pad=0.1')) 
        #ax2.text(0.02885*au/R_star-7,9e1,'d',ha='center',fontsize=15,bbox=dict(facecolor='white', edgecolor='none', boxstyle='square,pad=0.1'))
        #ax2.text(0.04856*au/R_star*0.9,6e1,'b',ha='center',fontsize=15,bbox=dict(facecolor='white', edgecolor='none', boxstyle='square,pad=0.1')) 
        ax2.text(xnom*0.9,6e1,'b',ha='center',fontsize=15,bbox=dict(facecolor='white', edgecolor='none', boxstyle='square,pad=0.1')) 
        ax2.text(0.02885*au/R_star*0.9,6e1,'d',ha='center',fontsize=15,bbox=dict(facecolor='white', edgecolor='none', boxstyle='square,pad=0.1'))
        #ax2.set_xlim([2,3e3])
        ax2.set_xlim([2,2e2])
        ax2.legend(handles=legend_elements, loc='upper left', fontsize=16, facecolor='white', edgecolor='white', framealpha=1)
        #ax2.axvline(x = xnom, ls='--', color='k', lw=2)
'''
if Exoplanet=='Trappist-1 b' or Exoplanet=='Trappist-1 b Reville':
    if STUDY == 'D_ORB':
        ax2.set_xlim([1,31])
        ax2.legend(handles=legend_elements, loc='upper left', fontsize=16, facecolor='white', edgecolor='white', framealpha=1)
        #ax2.set_xscale('linear')
    if STUDY == 'M_DOT':   
        #ax2.set_xlabel(xlabel,fontsize=20) 
        ax2.set_xlabel(r"Mass Loss rate [$\dot{M}_\odot$]",fontsize=20)
        ax2.set_xscale('log')
'''        
if Exoplanet == 'GJ1151 hypothetical 2' or Exoplanet == 'GJ1151 hypothetical 1':
    if STUDY == 'D_ORB':
        #x_pl1=12.927
        x_pl1=11.652574895873382
        print('orbital separation ranging from '+str(x_pl1)+' to '+str(xnom))
        ax2.axvline(x = x_pl1, ls='--', color='k', lw=2)
        #ax2.axvline(x = xnom, ls='--', color='k', lw=2)
        ax2.set_ylim(1e-3,1e1)
        ax2.set_xlim([2,50])
        ax2.set_xlabel(xlabel,fontsize=20)
        yticks2 = ax2.get_yticks()
        ax2.set_yticks(yticks2[2:-2]) 
        secax2 = ax2.secondary_yaxis('right', functions=(spi.identity,spi.identity))
        secax2.set_yticks(yticks2[2:-2]) 
        ax2.set_xscale('linear')
        ax2.legend(handles=legend_elements, loc='lower right', fontsize=16, facecolor='white', edgecolor='white', framealpha=1)
        ax2.set_ylim(1.01e-3,9.9)


planet_height_label=7e1
bbox=dict(facecolor='white', edgecolor='none', boxstyle='square,pad=0.1')
if 'Trappist' in Exoplanet:
    if STUDY == 'D_ORB':
        ax2.set_xlim([1,150])
        
        #ax2.set_xscale('linear')
        #c      
        ax2.axvline(x = 0.01521*au/R_star, ls='-.', color='k', lw=1.5)     
        #d
        ax2.axvline(x = 0.02144*au/R_star, ls='-.', color='k', lw=1.5)   
        #e       
        ax2.axvline(x = 0.02817*au/R_star, ls='-.', color='k', lw=1.5)        
        #f       
        ax2.axvline(x = 0.0371*au/R_star, ls='-.', color='k', lw=1.5)   
        #g
        ax2.axvline(x = 0.0451*au/R_star, ls='-.', color='k', lw=1.5)      
        #h
        ax2.axvline(x = 0.063*au/R_star, ls='-.', color='k', lw=1.5)
        ax2.text(0.01111*au/R_star,planet_height_label,'b',ha='center',fontsize=13,bbox=bbox)
        ax2.text(0.01521*au/R_star,planet_height_label,'c',ha='center',fontsize=13,bbox=bbox)    
        ax2.text(0.02144*au/R_star,planet_height_label,'d',ha='center',fontsize=13,bbox=bbox)
        ax2.text(0.02817*au/R_star,planet_height_label,'e',ha='center',fontsize=13,bbox=bbox)
        ax2.text(0.0371*au/R_star,planet_height_label,'f',ha='center',fontsize=13,bbox=bbox)
        ax2.text(0.0451*au/R_star,planet_height_label,'g',ha='center',fontsize=13,bbox=bbox)
        ax2.text(0.063*au/R_star,planet_height_label,'h',ha='center',fontsize=13,bbox=bbox) 
        ax2.text(0.063*au/R_star,planet_height_label,'h',ha='center',fontsize=13,bbox=bbox)
        #ax2.text(xnom,2e-3,'$\dot{M}_{*}$= '+str(M_star_dot_arr[0])+' $\dot{M}_{\odot}$',ha='center',fontsize=16,bbox=bbox)
        legend_elements = [
Line2D([0], [0], color='blue', linestyle='dashed', lw=lw, label='Parker Spiral'),
Line2D([0], [0], color='orange', linestyle='dotted', lw=lw, label='Dipole'),
#Line2D([0], [0], color='blue', linestyle='dashed', lw=2, label='Hybrid'),
Line2D([0], [0], color='black', lw=lw, label='PFSS'),
Line2D([], [], color='none', label=f'$B_{{*}}$ = {B_star} G'),
Line2D([], [], color='none', label=rf'$\dot{{M}}_{{*}} = {M_star_dot_arr[0]}\ \dot{{M}}_{{\odot}}$')
        ]
        ax2.legend(handles=legend_elements, loc='upper left', fontsize=16, facecolor='white', edgecolor='white', framealpha=1)
    if STUDY == 'M_DOT':   
        #ax2.set_xlabel(xlabel,fontsize=20) 
        ax2.set_xlabel(r"Mass Loss rate [$\dot{M}_\odot$]",fontsize=20)
        ax2.set_xscale('log')	


ax2.axvline(x = xnom, ls='--', color='k', lw=2)    
ax2.axhline(y = 1, ls='-.', color='grey', lw=2)   
ax2.set_ylim([1e-3,1e2])    
if STUDY == "D_ORB":    
    #ax2.axvspan(x[0], R_SS, facecolor='gray', alpha=0.7)
    ax2.axvline(R_SS, color='grey', alpha=0.9, linestyle='-', lw=3)
    #ax2.text(3.5, 2e-3, rf'$R_{{SS}}$ = {R_SS} $R_{{*}}$', fontsize=14, alpha=1,
    ax2.text(R_SS*0.8, 2e-3, rf'$R_{{SS}}$', fontsize=20, alpha=1,
             bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.2'))
if STUDY == 'M_DOT':   
    ax2.set_xscale('log')    
    ax2.set_xlim([1.01e-1,9.9e1])
ax2.margins(x=0)
#fig.tight_layout()
#plt.savefig(FOLDER + '/'+'COMPARISON_PDF'+'/'+'M_A_'+'model_comparison-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-hybrid-Bstar'+"{:.1f}".format(B_star) + "G" + "-Bplanet" + str(B_planet_arr[loc_pl]) + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'.pdf', bbox_inches='tight')
plt.savefig(FOLDER + '/'+'COMPARISON_PDF'+'/'+'M_A_'+'model_comparison-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-hybrid-Bstar'+"{:.1f}".format(B_star) + "G" + "-Bplanet" + '['+"{:.3f}".format(Bplanet_field)+']'+ "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'.pdf', bbox_inches='tight')







#print(FOLDER + '/' +'diagnostic-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-pfss-Bstar'+"{:.1f}".format(B_star) + "G" + "-Bplanet" + str(B_planet_arr[loc_pl]) + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'_B_sw.csv')
'''

'''

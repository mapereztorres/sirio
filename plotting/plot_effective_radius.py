import numpy as np

from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.patches as mpatches
#matplotlib.rc_file_defaults()
plt.style.use(['bmh','SPIworkflow/spi.mplstyle'])

plt.figure(figsize=(8,7.5))
ax = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
ax.plot(x, R_obs_normalized, color='k')
ax.plot(x, Rmp/Rp, color='r')

#ax.plot(x, Reff_lanza*np.ones(len(x)), color='blue')


ax.set_xlabel(xlabel,fontsize=20)
ax.set_ylabel(r"$R(R_{pl})$")
ax.set_facecolor("white")

ax.set_ylim([0,2])

ax.set_xlim(left=0)

ax.axvline(x = xnom, ls='--', color='k', lw=2)

if Bfield_geom_arr[ind] == 'pfss':
    ax.axvline(R_SS, color='grey', alpha=0.9, linestyle='-', lw=3)
    #ax.text(7.5, 0.1, rf'$R_{{SS}}$ = {R_SS} $R_{{*}}$', fontsize=14, alpha=1,bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.2'),rotation=90)
    ax.text(R_SS*0.8, 0.1, rf'$R_{{SS}}$', fontsize=18, alpha=1,bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.2'))
black_patch = mpatches.Patch(color='black', label='$R_{eff}$')
red_patch = mpatches.Patch(color='red', label='$R_{mp}$')
#blue_patch = mpatches.Patch(color='blue', label='$R_{mp\_Lanza}$')

if STUDY == 'M_DOT':
    ax.set_xscale('log')

#ax.legend(handles=[black_patch,red_patch,blue_patch],loc='upper left',fontsize=20,facecolor='white',edgecolor='white', framealpha=0)
ax.legend(handles=[black_patch,red_patch],loc='upper left',fontsize=20,facecolor='white',edgecolor='white', framealpha=0)
ax.set_xlim(right=x[-1]) #the limit on the right is the last element of the x array
secax = ax.secondary_yaxis('right', functions=(spi.identity,spi.identity))     
plt.savefig(FOLDER + '/' 'R_EFF_PDF'+'/'+ str(Exoplanet.replace(" ", "_"))
        +'-effective_radius_variation-'+STUDY+geometry+ "-Bplanet" +'['+"{:.3f}".format(Bplanet_field)+']' + "G" +'-'+'T_corona'+str(T_corona/1e6)+'MK'+'-'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'.pdf', bbox_inches='tight')
   
   


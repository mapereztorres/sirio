import numpy as np

from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.patches as mpatches
#matplotlib.rc_file_defaults()
plt.style.use(['bmh','SPIworkflow/spi.mplstyle'])

plt.figure(figsize=(8,7.5))
ax = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)

lw=3

ax.plot(x, R_obs_normalized, color='k',lw=lw+2)
ax.plot(x, Rmp/Rp, color='r',lw=lw+2)


ax.set_xlabel(xlabel,fontsize=30)

ax.set_facecolor("white")

ax.set_ylim([0,2])

ax.set_xlim(left=1)

ax.axvline(x = xnom, ls='--', color='k',lw=lw+2)

ax.set_ylabel(r"$R(R_{pl})$",fontsize=40)       

if Bfield_geom_arr[ind] == 'pfss':
    ax.axvline(R_SS, color='grey', alpha=0.9, linestyle='-',lw=lw+2)

    ax.text(R_SS*0.8, 0.1, rf'$R_{{SS}}$', fontsize=22, alpha=1,bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.2'))
black_patch = mpatches.Patch(color='black', label='$R_{eff}$')
red_patch = mpatches.Patch(color='red', label='$R_{mp}$')

if STUDY == 'M_DOT':
    ax.set_xscale('log')

ax.legend(handles=[black_patch,red_patch],loc='upper left',fontsize=25,facecolor='white',edgecolor='white', framealpha=0)
ax.set_xlim(right=x[-1]) #the limit on the right is the last element of the x array
secax = ax.secondary_yaxis('right', functions=(spi.identity,spi.identity))     


title_str = geometry.replace("-", " ")  
title_str = title_str.replace("Bstar", "") 
title_str = title_str.strip()  
title_str = title_str.title()  



if title_str=="Pfss":
    title_str="PFSS"

ax.set_title(title_str, fontsize=40,pad=20)

for ax in [ax,secax]:
    ax.tick_params(axis='both', which='major', labelsize=25, width=2, length=8)
    ax.tick_params(axis='both', which='minor', labelsize=25, width=1.5, length=5)  
    ax.xaxis.label.set_size(30)
    ax.yaxis.label.set_size(30)

    for spine in ax.spines.values():
        spine.set_linewidth(4)


ax.yaxis.label.set_size(50)   
plt.savefig(FOLDER + '/' 'R_EFF_PDF'+'/'+ str(Exoplanet.replace(" ", "_"))
        +'-effective_radius_variation-'+STUDY+geometry+ "-Bplanet" +'['+"{:.3f}".format(Bplanet_field)+']' + "G" +'-'+'T_corona'+str(T_corona/1e6)+'MK'+'-'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'.pdf', bbox_inches='tight')
   
   


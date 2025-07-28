import numpy as np

from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.patches as mpatches
#matplotlib.rc_file_defaults()
plt.style.use(['bmh','SPIworkflow/spi.mplstyle'])

## All diagnostic plots together 

fig, (ax1, ax2,ax3,ax4) = plt.subplots(4, 1, sharex=True)
x = d_orb / R_star
fig.subplots_adjust(hspace=0)
ax1.plot(x, v_orb/1e5*np.ones(len(x)), color='r', linestyle='dotted')
ax1.plot(x, v_alf/1e5*np.ones(len(x)), color='g', linestyle='dashdot')
ax1.plot(x, v_sw/1e5*np.ones(len(x)), color='b', linestyle='dashed')
ax1.plot(x, v_rel/1e5*np.ones(len(x)), color='k', linestyle='solid')
ax1.plot(x, v_sound/1e5*np.ones(len(x)), color='orange', linestyle=(0,(1,3.5)))
import sys
np.set_printoptions(threshold=sys.maxsize)
#print('v_rel/1e5 :')
#print(v_rel/1e5)
#print('x:')
#print(x)
#print(d_orb)
'''
red_patch = mpatches.Patch(color='red', label=r'v$_{\rm orb}$')
black_patch = mpatches.Patch(color='black', label=r'v$_{\rm rel}$')
green_patch = mpatches.Patch(color='green', label=r'v$_{\rm alf}$')
blue_patch = mpatches.Patch(color='blue', label=r'v$_{\rm sw}$')
orange_patch = mpatches.Patch(color='orange', label=r'v$_{\rm sound}$')      
#ax1.set_ylim(100,2000)
ax1.legend(handles=[blue_patch,red_patch,black_patch,orange_patch,green_patch],loc='upper left',fontsize=20,facecolor='white',edgecolor='white', framealpha=0)
'''
#ax1.axhline(y = v_sw_terminal/1e5)
legend_elements = [
Line2D([0], [0], color='r', linestyle='dotted', label=r'v$_{\rm orb}$'),
Line2D([0], [0], color='green', linestyle='dashdot', label=r'v$_{\rm alf}$'),
Line2D([0], [0], color='blue', linestyle='dashed',  label=r'v$_{\rm sw}$'),
Line2D([0], [0], color='black', linestyle='solid', label=r'v$_{\rm rel}$'),
Line2D([0], [0], color='orange', linestyle=(0,(1,3.5)), label=r'v$_{\rm sound}$'),
]
ax1.legend(handles=legend_elements, loc='upper right', fontsize=12, facecolor='white', edgecolor='white', framealpha=0)



ax2.plot(x, np.abs(B_r)*np.ones(len(x))*1e5, color='r', linestyle='dotted')
ax2.plot(x, B_phi*np.ones(len(x))*1e5, color='g', linestyle='dashdot')
ax2.plot(x, B_sw*np.ones(len(x))*1e5, color='b', linestyle='solid')
ax2.plot(x, B_sw*np.ones(len(x))*np.sqrt(geom_f)*1e5, color='k', linestyle='dashed')

'''
red_patch = mpatches.Patch(color='red', label=r'B$_{\rm r}$')
green_patch = mpatches.Patch(color='green', label=r'B$_{\rm phi}$')
blue_patch = mpatches.Patch(color='blue', label=r'B$_{\rm tot}$')
black_patch = mpatches.Patch(color='black', label=r'B$_{\rm perp}$')
ax2.legend(handles=[red_patch,green_patch,blue_patch,black_patch],loc='upper right',fontsize=20,facecolor='white',edgecolor='white', framealpha=0)
'''
legend_elements = [
Line2D([0], [0], color='red', linestyle='dotted', label=r'B$_{\rm r}$'),
Line2D([0], [0], color='green', linestyle='dashdot', label=r'B$_{\rm phi}$'),
Line2D([0], [0], color='blue', linestyle='solid',  label=r'B$_{\rm tot}$'),
Line2D([0], [0], color='black', linestyle='dashed', label=r'B$_{\rm perp}$'),
]
ax2.legend(handles=legend_elements, loc='upper right', fontsize=12, facecolor='white', edgecolor='white', framealpha=0)



ax3.plot(x, M_A*np.ones(len(x)), color='k', lw=lw)
#ax3.plot(x, geom_f*np.ones(len(x)), color='red', lw=lw)

#ax3.plot(x, eta*np.ones(len(x)), color='k', lw=lw)
#ax3.axyline(y = xnom, ls='--', color='k', lw=2)
#ax3.axvline(x = R_alfven)

ax4.plot(x, P_B_sw*np.ones(len(x))*1e8, color='b', linestyle='dashdot')
ax4.plot(x, P_dyn_sw*np.ones(len(x))*1e8, color='r', linestyle='solid')
ax4.plot(x, P_th_sw*np.ones(len(x))*1e8, color='g', linestyle=(0,(1,3.5)))
# LPM - CHECK IF USEFUL
#ax4.plot(x, P_sw*np.ones(len(x)), color='k', ls=(0, (1, 2)))
#ax4.plot(x, P_B_planet*np.ones(len(x)), color='magenta', ls=(0, (1, 2)))
'''
black_patch = mpatches.Patch(color='black', label=r'P$_{\rm sw}$') #r'Primary T$_{\rm eff}$')
red_patch = mpatches.Patch(color='red', label=r'P$_{\rm dyn_{\rm sw}}$')
green_patch = mpatches.Patch(color='green', label=r'P$_{\rm th_{\rm sw}}$')
blue_patch = mpatches.Patch(color='blue', label=r'P$_{\rm B_{\rm sw}}$')
magenta_patch = mpatches.Patch(color='magenta', label=r'P$_{\rm B_planet}$')             
ax4.legend(handles=[blue_patch,red_patch,green_patch],loc='upper right',fontsize=20,facecolor='white',edgecolor='white', framealpha=0)
'''
legend_elements = [
Line2D([0], [0],  color='r', linestyle='solid', label=r'P$_{\rm dyn_{\rm sw}}$'),
Line2D([0], [0],color='g', linestyle=(0,(1,3.5)), label=r'P$_{\rm th_{\rm sw}}$'),
Line2D([0], [0], color='blue', linestyle='dashdot', label=r'P$_{\rm B_{\rm sw}}$'),
]
ax4.legend(handles=legend_elements, loc='upper right', fontsize=12, facecolor='white', edgecolor='white', framealpha=0)




ax1.set_ylabel(r"v $[\rm km$ $s^{-1}] $")
ax2.set_ylabel(r"$B_{\rm sw}$ $[nT]$")
ax3.set_ylabel(r"$M_A$")
ax4.set_ylabel(r"P $[nPa] $")
        
ax1.set_facecolor("white")
ax2.set_facecolor("white")
ax3.set_facecolor("white")
ax4.set_facecolor("white")

ax1.axvline(x = xnom, ls='--', color='k', lw=2)
ax2.axvline(x = xnom, ls='--', color='k', lw=2)
ax3.axvline(x = xnom, ls='--', color='k', lw=2)
ax4.axvline(x = xnom, ls='--', color='k', lw=2)
     
     
  



#print('geometry', geometry)
#if STUDY == "D_ORB" and Bfield_geom_arr[ind] == 'pfss':
#    ax1.axvspan(x[0], R_SS, facecolor='grey', alpha=0.6)
#    ax2.axvspan(x[0], R_SS, facecolor='grey', alpha=0.6)
#    ax3.axvspan(x[0], R_SS, facecolor='grey', alpha=0.6)
#    ax4.axvspan(x[0], R_SS, facecolor='grey', alpha=0.6)
    
if STUDY == "D_ORB" and Bfield_geom_arr[ind] == 'pfss':
    for ax in [ax1, ax2, ax3, ax4]:
        ax.axvline(R_SS, color='grey', alpha=0.9, linestyle='-', lw=3)
    ax4.text(3.5, 1, rf'$R_{{SS}}$ = {R_SS}', fontsize=11, alpha=1,
         bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.2'))

if (M_A > 1).any():
    ax3.axhline(y = 1, ls='-.', color='grey', lw=2)   
ax4.set_xlabel(xlabel,fontsize=20)

if STUDY == "D_ORB":
    ax4.set_xlim([2,x[-1]])

fig.set_figwidth(8)
fig.set_figheight(20)
#ax1.set_xscale('log')
ax1.set_yscale('log')  
ax2.set_yscale('log')  
ax3.set_yscale('log')  
ax4.set_yscale('log')  

#ax2.set_ylim([1e-5,1e5])



print('Plots like in Turnpenny2018')

'''
ax1.set_ylim(1e1,1e5)
#ax2.set_ylim(1e0,1e10)  
#ax3.set_ylim(1e-3,1e1)
ax3.set_ylim([1e-7,1e1])
ax4.set_ylim(1e-2,1e14)
'''
ax1.set_ylim(1e1,1e4)
#ax2.set_ylim(1e0,1e10)  
ax3.set_ylim(1e-3,1e1)
ax4.set_ylim(1,1e14)
ax4.set_xscale('log')

if Exoplanet=='Trappist-1 b':  
    ax2.set_ylim(1e0,1e10)  
 
if Exoplanet=='Proxima b Turnpenney':  
    ax2.set_ylim(10**(-3),1e8)
 
    
yticks1 = ax1.get_yticks()
ax1.set_yticks(yticks1[2:])    
  
yticks2 = ax2.get_yticks()
ax2.set_yticks(yticks2[1:]) 
    
yticks3 = ax3.get_yticks()
ax3.set_yticks(yticks3[2:-2]) 

yticks4 = ax4.get_yticks()
ax4.set_yticks(yticks4[:-2]) 

      
secax1 = ax1.secondary_yaxis('right', functions=(spi.identity,spi.identity))
secax2 = ax2.secondary_yaxis('right', functions=(spi.identity,spi.identity))
secax3 = ax3.secondary_yaxis('right', functions=(spi.identity,spi.identity))
secax4 = ax4.secondary_yaxis('right', functions=(spi.identity,spi.identity))      
    
secax1.set_yticks(yticks1[2:])       
secax2.set_yticks(yticks2[1:])
secax3.set_yticks(yticks3[2:-2]) 
secax4.set_yticks(yticks4[:-2])
    
if Exoplanet=='Proxima b Turnpenney':  
    ax2.set_yticks(yticks2[2:-1])    
    secax2.set_yticks(yticks2[2:-1])
    ax3.set_yticks(yticks3[2:-1])    
    secax3.set_yticks(yticks3[2:-1])                 


planet_height_label=5e4
bbox=dict(facecolor='white', edgecolor='none', boxstyle='square,pad=0.1')
if Exoplanet=='Trappist-1 b' or Exoplanet=='Trappist-1 b Reville':
    if STUDY == 'D_ORB': 
        #ax2.set_xlim([1,3e3])
        #c      
        ax1.axvline(x = 0.01521*au/R_star, ls='-.', color='k', lw=1.5)
        ax2.axvline(x = 0.01521*au/R_star, ls='-.', color='k', lw=1.5)
        ax3.axvline(x = 0.01521*au/R_star, ls='-.', color='k', lw=1.5)
        ax4.axvline(x = 0.01521*au/R_star, ls='-.', color='k', lw=1.5)     
        #d
        ax1.axvline(x = 0.02144*au/R_star, ls='-.', color='k', lw=1.5)
        ax2.axvline(x = 0.02144*au/R_star, ls='-.', color='k', lw=1.5)
        ax3.axvline(x = 0.02144*au/R_star, ls='-.', color='k', lw=1.5)
        ax4.axvline(x = 0.02144*au/R_star, ls='-.', color='k', lw=1.5)      
        #e       
        ax1.axvline(x = 0.02817*au/R_star, ls='-.', color='k', lw=1.5)
        ax2.axvline(x = 0.02817*au/R_star, ls='-.', color='k', lw=1.5)
        ax3.axvline(x = 0.02817*au/R_star, ls='-.', color='k', lw=1.5)
        ax4.axvline(x = 0.02817*au/R_star, ls='-.', color='k', lw=1.5)           
        #f       
        ax1.axvline(x = 0.0371*au/R_star, ls='-.', color='k', lw=1.5)
        ax2.axvline(x = 0.0371*au/R_star, ls='-.', color='k', lw=1.5)
        ax3.axvline(x = 0.0371*au/R_star, ls='-.', color='k', lw=1.5)
        ax4.axvline(x = 0.0371*au/R_star, ls='-.', color='k', lw=1.5)       
        #g
        ax1.axvline(x = 0.0451*au/R_star, ls='-.', color='k', lw=1.5)
        ax2.axvline(x = 0.0451*au/R_star, ls='-.', color='k', lw=1.5)
        ax3.axvline(x = 0.0451*au/R_star, ls='-.', color='k', lw=1.5)
        ax4.axvline(x = 0.0451*au/R_star, ls='-.', color='k', lw=1.5)       
        #h
        ax1.axvline(x = 0.063*au/R_star, ls='-.', color='k', lw=1.5)
        ax2.axvline(x = 0.063*au/R_star, ls='-.', color='k', lw=1.5)
        ax3.axvline(x = 0.063*au/R_star, ls='-.', color='k', lw=1.5)
        ax4.axvline(x = 0.063*au/R_star, ls='-.', color='k', lw=1.5)       
        ax1.text(0.01111*au/R_star,planet_height_label,'b',ha='center',fontsize=11,bbox=bbox)
        ax1.text(0.01521*au/R_star,planet_height_label,'c',ha='center',fontsize=11,bbox=bbox)    
        ax1.text(0.02144*au/R_star,planet_height_label,'d',ha='center',fontsize=11,bbox=bbox)
        ax1.text(0.02817*au/R_star,planet_height_label,'e',ha='center',fontsize=11,bbox=bbox)
        ax1.text(0.0371*au/R_star,planet_height_label,'f',ha='center',fontsize=11,bbox=bbox)
        ax1.text(0.0451*au/R_star,planet_height_label,'g',ha='center',fontsize=11,bbox=bbox)
        ax1.text(0.063*au/R_star,planet_height_label,'h',ha='center',fontsize=11,bbox=bbox) 
    
if Exoplanet=='Proxima b Turnpenney' and STUDY == 'D_ORB':     
    #d      
    ax1.axvline(x = 0.02885*au/R_star, ls='-.', color='k', lw=1.5)
    ax2.axvline(x = 0.02885*au/R_star, ls='-.', color='k', lw=1.5)
    ax3.axvline(x = 0.02885*au/R_star, ls='-.', color='k', lw=1.5)
    ax4.axvline(x = 0.02885*au/R_star, ls='-.', color='k', lw=1.5)   
    ax1.text(0.04856*au/R_star,planet_height_label,'b',ha='center',fontsize=11,bbox=bbox) 
    ax1.text(0.02885*au/R_star,planet_height_label,'d',ha='center',fontsize=11,bbox=bbox)
#secax.set_xticks([0.01521*au/R_star])
#secax.set_xticklabels(['c'])
#secax = ax.secondary_xaxis('top', functions=(spi.identity, spi.identity))
#secax.set_xticks(np.array([0.01521*au/R_star]))

fig.set_size_inches(10, 10)

#diagnostic_string = "-Bplanet" + str(B_planet_arr[loc_pl]) + "G" +'-'+'T_corona'+str(T_corona/1e6)+'MK'+'-'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'.pdf' 
#diagnostic_string = "{:.1f}".format(B_star) + "G" + "-Bplanet" + str(B_planet_arr[loc_pl]) + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'
diagnostic_string = "{:.1f}".format(B_star) + "G" + "-Bplanet" + '['+"{:.3f}".format(Bplanet_field)+']' + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'
#if Bfield_geom_arr[ind] == 'open_parker_spiral':
    #out_diagnos =  FOLDER + '/' + STUDY + "_" + str(Exoplanet.replace(" ", "_")) + "-diagnostic" + "-Open-Bstar" + diagnostic_string 
#    geometry = "-Open-spiral-Bstar" 
#elif Bfield_geom_arr[ind]== 'closed_dipole':
    #out_diagnos =  FOLDER + '/' + STUDY + "_" + str(Exoplanet.replace(" ", "_")) + "-diagnostic" + "-Closed-Bstar" + diagnostic_string 
#    geometry = "-Closed-dipole-Bstar"
#else:
#    geometry = "-Closed-PFSS-Bstar"
out_diagnos =  FOLDER + '/'+'/DIAG_PDF'+ '/'+ STUDY +"-diagnostic-"  + "_" + str(Exoplanet.replace(" ", "_")) +  geometry + diagnostic_string +'.pdf' 
plt.savefig(out_diagnos,bbox_inches='tight')

df_B_tot= pd.DataFrame({
     STUDY: x,
    'Bsw': B_sw*np.ones(len(x))
})  
df_B_tot.to_csv(FOLDER + '/CSV/' +"diagnostic-" + STUDY + "_" + str(Exoplanet.replace(" ", "_")) +  geometry + diagnostic_string +'_B_sw.csv')

    
df_M_A = pd.DataFrame({
    STUDY: x,
    'M_A': M_A*np.ones(len(x))
})  
df_M_A.to_csv(FOLDER + '/CSV/' +"diagnostic-" + STUDY + "_" + str(Exoplanet.replace(" ", "_")) +  geometry + diagnostic_string +'_M_A.csv')
















fig, ax3 = plt.subplots(1, 1, sharex=True, figsize=(9, 8))

ax3.plot(x, M_A*np.ones(len(x)), color='k', lw=lw)
ax3.set_yscale('log')  
ax3.set_xscale('log')  
ax3.set_ylabel(r"$M_A$")
ax3.set_xlabel(xlabel,fontsize=20)
ax3.set_facecolor("white")
ax3.set_ylim(1e-3,1e1)
yticks3 = ax3.get_yticks()
ax3.set_yticks(yticks3[2:-2]) 
secax3 = ax3.secondary_yaxis('right', functions=(spi.identity,spi.identity))
secax3.set_yticks(yticks3[2:-2]) 
ax3.axvline(x = xnom, ls='--', color='k', lw=2)
ax3.axvspan(x[0], R_SS, facecolor='grey', alpha=0.6)
if (M_A > 1).any():
    ax3.axhline(y = 1, ls='-.', color='grey', lw=2)   
planet_height_label=1.3e1
bbox=dict(facecolor='white', edgecolor='none', boxstyle='square,pad=0.1')
#ax3.set_xlim([1,3e3])
if Exoplanet=='Trappist-1 b' or Exoplanet=='Trappist-1 b Reville':
    if STUDY == 'D_ORB': 
        #ax3.set_xlim([1,3e3])
        #c      
        ax3.axvline(x = 0.01521*au/R_star, ls='-.', color='k', lw=1.5)     
        #d
        ax3.axvline(x = 0.02144*au/R_star, ls='-.', color='k', lw=1.5)   
        #e       
        ax3.axvline(x = 0.02817*au/R_star, ls='-.', color='k', lw=1.5)        
        #f       
        ax3.axvline(x = 0.0371*au/R_star, ls='-.', color='k', lw=1.5)   
        #g
        ax3.axvline(x = 0.0451*au/R_star, ls='-.', color='k', lw=1.5)      
        #h
        ax3.axvline(x = 0.063*au/R_star, ls='-.', color='k', lw=1.5)
        ax3.text(0.01111*au/R_star,planet_height_label,'b',ha='center',fontsize=13,bbox=bbox)
        ax3.text(0.01521*au/R_star,planet_height_label,'c',ha='center',fontsize=13,bbox=bbox)    
        ax3.text(0.02144*au/R_star,planet_height_label,'d',ha='center',fontsize=13,bbox=bbox)
        ax3.text(0.02817*au/R_star,planet_height_label,'e',ha='center',fontsize=13,bbox=bbox)
        ax3.text(0.0371*au/R_star,planet_height_label,'f',ha='center',fontsize=13,bbox=bbox)
        ax3.text(0.0451*au/R_star,planet_height_label,'g',ha='center',fontsize=13,bbox=bbox)
        ax3.text(0.063*au/R_star,planet_height_label,'h',ha='center',fontsize=13,bbox=bbox) 

out_diagnos =  FOLDER + '/'+'/DIAG_PDF'+ '/'+ STUDY +"-diagnostic-"  + "_" + str(Exoplanet.replace(" ", "_")) +  geometry + diagnostic_string +'_only_MA.pdf' 
plt.savefig(out_diagnos,bbox_inches='tight')
plt.close(fig)





'''
# Create new figure
fig2, ax3_new = plt.subplots()

# Copy all elements from ax3 to ax3_new
for line in ax3.get_lines():
    ax3_new.plot(line.get_xdata(), line.get_ydata(), label=line.get_label(), 
                 color=line.get_color(), linestyle=line.get_linestyle(), linewidth=line.get_linewidth())

# Copy axis labels and limits
ax3_new.set_xlim(ax3.get_xlim())
ax3_new.set_ylim(ax3.get_ylim())
ax3_new.set_xlabel(ax3.get_xlabel())
ax3_new.set_ylabel(ax3.get_ylabel())
ax3_new.set_yscale(ax3.get_yscale())

ax3_new.set_yscale('log')  
ax3_new.set_xscale('log')  
# Optional: copy legend, title, etc. if needed
# ax2_new.legend()

# Save
#fig2.savefig("ax2_only.pdf", bbox_inches='tight', pad_inches=0.0)

'''


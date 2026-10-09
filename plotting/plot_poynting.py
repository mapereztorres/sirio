import numpy as np

from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.patches as mpatches
#matplotlib.rc_file_defaults()
plt.style.use(['bmh','SPIworkflow/spi.mplstyle'])


#Seeting up the plots
figure=plt.figure(figsize=(8,7.5))
ax2 = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
ax2.set_facecolor("white")




y_min = S_poynt*1e-7 # minimum flux (array), Saur/Turnpenney model
y_max = S_poynt*1e-7 # maximum flux (array)
y_inter = S_poynt*1e-7
y_min_reconnect = S_reconnect*1e-7
y_max_reconnect = S_reconnect*1e-7
y_inter_reconnect = S_reconnect*1e-7
y_min_Z = Flux_r_S_Z_min # minimum flux (array), Zarka model
y_max_Z = Flux_r_S_Z_max # maximum flux (array)
 


indices_Flux_larger_rms = np.argwhere(Flux_r_S_min > 3*RMS)
indices_Flux_smaller_rms = np.argwhere(Flux_r_S_max < 3*RMS)
if indices_Flux_larger_rms.size > 0:
    x_larger_rms = x[indices_Flux_larger_rms[0]]
    x_larger_rms=x_larger_rms[0]
    x_larger_rms="{:.2f}".format(x_larger_rms)    
    x_larger_rms=str(x_larger_rms)      
    
    x_last_larger=x[indices_Flux_larger_rms[-1]]
    x_last_larger=x_last_larger[0]
    x_last_larger="{:.2f}".format(x_last_larger)    
    x_last_larger=str(x_last_larger)
    #print('value of x where there is clear detection for the Alfvén Wing model: ( ',x_larger_rms+' , '+x_last_larger+' )')
    x_larger_rms=x_larger_rms+' , '+x_last_larger
else:
    x_larger_rms=np.nan
    x_larger_rms=str(x_larger_rms)


if indices_Flux_smaller_rms.size > 0:
    x_smaller_rms = x[indices_Flux_smaller_rms[0]]
    x_smaller_rms=x_smaller_rms[0]
    x_smaller_rms="{:.2f}".format(x_smaller_rms)    
    x_smaller_rms=str(x_smaller_rms)      
    
    x_last_smaller=x[indices_Flux_smaller_rms[-1]]
    x_last_smaller=x_last_smaller[0]
    x_last_smaller="{:.2f}".format(x_last_smaller)    
    x_last_smaller=str(x_last_smaller)
    #print('value of x where there is clear NON detection for the Alfvén Wing model: ( ',x_smaller_rms+' , '+x_last_smaller+' )')
    x_smaller_rms=x_smaller_rms+' , '+x_last_smaller
else:
    x_smaller_rms=np.nan
    x_smaller_rms=str(x_smaller_rms)


ax2.plot(x,y_inter,color='black',lw=1.5)




if STUDY == 'D_ORB':
    ax2.set_yscale('log') 
    # Draw vertical line at nominal orbital separation of planet
    xnom = r_orb/R_star
    #print('Planet '+Exoplanet+' at an orbital separation of '+ str(xnom))
    xlabel=r"Orbital separation / Stellar radius"
    if PLOT_M_A == True:
        ax0.axvline(x = xnom, ls='--', color='k', lw=2)
    ax2.axvline(x = xnom, ls='--', color='k', lw=2)
    ax2.set_xlabel(xlabel,fontsize=20)
    ax1 = ax2.twiny()
    ax1.set_xlabel(r"Orbital period (days)")
    
    if Bfield_geom_arr[ind] == 'pfss':
        ax2.axvspan(x[0], R_SS, facecolor='gray', alpha=0.7)
    def tick_function(X):
        V = spi.Kepler_P(M_star/M_sun,X*R_star/au)
        return ["%.1f" % z for z in V]
    xtickslocs = ax2.get_xticks()    
    new_tick_locations=xtickslocs[1:-1]
    
    ax1.set_xlim(ax2.get_xlim())
    ax1.set_xticks(new_tick_locations)
    ax1.set_xticklabels(tick_function(new_tick_locations))
    label_location='upper right'       
    ax2.text(lim_x[1]*0.1, 1e11, r'B$_{pl} = $'+"{:.2f}".format(Bplanet_field)+' G', fontsize = 16,bbox=dict(facecolor='white', alpha=1,edgecolor='white'))
    ax2.text(lim_x[1]*0.1, 2e11, r'T$_{c} = $'+"{:.1f}".format(T_corona/1e6)+' MK', fontsize = 16,bbox=dict(facecolor='white', alpha=1,edgecolor='white'))

elif STUDY == 'M_DOT':
    ax2.set_xscale('log') 
    ax2.set_yscale('log') 
    xnom = M_star_dot
    xlabel = r"Mass Loss rate [$\dot{M}_\odot$]"
    # Draw vertical line at nominal mass loss rate of the star
    ax2.axvline(x = xnom, ls='--', color='k', lw=2)
    ax2.set_xlabel(xlabel,fontsize=20)
    ax2.set_xlim([x[0],x[-1]])
    if magnetized_pl_arr[ind1]:
        ax2.text(1.5e-1, 10**((np.log10(YLIMHIGH)-1)*1.1), r'B$_{pl} = $'+"{:.2f}".format(Bplanet_field)+' G', fontsize = 16,bbox=dict(facecolor='white', alpha=1,edgecolor='white'))
    else:
        ax2.text(1.5e-1, 10**((np.log10(YLIMHIGH)-1)*1.1), r'B$_{pl} = $'+'0 G', fontsize = 16,bbox=dict(facecolor='white', alpha=1,edgecolor='white'))

    ax2.text(1.5e-1, 10**((np.log10(YLIMHIGH)-1)*0.85), r'T$_{c} = $'+"{:.1f}".format(T_corona/1e6)+' MK', fontsize = 16,bbox=dict(facecolor='white', alpha=1,edgecolor='white'))
        label_location='upper left'   

    
elif STUDY == 'B_PL':
    ax2.set_yscale('log'); 
    #xnom = B_planet_Sano
    xnom = Bplanet_field
    xlabel = r"Planetary magnetic field [Gauss]"
    # Draw vertical line at the reference planetary magnetic field
    ax2.axvline(x = xnom, ls='--', color='k', lw=2)
    ax2.set_xlabel(xlabel,fontsize=20)
    ax2.set_xlim([x[0],x[-1]])
    ax2.text(0.1, 1.05e1, r'T$_{c} = $'+"{:.1f}".format(T_corona/1e6)+' MK', fontsize = 16,bbox=dict(facecolor='white', alpha=0,edgecolor='white'))
    label_location='upper left'   
    if xnom<1: 
        ax2.set_xlim([0,1])
    
    
    
    
if (STUDY == 'D_ORB') or (STUDY == 'M_DOT'):
    if PLOT_M_A == True:
        ax0.set_yscale('log')                
        # Draw vertical line at the nomimal value of the x-axis
        ax0.axvline(x = xnom, ls='--', color='k', lw=2)
        if STUDY == 'M_DOT':
            ax0.set_xscale('log')
            if LIMS_MA == True:
                ax0.set_ylim((LIM_MA_LOW, LIM_MA_HIGH))
                


ax2.set_ylim([1e10,1e16])  
ax2.set_ylabel(r"S [W]")

orange_patch = mpatches.Patch(color='orange', label='Alfvén wing')
blue_patch = mpatches.Patch(facecolor='blue',label='Reconnection')

lim_x=ax2.get_xlim()
lim_y=ax2.get_ylim()



# Draw 3*RMS upper limit?
if DRAW_RMS == True:
    ax2.axhline(y = 3*RMS, ls='-.', color='grey', lw=2)

# Draw a little Earth at the planet position for visualization purposes?

#Print out relevant input and output parameters, including the expected flux received at Earth 
# from the SPI at the position of the planet
# To this end, first find out the position of the planet in the distance array
d_diff = np.abs((d_orb - r_orb) / R_star)
loc_pl = np.where(d_diff == d_diff.min())


# Print in the graph the value of the planetary magnetic field, in units of bfield_earth
if STUDY == 'B_PL':
    B_planet_ref = round(float(B_planet_Sano[0] /(bfield_earth * Tesla2Gauss) ), 2) 
else:
    B_planet_ref = round(float(B_planet_arr[loc_pl][0] / (bfield_earth*Tesla2Gauss) ), 2) 



common_string = "{:.1f}".format(B_star) + "G" + "-Bplanet" +'['+"{:.3f}".format(Bplanet_field)+']' + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'             

outfile =  STUDY +  "_" + str(Exoplanet.replace(" ", "_")) + geometry + common_string   

if freefree == True:
    outfile = outfile + '_freefree'


if PLOTOUT == True:
    plt.tight_layout()
    outfilePDF = os.path.join(FOLDER + '/S_poynting/' +'S_poynt_'+outfile+ ".pdf")
    plt.savefig(outfilePDF, bbox_inches='tight', pad_inches=0)
    plt.close()
else:
    plt.tight_layout()
    plt.show()
    
   

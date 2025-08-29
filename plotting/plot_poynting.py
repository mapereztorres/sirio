import numpy as np

from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.patches as mpatches
#matplotlib.rc_file_defaults()
plt.style.use(['bmh','SPIworkflow/spi.mplstyle'])

### Plot received flux density as a function of distance from the star
###
#plt.style.use(['bmh', '/home/torres/Dropbox/python/styles/paper.mplstyle'])

################################
# lw, PLOT_MA, plotout taken as defined in __init__.py
lw = LW 

# Kepler's third law, with d_orb_mark in units of R_star, 
# so that period_mark is in days.
period_mark = np.array([1, 10, 30, 60, 100, 200, 500, 1000, 2000])
d_orb_mark = (period_mark/yr)**(2/3) * (M_star/M_sun)**(1/3) * (au/R_star)

# Plotting is different, depending on the "STUDY" case
if STUDY == 'D_ORB':
    x = d_orb / R_star # (distance array, in units of R_star)
elif STUDY == 'M_DOT':
    x = M_star_dot_arr # (M_star_dot_arr array, in units of M_dot_sun)
elif STUDY == 'B_PL':
    x = B_planet_arr # (B_planet_arr array, in Gauss )

if (STUDY == 'D_ORB') or (STUDY == 'M_DOT'):
    if PLOT_M_A == True:
        figure=plt.figure(figsize=(8,11))
        ax0 = plt.subplot2grid((3,1),(0,0),rowspan=1,colspan=1)
        ax2 = plt.subplot2grid((3,1),(1,0),rowspan=2,colspan=1)
        ax0.plot(x, M_A, color='k', lw=lw)
        ax0.set_ylabel(r"$M_A$")
        ax0.set_facecolor("white")
    else:
        figure=plt.figure(figsize=(8,7.5))
        ax2 = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
        ax2.set_facecolor("white")             
        
    ax2.set_facecolor("white")	
    
elif STUDY == 'B_PL':
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

x_superalfv='nan'           
if any(ind > 1 for ind in M_A):
    #print('The planet enters super-Afvénic regime')
    M_A_superalfv_arr=np.where(M_A >1)
    M_A_superalfv_ind=M_A_superalfv_arr[0]
    M_A_superalfv_ind=M_A_superalfv_ind[0]
    x_superalfv=x[M_A_superalfv_ind]
    if PLOT_M_A == True:
        ax0.axvline(x = x_superalfv, color='grey',lw=2)
        ax0.axvspan(x_superalfv, x[-1], facecolor='grey', alpha=0.5)
    if x_superalfv!=x[0]: 
        ax2.axvline(x = x_superalfv, color='black',lw=2)

    ax2.axvspan(x_superalfv, x[-1], facecolor='black', alpha=0.6,hatch='x')
    ax2.set_xlim(1,x_superalfv)


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
        # ax2.axvline(x = R_SS, ls='--', color='k', lw=2)
        #ax2.axvspan(x[0], R_SS, facecolor='grey', alpha=0.5,hatch='x')
        ax2.axvspan(x[0], R_SS, facecolor='gray', alpha=0.7)
    def tick_function(X):
        V = spi.Kepler_P(M_star/M_sun,X*R_star/au)
        return ["%.1f" % z for z in V]
    xtickslocs = ax2.get_xticks()    
    new_tick_locations=xtickslocs[1:-1]

    ax1.set_xlim(ax2.get_xlim())
    ax1.set_xticks(new_tick_locations)
    ax1.set_xticklabels(tick_function(new_tick_locations))

elif STUDY == 'M_DOT':
    ax2.set_xscale('log') 
    ax2.set_yscale('log') 
    xnom = M_star_dot
    xlabel = r"Mass Loss rate [$\dot{M}_\odot$]"
    # Draw vertical line at nominal mass loss rate of the star
    ax2.axvline(x = xnom, ls='--', color='k', lw=2)
    ax2.set_xlabel(xlabel,fontsize=20)
    ax2.set_xlim([x[0],x[-1]])

    
elif STUDY == 'B_PL':
    ax2.set_yscale('log'); 
    #xnom = B_planet_Sano
    xnom = Bplanet_field
    xlabel = r"Planetary magnetic field [Gauss]"
    # Draw vertical line at the reference planetary magnetic field
    ax2.axvline(x = xnom, ls='--', color='k', lw=2)
    ax2.set_xlabel(xlabel,fontsize=20)
    ax2.set_xlim([x[0],x[-1]])
    
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


if STUDY == "D_ORB":
    label_location='upper right'       
    ax2.text(lim_x[1]*0.1, 1e11, r'B$_{pl} = $'+"{:.2f}".format(Bplanet_field)+' G', fontsize = 16,bbox=dict(facecolor='white', alpha=1,edgecolor='white'))
    ax2.text(lim_x[1]*0.1, 2e11, r'T$_{c} = $'+"{:.1f}".format(T_corona/1e6)+' MK', fontsize = 16,bbox=dict(facecolor='white', alpha=1,edgecolor='white'))


elif STUDY == "M_DOT":

        if magnetized_pl_arr[ind1]:
            ax2.text(1.5e-1, 10**((np.log10(YLIMHIGH)-1)*1.1), r'B$_{pl} = $'+"{:.2f}".format(Bplanet_field)+' G', fontsize = 16,bbox=dict(facecolor='white', alpha=1,edgecolor='white'))

        else:

            ax2.text(1.5e-1, 10**((np.log10(YLIMHIGH)-1)*1.1), r'B$_{pl} = $'+'0 G', fontsize = 16,bbox=dict(facecolor='white', alpha=1,edgecolor='white'))

        ax2.text(1.5e-1, 10**((np.log10(YLIMHIGH)-1)*0.85), r'T$_{c} = $'+"{:.1f}".format(T_corona/1e6)+' MK', fontsize = 16,bbox=dict(facecolor='white', alpha=1,edgecolor='white'))
        label_location='upper left'   
        
elif STUDY == "B_PL":  

    ax2.text(0.1, 1.05e1, r'T$_{c} = $'+"{:.1f}".format(T_corona/1e6)+' MK', fontsize = 16,bbox=dict(facecolor='white', alpha=0,edgecolor='white'))
    label_location='upper left'   


if STUDY == "B_PL" and xnom<1: 

    ax2.set_xlim([0,1])

# Draw 3*RMS upper limit?
if DRAW_RMS == True:
    ax2.axhline(y = 3*RMS, ls='-.', color='grey', lw=2)



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





if Exoplanet=='YZCet b Model A' or Exoplanet=='YZCet b Model B':
    if (STUDY == 'M_DOT') :
        ax2.axvline(x = 0.25, ls='--', color='k', lw=2)
        ax2.axvline(x = 5, ls='--', color='k', lw=2)

    if (STUDY == 'B_PL') :
        ax2.text(1.2,3e-3,'Model '+Exoplanet[-1])

secax = ax2.secondary_yaxis('right', functions=(spi.identity,spi.identity))


if Exoplanet=='Proxima b kavanagh'and STUDY=='D_ORB':
    ax2.set_xlim([1,140])




common_string = "{:.1f}".format(B_star) + "G" + "-Bplanet" +'['+"{:.3f}".format(Bplanet_field)+']' + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'             

outfile =  STUDY +  "_" + str(Exoplanet.replace(" ", "_")) + geometry + common_string   

if freefree == True:
    outfile = outfile + '_freefree'


if PLOTOUT == True:
    plt.tight_layout()
    #figure.subplots_adjust(top=0.98, bottom=0.1, left=0.1, right=0.95)
    outfilePDF = os.path.join(FOLDER + '/S_poynting/' +'S_poynt_'+outfile+ ".pdf")
    plt.savefig(outfilePDF, bbox_inches='tight', pad_inches=0)
    #outfilePNG = os.path.join(FOLDER + '/FLUX_PNG/' +'S_poynt_'+outfile +".png")
    #plt.savefig(outfilePNG, bbox_inches='tight')
    plt.close()
else:
    plt.tight_layout()
    plt.show()
    
   

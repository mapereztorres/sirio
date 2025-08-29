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
        plt.figure(figsize=(8,11))
        ax1 = plt.subplot2grid((3,1),(0,0),rowspan=1,colspan=1)
        ax2 = plt.subplot2grid((3,1),(1,0),rowspan=2,colspan=1)
        ax1.plot(x, M_A, color='k', lw=lw)
        ax1.set_ylabel(r"$M_A$")
        ax1.set_facecolor("white")
    else:
        plt.figure(figsize=(8,7.5))
        ax2 = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
        ax2.set_facecolor("white")             
        
    ax2.set_facecolor("white")	
    
elif STUDY == 'B_PL':
    plt.figure(figsize=(8,7.5))
    ax2 = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
    ax2.set_facecolor("white")	

y_min = Flux_r_S_min # minimum flux (array), Saur/Turnpenney model
y_max = Flux_r_S_max # maximum flux (array)
y_inter = Flux_r_S_inter
y_min_reconnect = Flux_reconnect_min
y_max_reconnect = Flux_reconnect_max
y_inter_reconnect = Flux_reconnect_inter
y_min_Z = Flux_r_S_Z_min # minimum flux (array), Zarka model
y_max_Z = Flux_r_S_Z_max # maximum flux (array)

y_min_sb = Flux_sb_min
y_max_sb = Flux_sb_max
y_inter_sb = Flux_sb_inter

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

ax2.fill_between(x, y_min, y_max, color="orange", alpha=0.7, hatch='///', edgecolor='orange')
ax2.fill_between(x, y_min_reconnect, y_max_reconnect, color="blue", alpha=0.7, hatch='\\\\\\', edgecolor='blue')
ax2.fill_between(x, y_min_sb, y_max_sb, color="green", alpha=0.7, hatch='xxx', edgecolor='green')
ax2.plot(x,y_inter,color='black',lw=1.5)
ax2.plot(x,y_inter_reconnect,color='black',lw=1.5)
ax2.plot(x,y_inter_sb,color='black',lw=1.5)
if STUDY == 'D_ORB':
    ax2.set_yscale('log') 
    # Draw vertical line at nominal orbital separation of planet
    xnom = r_orb/R_star
    #xlabel=r"Distance / Stellar radius"
    xlabel=r"Orbital separation / Stellar radius"
    if PLOT_M_A == True:
        ax1.axvline(x = xnom, ls='--', color='k', lw=2)
    ax2.axvline(x = xnom, ls='--', color='k', lw=2)
    ax2.set_xlabel(xlabel,fontsize=20)
    ax2.set_xlim(1,d_orb_max)

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
        ax1.set_yscale('log')                
        # Draw vertical line at the nomimal value of the x-axis
        ax1.axvline(x = xnom, ls='--', color='k', lw=2)
        if STUDY == 'M_DOT':
            ax1.set_xscale('log')
            if LIMS_MA == True:
                ax1.set_ylim((LIM_MA_LOW, LIM_MA_HIGH))
                



   
ax2.set_ylabel(r"Flux density [mJy]")

orange_patch = mpatches.Patch(facecolor='orange', edgecolor='orange', hatch='///', label="Alfvén wing")
blue_patch = mpatches.Patch(facecolor='blue', edgecolor='blue', hatch='\\\\\\', label='Reconnection')
green_patch = mpatches.Patch(facecolor='green', edgecolor='green', hatch='xxx', label='Stretch and break')


if STUDY == "D_ORB":

    label_location = 'upper right'

    ax2.text(d_orb_max * 0.7, 10 ** ((np.log10(YLIMHIGH) - 1) * 0.7),
             r'B$_{pl} = $' + "{:.2f}".format(Bplanet_field) + ' G', fontsize=16,
             bbox=dict(facecolor='white', alpha=1, edgecolor='white'))
    ax2.text(d_orb_max * 0.7, 10 ** ((np.log10(YLIMHIGH) - 1) * 0.45),
             r'T$_{c} = $' + "{:.1f}".format(T_corona / 1e6) + ' MK', fontsize=16,
             bbox=dict(facecolor='white', alpha=1, edgecolor='white'))


elif STUDY == "M_DOT":

    if magnetized_pl_arr[ind1]:
        ax2.text(1.5e-1, 10 ** ((np.log10(YLIMHIGH) - 1) * 0.8), r'B$_{pl} = $' + "{:.2f}".format(Bplanet_field) + ' G',
                 fontsize=16, bbox=dict(facecolor='white', alpha=1, edgecolor='white'))

    else:

        ax2.text(1.5e-1, 10 ** ((np.log10(YLIMHIGH) - 1) * 0.8), r'B$_{pl} = $' + '0 G', fontsize=16,
                 bbox=dict(facecolor='white', alpha=1, edgecolor='white'))


    ax2.text(1.5e-1, 10 ** ((np.log10(YLIMHIGH) - 1) * 0.45),
             r'T$_{c} = $' + "{:.1f}".format(T_corona / 1e6) + ' MK', fontsize=16,
             bbox=dict(facecolor='white', alpha=1, edgecolor='white'))
    label_location = 'upper left'

elif STUDY == "B_PL":

    ax2.text(0.1, 0.5e1, r'T$_{c} = $' + "{:.1f}".format(T_corona / 1e6) + ' MK', fontsize=16,
             bbox=dict(facecolor='white', alpha=1, edgecolor='white'))
    label_location = 'upper left'
    ax2.legend(handles=[blue_patch, orange_patch, green_patch], loc=label_location, fontsize=22, facecolor='white',
           edgecolor='white', framealpha=1)

if STUDY == "B_PL" and xnom < 1:
    # ax2.set_xscale('log')
    ax2.set_xlim([0, 1])
    
ax2.legend(handles=[green_patch,blue_patch,orange_patch],loc=label_location,fontsize=25,facecolor='white',edgecolor='white', framealpha=1)



# Draw 3*RMS upper limit?
if DRAW_RMS == True:
    ax2.axhline(y = 3*RMS, ls='-.', color='grey', lw=2)

# Draw a little Earth at the planet position for visualization purposes?
if (DRAW_EARTH == True) and (STUDY == 'D_ORB'):
    paths = ['./pics/earth.png']
    x_earth = [r_orb / R_star]
    y = [3*RMS]
    for x0, y0, path in zip(x_earth, y, paths):
        ab_earth = AnnotationBbox(spi.getImage(path), (x0, y0), frameon=False)
        ax2.add_artist(ab_earth)            

#Print out relevant input and output parameters, including the expected flux received at Earth 
# from the SPI at the position of the planet
# To this end, first find out the position of the planet in the distance array
d_diff = np.abs((d_orb - r_orb) / R_star)
loc_pl = np.where(d_diff == d_diff.min())
#if STUDY == 'D_ORB':
#print('Position in d_orb array where the planet is located', loc_pl)

# Print in the graph the value of the planetary magnetic field, in units of bfield_earth
if STUDY == 'B_PL':
    B_planet_ref = round(float(B_planet_Sano[0] /(bfield_earth * Tesla2Gauss) ), 2) 
else:
    B_planet_ref = round(float(B_planet_arr[loc_pl][0] / (bfield_earth*Tesla2Gauss) ), 2) 

# 
x_superalfv='nan'           
if any(ind > 1 for ind in M_A):
    #print('The planet enters super-Afvénic regime')
    M_A_superalfv_arr=np.where(M_A >1)
    M_A_superalfv_ind=M_A_superalfv_arr[0]
    M_A_superalfv_ind=M_A_superalfv_ind[0]
    #mdot_superalfv=M_star_dot_arr[M_A_superalfv_ind]
    x_superalfv=x[M_A_superalfv_ind]
    if PLOT_M_A == True:
        ax1.axvline(x = x_superalfv, color='grey',lw=2)
        ax1.axvspan(x_superalfv, x[-1], facecolor='grey', alpha=0.5)
    if x_superalfv!=x[0]: 
        ax2.axvline(x = x_superalfv, color='black',lw=2)
    #ax2.axvspan(x_superalfv, x[-1], facecolor='grey', alpha=0.5)
    ax2.axvspan(x_superalfv, x[-1], facecolor='black', alpha=0.6,hatch='x')
    #print(f'For the study {STUDY}, planet enters a superalfvénic regime at value {STUDY}',x_superalfv)



if (STUDY == 'M_DOT') :
    if Exoplanet=='YZCet b Model A' or Exoplanet=='YZCet b Model B':
        ax2.axvline(x = 0.25, ls='--', color='k', lw=2)
        ax2.axvline(x = 5, ls='--', color='k', lw=2)
        ax2.text(0.25, 1e2, 'A')
        ax2.text(5, 1e2, 'B')



secax = ax2.secondary_yaxis('right', functions=(spi.identity,spi.identity))
ax2.set_ylim([YLIMLOW, YLIMHIGH])       

for ax in [ax2,secax]:
    ax.tick_params(axis='both', which='major', labelsize=35, width=2, length=8)  # bigger ticks
    ax.tick_params(axis='both', which='minor', labelsize=35, width=1.5, length=5)  
    ax.xaxis.label.set_size(30)
    ax.yaxis.label.set_size(30)

    # Make axis spines (borders) thicker
    for spine in ax.spines.values():
        spine.set_linewidth(4)

if STUDY == 'D_ORB':
    ax2.set_ylim([1e-5, 1e4])      
    ax2.set_xscale('log') 
    ax2.set_xlim([9,2e2])
    
if STUDY == "B_PL" and xnom<1: 
    #ax2.set_xscale('log') 
    ax2.set_xlim([0,1])


    
    

common_string = "{:.1f}".format(B_star) + "G" + "-Bplanet" +'['+"{:.3f}".format(Bplanet_field)+']' + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'             


outfile =  STUDY +  "_" + str(Exoplanet.replace(" ", "_")) + geometry + common_string   
# Variable to send output to files (PLOTOUT= True), or show them in
# the terminal (PLOTOUT = False) 

if freefree == True:
    outfile = outfile + '_freefree'


if PLOTOUT == True:
    plt.tight_layout()
    outfilePDF = os.path.join(FOLDER + '/FLUX_PDF/' +'Flux_'+outfile+ ".pdf")
    plt.savefig(outfilePDF)
    outfilePNG = os.path.join(FOLDER + '/FLUX_PNG/' +'Flux_'+outfile +".png")
    plt.savefig(outfilePNG)
    plt.close()
else:
    plt.tight_layout()
    plt.show()
    
   

df_alfven = pd.DataFrame({
     STUDY: x,
    'FLUX': y_inter
})  
df_alfven.to_csv(FOLDER + '/CSV/' +outfile+'_alfven_wing_model.csv')

df_reconnect = pd.DataFrame({
    STUDY: x,
    'FLUX': y_inter_reconnect
})  
df_reconnect.to_csv(FOLDER + '/CSV/' +outfile+'_reconnection_model.csv')  

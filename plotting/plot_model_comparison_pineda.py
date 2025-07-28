import numpy as np

from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
#matplotlib.rc_file_defaults()
plt.style.use(['bmh','SPIworkflow/spi.mplstyle'])

path = FOLDER.replace("YZCet_b_Model_B", "YZ_Cet_Pineda")
print(path)
#print(a)
if not(os.path.isdir(path)):
    os.system('mkdir '+path)

#if not(os.path.isdir('/home/luis/github/spirou/OUTPUT/YZ_Cet_Pineda')):
#    os.system('mkdir /home/luis/github/spirou/OUTPUT/YZ_Cet_Pineda')
print('Generating Pineda-like plots')
##comparison plots for the fluxes
outfile =  STUDY +  "_" + str(Exoplanet.replace(" ", "_")) + geometry + common_string
df_alfven.to_csv(FOLDER + '/CSV/' +outfile+'_alfven_wing_model.csv')

infile = STUDY + "_" + str(Exoplanet.replace(" ", "_")) + geometry + common_string + '_alfven_wing_model.csv'
alfven_wing_parker_A = pd.read_csv(FOLDER + '/CSV/' + infile)

FOLDER_A=FOLDER
FOLDER_A=FOLDER_A.replace("YZCet_b_Model_B", "YZCet_b_Model_A")
#alfven_wing_parker_A = pd.read_csv(FOLDER_A + '/CSV/' + STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-open-parker-spiral' + common_string + '_alfven_wing_model.csv' )
#reconnection_parker_A = pd.read_csv(FOLDER_A + '/CSV/' + STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-open-parker-spiral' + common_string + '_reconnection_model.csv')
alfven_wing_parker_A = pd.read_csv(FOLDER_A + '/CSV/D_ORB_YZCet_b_Model_A-open-parker-spiral-Bstar220.0G-Bplanet[0.197]G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_alfven_wing_model.csv' )
reconnection_parker_A= pd.read_csv(FOLDER_A + '/CSV/D_ORB_YZCet_b_Model_A-open-parker-spiral-Bstar220.0G-Bplanet[0.197]G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_reconnection_model.csv' )
FOLDER_B=FOLDER
#FOLDER_B=FOLDER_A.replace("YZCet_b_Model_B", "YZCet_b_Model_B")
#alfven_wing_parker_B = pd.read_csv(FOLDER_B + '/CSV/' + STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-open-parker-spiral' + common_string + '_alfven_wing_model.csv' )
#reconnection_parker_B = pd.read_csv(FOLDER_B + '/CSV/' + STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-open-parker-spiral' + common_string + '_reconnection_model.csv')
alfven_wing_parker_B = pd.read_csv(FOLDER_B + '/CSV/D_ORB_YZCet_b_Model_B-open-parker-spiral-Bstar220.0G-Bplanet[0.197]G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_alfven_wing_model.csv' )
reconnection_parker_B= pd.read_csv(FOLDER_B + '/CSV/D_ORB_YZCet_b_Model_B-open-parker-spiral-Bstar220.0G-Bplanet[0.197]G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_reconnection_model.csv' )

#alfven_wing_parker_A = pd.read_csv('/home/luis/github/spirou/OUTPUT/YZCet_b_Model_A/CSV/D_ORB_YZCet_b_Model_A-open-parker-spiral-Bsta#r220.0G-Bplanet[0.197]G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_freefree_alfven_wing_model.csv')
#reconnection_parker_A = pd.read_csv('/home/luis/github/spirou/OUTPUT/YZCet_b_Model_A/CSV/D_ORB_YZCet_b_Model_A-open-parker-spiral-Bstar220.0G-Bplanet[0.197]G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_freefree_reconnection_model.csv')


#alfven_wing_pfss_parker_B = pd.read_csv('/home/luis/github/spirou/OUTPUT/YZCet_b_Model_B/CSV/D_ORB_YZCet_b_Model_B-open-parker-spiral-Bstar220.0G-Bplanet[0.197]G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_freefree_alfven_wing_model.csv')
#reconnection_pfss_parker_B = pd.read_csv('/home/luis/github/spirou/OUTPUT/YZCet_b_Model_B/CSV/D_ORB_YZCet_b_Model_B-open-parker-spiral-Bstar220.0G-Bplanet[0.197]G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_freefree_reconnection_model.csv')
alfven_wing_pfss_parker_B = pd.read_csv(FOLDER_B + '/CSV/D_ORB_YZCet_b_Model_B-pfss-Bstar220.0G-Bplanet[0.197]G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_freefree_alfven_wing_model.csv')
reconnection_pfss_parker_B = pd.read_csv(FOLDER_B + '/CSV/D_ORB_YZCet_b_Model_B-pfss-Bstar220.0G-Bplanet[0.197]G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_alfven_wing_model.csv' )


plt.figure(figsize=(10, 7.5))
lw=2
ax2 = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)
ax2.set_facecolor("white")    
ax2.set_ylim([1e-2, 5e2])     
ax2.set_yscale('log')
ax2.plot(alfven_wing_parker_A[STUDY], alfven_wing_parker_A['FLUX'], color='orange')
ax2.plot(reconnection_parker_A[STUDY], reconnection_parker_A['FLUX'], color='blue')
ax2.plot(alfven_wing_pfss_parker_B[STUDY], alfven_wing_pfss_parker_B['FLUX'], color='orange', linestyle='dotted')
ax2.plot(reconnection_pfss_parker_B[STUDY], reconnection_pfss_parker_B['FLUX'], color='blue', linestyle='dotted')

legend_elements = [
Line2D([0], [0], color='orange', lw=lw, label='Alfvén wing (Model A)'),
Line2D([0], [0], color='blue', lw=lw, label='Reconnection (Model A)'),
Line2D([0], [0], color='orange', linestyle='dotted', lw=lw, label='Alfvén wing (Model B)'),
Line2D([0], [0], color='blue', linestyle='dotted', lw=lw, label='Reconnection (Model B)'),

]


ax2.set_xlabel(xlabel,fontsize=20)
ax2.set_ylabel(r"Flux density [mJy]")
ax2.legend(handles=legend_elements, loc='lower left', fontsize=12, facecolor='white', edgecolor='white', framealpha=0)
if STUDY == "D_ORB":
    ax2.set_xlim(left=1,right=60)       
ax2.axvspan(x[0], R_SS, facecolor='gray', alpha=0.7)
secax = ax2.secondary_yaxis('right', functions=(spi.identity,spi.identity))
ax2.axvline(x = xnom, ls='--', color='k', lw=2)
if STUDY == 'MDOT':
    ax2.set_xscale('log')
plt.savefig(path+'/'+'Flux'+'_model_comparison-Pineda.pdf', bbox_inches='tight')



##comparison plots for the stellar wind magnetic field
#bsw_parker_A = pd.read_csv('/home/luis/github/spirou/OUTPUT/YZCet_b_Model_A/CSV/diagnostic-D_ORB_YZCet_b_Model_A-open-parker-spiral-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_B_sw.csv')
#bsw_pfss_B = pd.read_csv('/home/luis/github/spirou/OUTPUT/YZCet_b_Model_B/CSV/diagnostic-D_ORB_YZCet_b_Model_B-pfss-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_B_sw.csv')
#bsw_BASE = pd.read_csv('/home/luis/github/spirou/OUTPUT/YZCet_b_Model_A/CSV/diagnostic-D_ORB_YZCet_b_Model_A-pfss-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_B_sw.csv') #M_A_PFSS_A

bsw_parker_A = pd.read_csv(FOLDER_A + '/CSV/diagnostic-D_ORB_YZCet_b_Model_A-open-parker-spiral-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_B_sw.csv')
bsw_pfss_B = pd.read_csv(FOLDER_B + '/CSV/diagnostic-D_ORB_YZCet_b_Model_B-pfss-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_B_sw.csv')
bsw_BASE = pd.read_csv(FOLDER_A + '/CSV/diagnostic-D_ORB_YZCet_b_Model_A-pfss-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_B_sw.csv') #M_A_PFSS_A


#alfven_wing_parker_B = pd.read_csv(FOLDER_B + '/CSV/' + STUDY + "_" + str(Exoplanet.replace(" ", "_")) + geometry + common_string + '_alfven_wing_model.csv' )
#reconnection_parker_B = pd.read_csv(FOLDER_B + '/CSV/' + STUDY + "_" + str(Exoplanet.replace(" ", "_")) + geometry + common_string + '_reconnection_model.csv')
plt.figure(figsize=(11, 7.5))
lw=2
ax2 = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)
ax2.set_facecolor("white")  

ax2.set_yscale('log')

ax2.plot(bsw_parker_A[STUDY], bsw_parker_A['Bsw'], color='blue', linestyle='solid')
ax2.plot(bsw_pfss_B[STUDY], bsw_pfss_B['Bsw'], color='purple', linestyle='dotted',lw=3)
ax2.plot(bsw_BASE[STUDY], bsw_BASE['Bsw'], color='orange', linestyle='dashdot')

print('BSW: Model A / MODEL B: ', bsw_parker_A['Bsw']/bsw_pfss_B['Bsw'])

legend_elements = [
Line2D([0], [0], color='blue', linestyle='solid', lw=lw, label='Model A (Parker spiral)'),
Line2D([0], [0], color='purple', linestyle='dotted', lw=lw, label='Model B (PFSS)'),
Line2D([0], [0], color='orange', linestyle='dashdot', lw=lw, label='PFSS Base (A PFSS)', alpha=0.5),
]

ax2.legend(handles=legend_elements, loc='upper right', fontsize=12, facecolor='white', edgecolor='white', framealpha=1)
if STUDY == "D_ORB":
    ax2.set_xlim(left=1,right=60)      
ax2.set_xlabel(xlabel,fontsize=20)
ax2.set_ylabel(r"$B_{\rm sw}$ $[G]$")
ax2.axvspan(x[0], R_SS, facecolor='gray', alpha=0.7)
secax = ax2.secondary_yaxis('right', functions=(spi.identity,spi.identity))
#ax2.axvline(x = xnom, ls='--', color='k', lw=2)
ax2.axvline(x = 0.01634*au/R_star, ls='--', color='k', lw=2)
ax2.axvline(x = 0.02156*au/R_star, ls='-.', color='k', lw=2)
ax2.axvline(x = 0.02851*au/R_star, ls='-.', color='k', lw=2)

#ax2.text(xnom,5e2,'b',ha='center',fontsize=11)
#ax2.text(0.01634*au/R_star,max(bsw_parker_A['Bsw'])*1.01,'c',ha='center',fontsize=13)
#ax2.text(0.02156*au/R_star,max(bsw_parker_A['Bsw'])*1.01,'c',ha='center',fontsize=13)
#ax2.text(0.02851*au/R_star,max(bsw_parker_A['Bsw'])*1.01,'d',ha='center',fontsize=13)
#ax2.text(0.01634*au/R_star,4e2,'b',ha='center',fontsize=13)
#ax2.text(0.02156*au/R_star,4e2,'c',ha='center',fontsize=13)
#ax2.text(0.02851*au/R_star,4e2,'d',ha='center',fontsize=13)
ax2.text(0.01634*au/R_star-1,2e2,'b',ha='center',fontsize=13)
ax2.text(0.02156*au/R_star-1,2e2,'c',ha='center',fontsize=13)
ax2.text(0.02851*au/R_star-1,2e2,'d',ha='center',fontsize=13)

plt.savefig(path+'/'+'B_sw_'+'_model_comparison-Pineda.pdf', bbox_inches='tight')




##comparison plots for M_A
#M_A_parker_A = pd.read_csv('/home/luis/github/spirou/OUTPUT/YZCet_b_Model_A/CSV/diagnostic-D_ORB_YZCet_b_Model_A-open-parker-spiral-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_M_A.csv')
#M_A_PFSS_B = pd.read_csv('/home/luis/github/spirou/OUTPUT/YZCet_b_Model_B/CSV/diagnostic-D_ORB_YZCet_b_Model_B-pfss-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_M_A.csv')
#M_A_BASE = pd.read_csv('/home/luis/github/spirou/OUTPUT/YZCet_b_Model_A/CSV/diagnostic-D_ORB_YZCet_b_Model_A-pfss-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_M_A.csv') #M_A_PFSS_A
M_A_parker_A = pd.read_csv(FOLDER_A + '/CSV/diagnostic-D_ORB_YZCet_b_Model_A-open-parker-spiral-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_M_A.csv')
M_A_PFSS_B = pd.read_csv(FOLDER_B + '/CSV/diagnostic-D_ORB_YZCet_b_Model_B-pfss-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_M_A.csv')
M_A_BASE = pd.read_csv(FOLDER_A + '/CSV/diagnostic-D_ORB_YZCet_b_Model_A-pfss-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_M_A.csv') #M_A_PFSS_A

#plt.figure(figsize=(8, 7.5))
plt.figure(figsize=(11, 7.5))
lw=2
ax2 = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)
ax2.set_facecolor("white")  

ax2.set_ylim([0, 1])   
ax2.plot(M_A_parker_A[STUDY], M_A_parker_A['M_A'], color='blue', linestyle='solid')
ax2.plot(M_A_PFSS_B[STUDY], M_A_PFSS_B['M_A'], color='purple', linestyle='dotted')
ax2.plot(M_A_BASE[STUDY], M_A_BASE['M_A'], color='orange', linestyle='dashdot')

legend_elements = [
Line2D([0], [0], color='blue', linestyle='solid', lw=lw, label='Model A (Parker spiral)'),
Line2D([0], [0], color='purple', linestyle='dotted', lw=lw, label='Model B (PFSS)'),
Line2D([0], [0], color='orange', linestyle='dashdot', lw=lw, label='PFSS Base (A PFSS)'),
]

ax2.legend(handles=legend_elements, loc='upper right', fontsize=12, facecolor='white', edgecolor='white', framealpha=1)
if STUDY == "D_ORB":
    ax2.set_xlim(left=1,right=60)      
ax2.set_xlabel(xlabel,fontsize=20)
ax2.set_ylabel(r"$M_{\rm A}$")
#ax2.axvspan(x[0], R_SS, facecolor='gray', alpha=0.7)
ax2.axvline(x=R_SS, color='gray', alpha=0.7)
secax = ax2.secondary_yaxis('right', functions=(spi.identity,spi.identity))
ax2.axvline(x = xnom, ls='-.', color='k', lw=2)
ax2.axvline(x = 0.02156*au/R_star, ls='-.', color='k', lw=2)
ax2.axvline(x = 0.02851*au/R_star, ls='-.', color='k', lw=2)
#ax2.text(xnom,1*1.01,'b',ha='center',fontsize=13)
#ax2.text(0.02156*au/R_star,1*1.01,'c',ha='center',fontsize=13)
#ax2.text(0.02851*au/R_star,1*1.01,'d',ha='center',fontsize=13)
ax2.text(xnom-1,1*0.95,'b',ha='center',fontsize=13)
ax2.text(0.02156*au/R_star-1,1*0.95,'c',ha='center',fontsize=13)
ax2.text(0.02851*au/R_star-1,1*0.95,'d',ha='center',fontsize=13)
plt.savefig(path+'/'+'M_A_'+'_model_comparison-Pineda.pdf', bbox_inches='tight')





##comparison plots for v/v_A
#M_A_parker_A = pd.read_csv('/home/luis/github/spirou/OUTPUT/YZCet_b_Model_A/CSV/diagnostic-D_ORB_YZCet_b_Model_A-open-parker-spiral-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_M_A.csv')
#M_A_PFSS_B = pd.read_csv('/home/luis/github/spirou/OUTPUT/YZCet_b_Model_B/CSV/diagnostic-D_ORB_YZCet_b_Model_B-pfss-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_M_A.csv')

M_A_parker_A = pd.read_csv(FOLDER_A + '/CSV/diagnostic-D_ORB_YZCet_b_Model_A-open-parker-spiral-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_M_A.csv')
M_A_PFSS_B = pd.read_csv(FOLDER_B +'/CSV/diagnostic-D_ORB_YZCet_b_Model_B-pfss-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_M_A.csv')

plt.figure(figsize=(8, 7.5))
lw=2
ax2 = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)
ax2.set_facecolor("white")  

#ax2.set_yscale('log')
ax2.plot(M_A_parker_A[STUDY], M_A_parker_A['M_A']**(-1), color='blue', linestyle='dashed')
ax2.plot(M_A_PFSS_B[STUDY], M_A_PFSS_B['M_A']**(-1), color='orange', linestyle='dotted')

legend_elements = [
Line2D([0], [0], color='blue', linestyle='dashed', lw=lw, label='Model A (Parker spiral)'),
Line2D([0], [0], color='orange', linestyle='dotted', lw=lw, label='Model B (PFSS)'),
]

ax2.legend(handles=legend_elements, loc='lower left', fontsize=12, facecolor='white', edgecolor='white', framealpha=0)
if STUDY == "D_ORB":
    ax2.set_xlim(left=1,right=60)      
ax2.set_xlabel(xlabel,fontsize=20)
ax2.set_ylabel(r"$v/v_{\rm A}$")
ax2.axvspan(x[0], R_SS, facecolor='gray', alpha=0.7)
secax = ax2.secondary_yaxis('right', functions=(spi.identity,spi.identity))
ax2.axvline(x = xnom, ls='--', color='k', lw=2)
plt.savefig(path+'/'+'v_A_v_'+'_model_comparison-Pineda.pdf', bbox_inches='tight')
































fig, (ax1,ax2) = plt.subplots(2, 1, sharex=True, figsize=(11, 14))
#fig, ax1 = plt.subplots(1, 1, sharex=True, figsize=(11, 7))
fig.subplots_adjust(hspace=0)
##comparison plots for the stellar wind magnetic field
#bsw_parker_A = pd.read_csv('/home/luis/github/spirou/OUTPUT/YZCet_b_Model_A/CSV/diagnostic-D_ORB_YZCet_b_Model_A-open-parker-spiral-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_B_sw.csv')
#bsw_pfss_B = pd.read_csv('/home/luis/github/spirou/OUTPUT/YZCet_b_Model_B/CSV/diagnostic-D_ORB_YZCet_b_Model_B-pfss-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_B_sw.csv')
#bsw_BASE = pd.read_csv('/home/luis/github/spirou/OUTPUT/YZCet_b_Model_A/CSV/diagnostic-D_ORB_YZCet_b_Model_A-pfss-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_B_sw.csv') #M_A_PFSS_A

bsw_parker_A = pd.read_csv(FOLDER_A+'/CSV/diagnostic-D_ORB_YZCet_b_Model_A-open-parker-spiral-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_B_sw.csv')
bsw_pfss_B = pd.read_csv(FOLDER_B+'/CSV/diagnostic-D_ORB_YZCet_b_Model_B-pfss-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_B_sw.csv')
bsw_BASE = pd.read_csv(FOLDER_A+'/CSV/diagnostic-D_ORB_YZCet_b_Model_A-pfss-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_B_sw.csv') #M_A_PFSS_A

#plt.figure(figsize=(11, 7.5))
lw=2
#ax1 = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)
ax1.set_facecolor("white")  

ax1.set_yscale('log')
ax1.plot(bsw_parker_A[STUDY], bsw_parker_A['Bsw'], color='blue', linestyle='solid')
ax1.plot(bsw_pfss_B[STUDY], bsw_pfss_B['Bsw'], color='purple', linestyle='dotted',lw=5)
ax1.plot(bsw_BASE[STUDY], bsw_BASE['Bsw'], color='orange', linestyle='dashdot',lw=2.9)

print('BSW: Model A / MODEL B: ', bsw_parker_A['Bsw']/bsw_pfss_B['Bsw'])

legend_elements = [
Line2D([0], [0], color='blue', linestyle='solid', lw=lw, label='Model A (Parker spiral)'),
Line2D([0], [0], color='purple', linestyle='dotted', lw=lw, label='Model B (PFSS)'),
Line2D([0], [0], color='orange', linestyle='dashdot', lw=lw, label='PFSS Base (A PFSS)', alpha=0.5),
]

ax1.legend(handles=legend_elements, loc='upper right', fontsize=14, facecolor='white', edgecolor='white', framealpha=1)
if STUDY == "D_ORB":
    ax2.set_xlim(left=1,right=60)      
ax1.set_xlabel(xlabel,fontsize=20)
ax1.set_ylabel(r"$B_{\rm sw}$ $[G]$")
#ax1.axvspan(x[0], R_SS, facecolor='gray', alpha=0.7)
ax1.axvline(x=R_SS, color='gray', alpha=0.7)
secax_ax1 = ax1.secondary_yaxis('right', functions=(spi.identity,spi.identity))
#ax2.axvline(x = xnom, ls='--', color='k', lw=2)
ax1.axvline(x = 0.01634*au/R_star, ls='-.', color='k', lw=2)
ax1.axvline(x = 0.02156*au/R_star, ls='-.', color='k', lw=2)
ax1.axvline(x = 0.02851*au/R_star, ls='-.', color='k', lw=2)

#ax2.text(xnom,5e2,'b',ha='center',fontsize=11)
#ax2.text(0.01634*au/R_star,max(bsw_parker_A['Bsw'])*1.01,'c',ha='center',fontsize=13)
#ax2.text(0.02156*au/R_star,max(bsw_parker_A['Bsw'])*1.01,'c',ha='center',fontsize=13)
#ax2.text(0.02851*au/R_star,max(bsw_parker_A['Bsw'])*1.01,'d',ha='center',fontsize=13)
#ax2.text(0.01634*au/R_star,4e2,'b',ha='center',fontsize=13)
#ax2.text(0.02156*au/R_star,4e2,'c',ha='center',fontsize=13)
#ax2.text(0.02851*au/R_star,4e2,'d',ha='center',fontsize=13)
ax1.text(0.01634*au/R_star-1,2e2,'b',ha='center',fontsize=18)
ax1.text(0.02156*au/R_star-1,2e2,'c',ha='center',fontsize=18)
ax1.text(0.02851*au/R_star-1,2e2,'d',ha='center',fontsize=18)

#plt.savefig('/home/luis/github/spirou/OUTPUT/YZ_Cet_Pineda'+'/'+'B_sw_'+'_model_comparison-Pineda.pdf', bbox_inches='tight')
#plt.savefig(path+'/'+'B_sw_'+'_model_comparison-Pineda.pdf', bbox_inches='tight')
ax1.set_ylim([5e-3,4e2])



##comparison plots for M_A
#M_A_parker_A = pd.read_csv('/home/luis/github/spirou/OUTPUT/YZCet_b_Model_A/CSV/diagnostic-D_ORB_YZCet_b_Model_A-open-parker-spiral-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_M_A.csv')
#M_A_PFSS_B = pd.read_csv('/home/luis/github/spirou/OUTPUT/YZCet_b_Model_B/CSV/diagnostic-D_ORB_YZCet_b_Model_B-pfss-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_M_A.csv')
#M_A_BASE = pd.read_csv('/home/luis/github/spirou/OUTPUT/YZCet_b_Model_A/CSV/diagnostic-D_ORB_YZCet_b_Model_A-pfss-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_M_A.csv') #M_A_PFSS_A


M_A_parker_A = pd.read_csv(FOLDER_A+'/CSV/diagnostic-D_ORB_YZCet_b_Model_A-open-parker-spiral-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_M_A.csv')
M_A_PFSS_B = pd.read_csv(FOLDER_B+'/CSV/diagnostic-D_ORB_YZCet_b_Model_B-pfss-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_M_A.csv')
M_A_BASE = pd.read_csv(FOLDER_A+'/CSV/diagnostic-D_ORB_YZCet_b_Model_A-pfss-Bstar220.0G-Bplanet'+'['+"{:.3f}".format(Bplanet_field)+']'+'G-1.0e-03-1.0e-03-T_corona1.5MKSPI_at_1.0R_star_M_A.csv') #M_A_PFSS_A


#plt.figure(figsize=(8, 7.5))
#plt.figure(figsize=(11, 7.5))
lw=2
#ax2 = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)
ax2.set_facecolor("white")  

ax2.set_ylim([0, 1])   
ax2.plot(M_A_parker_A[STUDY], M_A_parker_A['M_A'], color='blue', linestyle='solid',lw=3)
ax2.plot(M_A_PFSS_B[STUDY], M_A_PFSS_B['M_A'], color='purple', linestyle='dotted',lw=6)
ax2.plot(M_A_BASE[STUDY], M_A_BASE['M_A'], color='orange', linestyle='dashdot',lw=2.9)

legend_elements = [
Line2D([0], [0], color='blue', linestyle='solid', lw=lw, label='Model A (Parker spiral)'),
Line2D([0], [0], color='purple', linestyle='dotted', lw=lw, label='Model B (PFSS)'),
Line2D([0], [0], color='orange', linestyle='dashdot', lw=lw, label='PFSS Base (A PFSS)'),
]

#ax2.legend(handles=legend_elements, loc='upper right', fontsize=12, facecolor='white', edgecolor='white', framealpha=1)
if STUDY == "D_ORB":
    ax2.set_xlim(left=1,right=60)      
ax2.set_xlabel(xlabel,fontsize=20)
ax2.set_ylabel(r"$M_{\rm A}$")
#ax2.axvspan(x[0], R_SS, facecolor='gray', alpha=0.7)
ax2.axvline(x=R_SS, color='gray', alpha=0.7)
secax_ax2 = ax2.secondary_yaxis('right', functions=(spi.identity,spi.identity))
ax2.axvline(x = xnom, ls='-.', color='k', lw=2)
ax2.axvline(x = 0.02156*au/R_star, ls='-.', color='k', lw=2)
ax2.axvline(x = 0.02851*au/R_star, ls='-.', color='k', lw=2)
#ax2.text(xnom,1*1.01,'b',ha='center',fontsize=13)
#ax2.text(0.02156*au/R_star,1*1.01,'c',ha='center',fontsize=13)
#ax2.text(0.02851*au/R_star,1*1.01,'d',ha='center',fontsize=13)
#ax2.text(xnom-1,1*0.95,'b',ha='center',fontsize=13)
#ax2.text(0.02156*au/R_star-1,1*0.95,'c',ha='center',fontsize=13)
#ax2.text(0.02851*au/R_star-1,1*0.95,'d',ha='center',fontsize=13)

#plt.show()
#plt.savefig('/home/luis/github/spirou/OUTPUT/YZ_Cet_Pineda'+'/'+'B_sw_M_A_'+'_model_comparison-Pineda.pdf', bbox_inches='tight')
plt.savefig(path+'/'+'B_sw_M_A_'+'_model_comparison-Pineda.pdf', bbox_inches='tight')
plt.savefig(path, bbox_inches='tight')



#FOLDER = 'OUTPUT/' + str(Exoplanet.replace(" ", "_"))
#path=FOLDER + '/'+'COMPARISON_PDF'+'/'+'M_A_'+'model_comparison-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-hybrid-Bstar'+"{:.1f}".format(B_star) + "G" + "-Bplanet" + '['+"{:.3f}".format(Bplanet_field)+']'+ "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'.pdf'

# Hide ax1
ax2.remove()

# Resize ax2 to fill the figure space
#fig.subplots_adjust(top=0.95, bottom=0.05, left=0.1, right=0.95)
#ax1.set_position([0.1, 0.1, 0.85, 0.85])  # [left, bottom, width, height] in figure coords

# Save the figure with only ax2 visible
plt.savefig(path+'/'+'B_sw_'+'_model_comparison-Pineda.pdf', bbox_inches='tight')
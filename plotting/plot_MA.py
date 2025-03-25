import numpy as np

from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.patches as mpatches
#matplotlib.rc_file_defaults()
plt.style.use(['bmh','SPIworkflow/spi.mplstyle'])


#Individual MA plot
print('##########################101############################')
print(Exoplanet)
print('##########################202############################')
height=30
#if Exoplanet == 'Trappist-1 b Turn-regular':
if Exoplanet == 'Trappist-1 b':
	print('##########################303############################')
	print('##########################404############################')
	print(Exoplanet)
	print('##########################505############################')
	print('##########################606############################')
	fig, ax3 = plt.subplots(1, 1,figsize=(6, 6), sharex=True)
	ax3.plot(x, M_A*np.ones(len(x)), color='k', lw=lw)
	xlabel=r"Orbital separation / Stellar radius"
	ax3.set_xlabel(xlabel,fontsize=15)
	ax3.set_ylabel(r"$M_A$")
	ax3.set_facecolor("white")
	ax3.axvline(x = xnom, ls='--', color='k', lw=2)
	ax3.axhline(y = 1, ls='-.', color='grey', lw=2) 
	#print('geometry', geometry)
	if STUDY == "D_ORB" and Bfield_geom_arr[ind] == 'pfss':
	    ax3.axvspan(x[0], R_SS, facecolor='grey', alpha=0.6)
		

	#b
	ax3.axvline(x = 0.01111*au/R_star, ls='-.', color='k', lw=1.5)      
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
	ax3.text(0.01111*au/R_star,height,'b',ha='center',fontsize=11)
	ax3.text(0.01521*au/R_star,height,'c',ha='center',fontsize=11)
	ax3.text(0.02144*au/R_star,height,'d',ha='center',fontsize=11)
	ax3.text(0.02817*au/R_star,height,'e',ha='center',fontsize=11)
	ax3.text(0.0371*au/R_star,height,'f',ha='center',fontsize=11)
	ax3.text(0.0451*au/R_star,height,'g',ha='center',fontsize=11)
	ax3.text(0.063*au/R_star,height,'h',ha='center',fontsize=11) 

	if Exoplanet=='Proxima b Turnpenney' and STUDY == 'D_ORB':     
	    #d      
	    ax1.axvline(x = 0.02885*au/R_star, ls='-.', color='k', lw=1.5)
	    ax2.axvline(x = 0.02885*au/R_star, ls='-.', color='k', lw=1.5)
	    ax3.axvline(x = 0.02885*au/R_star, ls='-.', color='k', lw=1.5)
	    ax4.axvline(x = 0.02885*au/R_star, ls='-.', color='k', lw=1.5)   
	    ax1.text(0.04856*au/R_star,1.1e5,'b',ha='center',fontsize=11) 
	    ax1.text(0.02885*au/R_star,1.1e5,'d',ha='center',fontsize=11)



	ax3.set_yscale('log')   
	ax3.set_ylim(1e-3,2e1)
	ax3.set_xlim(1,150)
        
	yticks3 = ax3.get_yticks()
	ax3.set_yticks(yticks3[2:-2]) 

	#secax3 = ax3.secondary_yaxis('right', functions=(spi.identity,spi.identity))

	#secax3.set_yticks(yticks3[2:-2]) 

	out_M_A =  FOLDER + '/DIAG_PDF'+ '/'+ "-MA-variation"  + "_" + str(Exoplanet.replace(" ", "_")) +  geometry + diagnostic_string +'.pdf' 
	print('saving: ',out_M_A)
	#plt.show()
	plt.savefig(out_M_A,bbox_inches='tight')

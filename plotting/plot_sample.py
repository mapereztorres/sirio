import os
import shutil
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from scipy.special import lambertw
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
#matplotlib.rc_file_defaults()
plt.style.use(['bmh','SPIworkflow/spi.mplstyle'])

from output import OutputWriter

## Getting parameters to predict SPI radio emission
#
# Observing parameters: Observin frequency, assumed rms noise, Delta_nu
# Stellar parameters: T_corona, n_base_corona, B_field, Isothermality of plasma
# Geometry of the sub-Alfvénic interaction: ALPHA_SPI, THETA_M
#
from setup import *

# Import useful constants and functions to be used in the code
from SPIworkflow.constants import *
import SPIworkflow.SPIutils as spi
import SPIworkflow.freefree as ff
from SPIworkflow.load_data import get_spi_data, create_data_tables, load_target, table_lists

import importlib
import re

geom_list = ['closed_dipole', 'pfss']
colors = {
    'Flux_r_S': 'orange',
    'Flux_reconnect': 'blue',
    'Flux_sb': 'green'
}

def extract_number(s):
    if isinstance(s, str):
        match = re.search(r'(\d+(\.\d+)?)', s)  # find first number in the string
        return float(match.group(1)) if match else np.nan
    return np.nan  # if not a string, return NaN

# Apply to your DataFrame column





ylim_low = 1.001e-4
ylim_high = 2e4

ylim_low_prot = 1
ylim_high_prot = 200

M_star_dot=M_DOT_DEFAULT
#plot_data=pd.read_csv("OUTPUT/flux_data_"+str(geom)+'M_star_dot_'+str(M_star_dot)+".csv", index=False)
protmin_min = 1
P_rot_min = 1
P_rot_max = 200
#flux_data_closed_dipole


originaltable=pd.read_csv('./INPUT/targets_'+TABLE+'.csv')
originaltable["Numeric_SpType"] = originaltable['star_sp_type'].apply(extract_number)
originaltable=originaltable.rename(columns={"planet_name": "name"})

#min_separation = originaltable['a(au)'].min()
#max_separation = originaltable['a(au)'].max()    
min_separation = 3
max_separation = 250

min_star_mass = originaltable['mass_star(m_sun)'].min()*0.8   
max_star_mass = originaltable['mass_star(m_sun)'].max()*1.2   

min_sp = originaltable["Numeric_SpType"].min()-0.5     
max_sp = originaltable["Numeric_SpType"].max()*1.2     

min_dist = originaltable["d_star(pc)"].min()*0.8  
max_dist = originaltable["d_star(pc)"].max()*1.2  
#min_dist = 1
#max_dist = 80

originaltable['cycl_freq']= 2.8*originaltable['bfield_star(gauss)']

min_obs_freq = originaltable['cycl_freq'].min()*0.8
max_obs_freq = originaltable['cycl_freq'].max()*1.2  
#min_obs_freq = 350
#max_obs_freq = 2000


#min_planet_mass= originaltable['mass_planet(m_earth)'].min()*0.8
#max_planet_mass= originaltable['mass_planet(m_earth)'].min()*1.2
min_planet_mass = 0.1
max_planet_mass = 100 

plot_data=pd.read_csv('OUTPUT/flux_data_closed_dipole_M_star_dot_1.0.csv')

latex_table = "\\begin{table}[h!]\n\\centering\n\\begin{tabular}{|c|l|}\n\\hline\n"
latex_table += "No. & Name \\\\\n\\hline\n"

for i, name in enumerate(plot_data['name'], start=1):
    latex_table += f"{i} & {name} \\\\\n"
latex_table += "\\hline\n\\end{tabular}\n\\caption{List of Names}\n\\end{table}"

print(latex_table)    
    

flux_cases = [
    ('Flux_r_S', 'Alfvén wing', 'Alfven_wing'),
    #('Flux_reconnect', 'Reconnection', 'reconnection'),
    ('Flux_sb', 'Stretch and break', 'sb')
]



# Define plotting configurations
plot_vars = [
    {
        "col": "mass_star(m_sun)",
        "xlabel": "Star mass / Sun mass",
        "xscale": "linear",
        "xlim": (min_star_mass, max_star_mass),
        "suffix": "starmass"
    },
    {
        "col": "x",
        "xlabel": "Orbital separation / Stellar radius",
        "xscale": "log",
        "xlim": (min_separation, max_separation),
        "suffix": "separation"
    },
    {
        "col": "Numeric_SpType",
        "xlabel": "Spectral type (M)",
        "xscale": "linear",
        "xlim": (min_sp, max_sp),
        "suffix": "star_sp_type"
    },
    {
        "col": "d_star(pc)",
        "xlabel": "Distance (pc)",
        "xscale": "log",
        "xlim": (min_dist, max_dist),
        "suffix": "distance"
    },
    {
        "col": "mass_planet(m_earth)",
        "xlabel": "Planet mass ($M_\\oplus$)",
        "xscale": "log",
        "xlim": (min_planet_mass,max_planet_mass),
        "suffix": "planetmass"
    },
    {
        "col": "P_rot",
        "xlabel": "$P_{\\rm rot}$ (days)",
        "xscale": "log",
        "xlim": (P_rot_min,P_rot_max),
        "suffix": "P_rot"
    },
    {
        "col": "obs_freq",
        "xlabel": "$\\nu_{\\rm obs}$ (MHz)",
        "xscale": "log",
        "xlim": (min_obs_freq, max_obs_freq),
        "suffix": "freq"
    }
]


for geom in geom_list:
    #plot_data = flux_data[geom]
    #plot_data['Flux_r_S'] = plot_data['Flux_r_S'].apply(lambda x: x[0] if isinstance(x, list) else x)
    #plot_data['Flux_r_S'] = pd.to_numeric(plot_data['Flux_r_S'], errors="coerce")

    if 'flux_data' in locals() or 'flux_data' in globals():
        plot_data = pd.DataFrame(flux_data[geom])  
    else:
        plot_data=pd.read_csv("OUTPUT/flux_data_"+str(geom)+'_M_star_dot_'+str(M_star_dot)+".csv")
        #plot_data=pd.read_csv("OUTPUT/flux_data_"+str(geom)+".csv")
    
    print('Number of interacting planets in the '+geom+' geometry for the Alfvén wing model:', len(plot_data['Flux_r_S'][~np.isnan(plot_data['Flux_r_S'])]))
    print('Number of interacting planets in the '+geom+' geometry for the reconnection model:', len(plot_data['Flux_reconnect'][~np.isnan(plot_data['Flux_reconnect'])]))
    print('Number of interacting planets in the '+geom+' geometry for the Stretch and Break model:', len(plot_data['Flux_sb'][~np.isnan(plot_data['Flux_sb'])]))
    
    
    plot_data['Flux_r_S'] = plot_data['Flux_r_S'].apply(lambda x: x[0] if isinstance(x, list) else x)
    plot_data['Flux_r_S'] = pd.to_numeric(plot_data['Flux_r_S'], errors="coerce")
    plot_data['Flux_r_S'] = plot_data['Flux_r_S']*1000
    plot_data['Flux_reconnect'] = plot_data['Flux_reconnect'] * 1000
    plot_data['Flux_sb'] = plot_data['Flux_sb'] * 1000

    
    #plot_data = plot_data[(plot_data['Flux_r_S'] > 1e-3) & (plot_data['Flux_r_S'] < 1e2)]
    
    #plot_data=plot_data.rename(columns={"name": "planet_name"})
    
    plot_data=plot_data.merge(originaltable,on='name',how='left')
    
    


    for flux_col, flux_label, file_suffix in flux_cases:
        for pv in plot_vars:
            if pv["col"] not in plot_data.columns:
                print(f"Warning: Column {pv['col']} not in data, skipping.")
                continue
                
                
            # Filter points within limits
            mask = (
                pd.notna(plot_data[pv["col"]]) &
                pd.notna(plot_data[flux_col]) &
                np.isfinite(plot_data[pv["col"]]) &
                np.isfinite(plot_data[flux_col]) &
                (plot_data[flux_col] > ylim_low) &
                (plot_data[flux_col] < ylim_high)
            )
            if pv["xlim"] is not None:
                mask &= (
                    (plot_data[pv["col"]] > pv["xlim"][0]) &
                    (plot_data[pv["col"]] < pv["xlim"][1])
                )

            subset = plot_data[mask]

            # Skip empty plots
            if subset.empty:
                print(f"No valid data for {pv['suffix']} in {geom} - {flux_label}")
                continue
                
            fig, ax = plt.subplots(figsize=(10, 12), constrained_layout=True)
            ax.set_yscale("log")
            if pv["xscale"]:
                ax.set_xscale(pv["xscale"])

            ax.set_ylim([ylim_low, ylim_high])
            if pv["xlim"] is not None:
                ax.set_xlim(pv["xlim"])

            # Scatter plot
            ax.scatter(plot_data[pv["col"]], plot_data[flux_col], alpha=0.6)

            # Annotate with numbers
            '''
            number_map = {name: i + 1 for i, name in enumerate(plot_data["name"])}
            for xi, frs, name in zip(plot_data[pv["col"]], plot_data[flux_col], plot_data["name"]):
                if pd.notna(xi) and pd.notna(frs) and np.isfinite(xi) and np.isfinite(frs):
                    if ylim_low < frs < ylim_high:
                        ax.text(xi, frs, str(number_map[name]), fontsize=8, ha="right", va="bottom")
                        
            '''            
            number_map = {name: i + 1 for i, name in enumerate(plot_data["name"])}
            for xi, frs, name in zip(subset[pv["col"]], subset[flux_col], subset["name"]):
                ax.text(xi, frs, str(number_map[name]), fontsize=8, ha="right", va="bottom")

            # Labels and title
            ax.set_xlabel(pv["xlabel"])
            ax.set_ylabel("Flux ($\\mu Jy$)")
            ax.set_title(f"{geom.replace('_', ' ').title()} - {flux_label}")

            # Save file
            #fname = f"OUTPUT/alltargets_flux_{pv['suffix']}_{file_suffix}_{geom}.pdf"
            #fname="OUTPUT/alltargets_flux_data_"+str(geom)+'_M_star_dot_'+str(M_star_dot)+".pdf"
            #fname = f"OUTPUT/"+TABLE+"_alltargets_flux_{pv['suffix']}_{file_suffix}_{geom}"+'_M_star_dot_'+str(M_star_dot)+".pdf" 
            #alltargets_flux_distance_Alfven_wing_closed_dipole
            fname="OUTPUT/alltargets_flux_distance_"+str(geom)+'_M_star_dot_'+str(M_star_dot)+".pdf"
            plt.tight_layout()
            plt.savefig(fname, bbox_inches="tight")
            plt.close()

            
    
    fig, ax = plt.subplots(figsize=(10, 12), constrained_layout=True)
    #ax.set_ylim([1.001e-3, 1e2])
    #ax.set_xlim([10, 50])
    ax.set_yscale('log')    
    #ax.set_xscale('log')    

    ax.set_ylim([ylim_low_prot, ylim_high_prot])
    ax.set_xlim([min_star_mass, max_star_mass])
    star_data = plot_data.drop_duplicates(subset="star_name", keep="first")
    #print(star_data["star_name"])
    #ax.set_xlim([1, 150])
    # Scatter plots

    ax.scatter(star_data['mass_star(m_sun)'], star_data['P_rot'], s=60, alpha=0.6)
    #ax.scatter(plot_data['x'], plot_data['Flux_reconnect'], color=colors['Flux_reconnect'], marker='x', s=60)
    #ax.scatter(plot_data['x'], plot_data['Flux_sb'], color=colors['Flux_sb'], marker='^', s=60)
    # Assign numbers to names

    number_map = {name: i+1 for i, name in enumerate(star_data['star_name'])}
    # Annotate points with numbers
    for xi, frs, name in zip(star_data['mass_star(m_sun)'], star_data['P_rot'], star_data['star_name']):
        if pd.notna(xi) and pd.notna(frs) and np.isfinite(xi) and np.isfinite(frs):
            if ylim_low < frs < ylim_high:
                ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
                #if 1 < xi < 150:
                #    ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
    # Create a legend mapping numbers to names
    mapping_handles = [Line2D([0], [0], color='none', label=f"{num} = {name}") 
                       for name, num in number_map.items()]
    #ax.legend(handles=mapping_handles, loc='outside upper right', fontsize=6, title='Targets')
    #ax.legend(handles=mapping_handles, loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=6, title='Targets')
    # Axis labels & title
    ax.set_xlabel('Stellar Mass / Sun Mass')
    ax.set_ylabel('Rotation period (days)')
    #ax.set_title(geom.replace('_', ' ').title())
    #ax.set_xlabel(xlabel,fontsize=20) 
    plt.tight_layout()
    #plt.savefig('OUTPUT/alltargets_P_rot_star_starmass_' + geom + '.pdf', bbox_inches='tight')
    plt.savefig('OUTPUT/'+TABLE+'_alltargets_P_rot_star_starmass_' + geom +'_M_star_dot_'+str(M_star_dot)+'.pdf', bbox_inches='tight')
    plt.close() 
 
    
    
    fig, ax = plt.subplots(figsize=(10, 12), constrained_layout=True)
    #ax.set_ylim([1.001e-3, 1e2])
    #ax.set_xlim([10, 50])
    #ax.set_yscale('log')    
    #ax.set_xscale('log')    
    #ylim_low = 1
    #ylim_high = 200
    #ax.set_ylim([ylim_low, ylim_high])
    star_data = plot_data.drop_duplicates(subset="star_name", keep="first")
    #print(star_data["star_name"])
    #ax.set_xlim([1, 150])
    # Scatter plots
    #ax.scatter(plot_data['x'], plot_data['Flux_r_S'], color=colors_freq, s=60)
    ax.scatter(star_data['Numeric_SpType'], star_data['mass_star(m_sun)'], s=60, alpha=0.6)
    #ax.scatter(plot_data['x'], plot_data['Flux_reconnect'], color=colors['Flux_reconnect'], marker='x', s=60)
    #ax.scatter(plot_data['x'], plot_data['Flux_sb'], color=colors['Flux_sb'], marker='^', s=60)
    # Assign numbers to names

    number_map = {name: i+1 for i, name in enumerate(star_data['star_name'])}
    # Annotate points with numbers
    for xi, frs, name in zip(star_data['Numeric_SpType'], star_data['mass_star(m_sun)'], star_data['star_name']):
        if pd.notna(xi) and pd.notna(frs) and np.isfinite(xi) and np.isfinite(frs):
            ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
                #if 1 < xi < 150:
                #    ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
    # Create a legend mapping numbers to names
    mapping_handles = [Line2D([0], [0], color='none', label=f"{num} = {name}") 
                       for name, num in number_map.items()]
    #ax.legend(handles=mapping_handles, loc='outside upper right', fontsize=6, title='Targets')
    #ax.legend(handles=mapping_handles, loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=6, title='Targets')
    # Axis labels & title
    ax.set_xlabel('Spectral type (M)')
    ax.set_ylabel('Stellar Mass / Sun Mass')
    #ax.set_title(geom.replace('_', ' ').title())
    #ax.set_xlabel(xlabel,fontsize=20) 
    plt.tight_layout()
    plt.savefig('OUTPUT/alltargets_starmass_spectraltype_' + geom + '.pdf', bbox_inches='tight')
    plt.close()       
  


#plot_data['name']

latex_table = "\\begin{table}[h!]\n\\centering\n\\begin{tabular}{|c|l|}\n\\hline\n"
latex_table += "No. & Name \\\\\n\\hline\n"

for i, name in enumerate(plot_data['name'], start=1):
    latex_table += f"{i} & {name} \\\\\n"
    latex_table += "\\hline\n\\end{tabular}\n\\caption{List of Names}\n\\end{table}"

print(latex_table)    
    

    

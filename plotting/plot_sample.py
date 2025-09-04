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

geom_list = ['open_parker_spiral', 'closed_dipole', 'pfss']
colors = {
    'Flux_r_S': 'orange',
    'Flux_reconnect': 'blue',
    'Flux_sb': 'green'
}

ylim_low = 1.001e-6
ylim_high = 1e-1

M_star_dot=int(M_DOT_DEFAULT)
#plot_data=pd.read_csv("OUTPUT/flux_data_"+str(geom)+'M_star_dot_'+str(M_star_dot)+".csv", index=False)



#flux_data_closed_dipole


originaltable=pd.read_csv("INPUT/table.csv")
originaltable=originaltable.rename(columns={"planet_name": "name"})
    
def extract_number(s):
    if isinstance(s, str):
        match = re.search(r'(\d+(\.\d+)?)', s)  # find first number in the string
        return float(match.group(1)) if match else np.nan
    return np.nan  # if not a string, return NaN

# Apply to your DataFrame column

originaltable["Numeric_SpType"] = originaltable['star_sp_type'].apply(extract_number)

for geom in geom_list:
    #plot_data = flux_data[geom]
    #plot_data['Flux_r_S'] = plot_data['Flux_r_S'].apply(lambda x: x[0] if isinstance(x, list) else x)
    #plot_data['Flux_r_S'] = pd.to_numeric(plot_data['Flux_r_S'], errors="coerce")

    if 'flux_data' in locals() or 'flux_data' in globals():
        plot_data = pd.DataFrame(flux_data[geom])  
    else:
        plot_data=pd.read_csv("OUTPUT/flux_data_"+str(geom)+'_M_star_dot_'+str(M_star_dot)+".csv")
        #plot_data=pd.read_csv("OUTPUT/flux_data_"+str(geom)+".csv")
    
    
    
    
    plot_data['Flux_r_S'] = plot_data['Flux_r_S'].apply(lambda x: x[0] if isinstance(x, list) else x)
    plot_data['Flux_r_S'] = pd.to_numeric(plot_data['Flux_r_S'], errors="coerce")
    #plot_data = plot_data[(plot_data['Flux_r_S'] > 1e-3) & (plot_data['Flux_r_S'] < 1e2)]
    
    #plot_data=plot_data.rename(columns={"name": "planet_name"})
    
    plot_data=plot_data.merge(originaltable,on='name',how='left')
    
    
    
    def get_color(freq):
        if 400 <= freq < 750:
            return 'red'
        elif 750 <= freq < 1000:
            return 'orange'
        elif 1000 <= freq < 1500:
            return 'darkgreen'
        elif 1500 <= freq < 2000:
            return 'blue'
        else:
            return 'blueviolet'

    
    
    
    #plot_data = plot_data[(plot_data['Flux_r_S'] > 1e-3) & (plot_data['Flux_r_S'] < 1e2)]
        
        
        
        
    # Apply color mapping to your data
    colors_freq = [get_color(freq) for freq in plot_data['obs_freq']]            
    d_values = plot_data['d_star(pc)']  
    min_size, max_size = 150, 10  # adjust for visibility
    sizes = (d_values - d_values.min()) / (d_values.max() - d_values.min()) * (max_size - min_size) + min_size



    
    #### Plotting flux vs orbital separation for the Alfvén wing model for all targets     
    fig, ax = plt.subplots(figsize=(10, 12), constrained_layout=True)
    #ax.set_ylim([1.001e-3, 1e2])
    #ax.set_xlim([10, 50])
    ax.set_yscale('log')    
    ax.set_xscale('log')    
    ylim_low = 1.001e-6
    ylim_high = 1e-1
    ax.set_ylim([ylim_low, ylim_high])
    #ax.set_xlim([1, 150])
    # Scatter plots
    #ax.scatter(plot_data['x'], plot_data['Flux_r_S'], color=colors_freq, s=60)
    ax.scatter(plot_data['x'], plot_data['Flux_r_S'], color=colors_freq, s=sizes, alpha=0.6)
    #ax.scatter(plot_data['x'], plot_data['Flux_reconnect'], color=colors['Flux_reconnect'], marker='x', s=60)
    #ax.scatter(plot_data['x'], plot_data['Flux_sb'], color=colors['Flux_sb'], marker='^', s=60)
    # Assign numbers to names
    number_map = {name: i+1 for i, name in enumerate(plot_data['name'])}
    # Annotate points with numbers
    for xi, frs, name in zip(plot_data['x'], plot_data['Flux_r_S'], plot_data['name']):
        if pd.notna(xi) and pd.notna(frs) and np.isfinite(xi) and np.isfinite(frs):
            if ylim_low < frs < ylim_high:
                ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
                #if 1 < xi < 150:
                #    ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
    # Create a legend mapping numbers to names
    mapping_handles = [Line2D([0], [0], color='none', label=f"{num} = {name}") 
                       for name, num in number_map.items()]
    #ax.legend(handles=mapping_handles, loc='outside upper right', fontsize=6, title='Targets')
    ax.legend(handles=mapping_handles, loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=6, title='Targets')
    # Axis labels & title
    ax.set_xlabel('Orbital separation / Stellar radius')
    ax.set_ylabel('Flux (mJy)')
    ax.set_title(geom.replace('_', ' ').title())
    #ax.set_xlabel(xlabel,fontsize=20) 
    plt.tight_layout()
    plt.savefig('OUTPUT/alltargets_flux_separation_Alfven_wing_' + geom + '.pdf', bbox_inches='tight')
    plt.close()   
    
    
    
    
    

    #### Plotting flux vs orbital separation for the reconnection model for all targets     
    fig, ax = plt.subplots(figsize=(10, 12), constrained_layout=True)
    #ax.set_ylim([1.001e-3, 1e2])
    #ax.set_xlim([10, 50])
    ax.set_yscale('log')    
    ax.set_xscale('log')    
    ylim_low = 1.001e-6
    ylim_high = 1e-1
    ax.set_ylim([ylim_low, ylim_high])
    #ax.set_xlim([1, 150])
    # Scatter plots
    #ax.scatter(plot_data['x'], plot_data['Flux_r_S'], color=colors_freq, s=60)
    #ax.scatter(plot_data['x'], plot_data['Flux_r_S'], color=colors_freq, s=60, alpha=0.6)
    ax.scatter(plot_data['x'], plot_data['Flux_reconnect'], color=colors['Flux_reconnect'], marker='x', s=60)
    #ax.scatter(plot_data['x'], plot_data['Flux_sb'], color=colors['Flux_sb'], marker='^', s=60)
    # Assign numbers to names
    number_map = {name: i+1 for i, name in enumerate(plot_data['name'])}
    # Annotate points with numbers
    for xi, frs, name in zip(plot_data['x'], plot_data['Flux_reconnect'], plot_data['name']):
        if pd.notna(xi) and pd.notna(frs) and np.isfinite(xi) and np.isfinite(frs):
            if ylim_low < frs < ylim_high:
                ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
                #if 1 < xi < 150:
                #    ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
    # Create a legend mapping numbers to names
    mapping_handles = [Line2D([0], [0], color='none', label=f"{num} = {name}") 
                       for name, num in number_map.items()]
    #ax.legend(handles=mapping_handles, loc='outside upper right', fontsize=6, title='Targets')
    ax.legend(handles=mapping_handles, loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=6, title='Targets')
    # Axis labels & title
    ax.set_xlabel('Orbital separation / Stellar radius')
    ax.set_ylabel('Flux (mJy)')
    ax.set_title(geom.replace('_', ' ').title())
    #ax.set_xlabel(xlabel,fontsize=20) 
    plt.tight_layout()
    plt.savefig('OUTPUT/alltargets_flux_separation_reconnection_' + geom + '.pdf', bbox_inches='tight')
    plt.close()   
    
    
    



    #### Plotting flux vs orbital separation for the stretch and break model for all targets     
    fig, ax = plt.subplots(figsize=(10, 12), constrained_layout=True)
    #ax.set_ylim([1.001e-3, 1e2])
    #ax.set_xlim([10, 50])
    ax.set_yscale('log')    
    ax.set_xscale('log')    

    ax.set_ylim([ylim_low, ylim_high])
    #ax.set_xlim([1, 150])
    # Scatter plots
    #ax.scatter(plot_data['x'], plot_data['Flux_r_S'], color=colors_freq, s=60)
    #ax.scatter(plot_data['x'], plot_data['Flux_r_S'], color=colors_freq, s=60, alpha=0.6)
    #ax.scatter(plot_data['x'], plot_data['Flux_reconnect'], color=colors['Flux_reconnect'], marker='x', s=60)
    ax.scatter(plot_data['x'], plot_data['Flux_sb'], color=colors['Flux_sb'], marker='^', s=60)
    # Assign numbers to names
    number_map = {name: i+1 for i, name in enumerate(plot_data['name'])}
    # Annotate points with numbers
    for xi, frs, name in zip(plot_data['x'], plot_data['Flux_sb'], plot_data['name']):
        if pd.notna(xi) and pd.notna(frs) and np.isfinite(xi) and np.isfinite(frs):
            if ylim_low < frs < ylim_high:
                ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
                #if 1 < xi < 150:
                #    ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
    # Create a legend mapping numbers to names
    mapping_handles = [Line2D([0], [0], color='none', label=f"{num} = {name}") 
                       for name, num in number_map.items()]
    #ax.legend(handles=mapping_handles, loc='outside upper right', fontsize=6, title='Targets')
    ax.legend(handles=mapping_handles, loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=6, title='Targets')
    # Axis labels & title
    ax.set_xlabel('Orbital separation / Stellar radius')
    ax.set_ylabel('Flux (mJy)')
    ax.set_title(geom.replace('_', ' ').title())
    #ax.set_xlabel(xlabel,fontsize=20) 
    plt.tight_layout()
    plt.savefig('OUTPUT/alltargets_flux_separation_stretchandbreak_' + geom + '.pdf', bbox_inches='tight')
    plt.close()   
    
    
    ########### STAR MASS ###########
    #### Plotting flux vs orbital separation for the Alfvén wing model for all targets     
    fig, ax = plt.subplots(figsize=(10, 12), constrained_layout=True)
    #ax.set_ylim([1.001e-3, 1e2])
    #ax.set_xlim([10, 50])
    ax.set_yscale('log')    
    #ax.set_xscale('log')    
    ylim_low = 1.001e-6
    ylim_high = 1e-1
    ax.set_ylim([ylim_low, ylim_high])
    #ax.set_xlim([1, 150])
    # Scatter plots
    #ax.scatter(plot_data['x'], plot_data['Flux_r_S'], color=colors_freq, s=60)
    ax.scatter(plot_data['mass_star(m_sun)'], plot_data['Flux_r_S'], color=colors_freq, s=sizes, alpha=0.6)
    #ax.scatter(plot_data['x'], plot_data['Flux_reconnect'], color=colors['Flux_reconnect'], marker='x', s=60)
    #ax.scatter(plot_data['x'], plot_data['Flux_sb'], color=colors['Flux_sb'], marker='^', s=60)
    # Assign numbers to names
    number_map = {name: i+1 for i, name in enumerate(plot_data['name'])}
    # Annotate points with numbers
    for xi, frs, name in zip(plot_data['mass_star(m_sun)'], plot_data['Flux_r_S'], plot_data['name']):
        if pd.notna(xi) and pd.notna(frs) and np.isfinite(xi) and np.isfinite(frs):
            if ylim_low < frs < ylim_high:
                ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
                #if 1 < xi < 150:
                #    ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
    # Create a legend mapping numbers to names
    mapping_handles = [Line2D([0], [0], color='none', label=f"{num} = {name}") 
                       for name, num in number_map.items()]
    #ax.legend(handles=mapping_handles, loc='outside upper right', fontsize=6, title='Targets')
    ax.legend(handles=mapping_handles, loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=6, title='Targets')
    # Axis labels & title
    ax.set_xlabel('STELLAR MASS / SUN MASS')
    ax.set_ylabel('Flux (mJy)')
    ax.set_title(geom.replace('_', ' ').title())
    #ax.set_xlabel(xlabel,fontsize=20) 
    plt.tight_layout()
    plt.savefig('OUTPUT/alltargets_flux_starmass_Alfven_wing_' + geom + '.pdf', bbox_inches='tight')
    plt.close()   
    
    
        
    ########### STAR PLANET MASS RATIO ###########
    plot_data['STAR_PLANET_MASS_RATIO']=(plot_data['mass_planet(m_earth)']*M_earth)/(plot_data['mass_star(m_sun)']*M_jup)
    #### Plotting flux vs orbital separation for the Alfvén wing model for all targets     
    fig, ax = plt.subplots(figsize=(10, 12), constrained_layout=True)
    #ax.set_ylim([1.001e-3, 1e2])
    #ax.set_xlim([10, 50])
    ax.set_yscale('log')    
    ax.set_xscale('log')    
    ylim_low = 1.001e-6
    ylim_high = 1e-1
    ax.set_ylim([ylim_low, ylim_high])
    #ax.set_xlim([1, 150])
    # Scatter plots
    #ax.scatter(plot_data['x'], plot_data['Flux_r_S'], color=colors_freq, s=60)
    ax.scatter(plot_data['STAR_PLANET_MASS_RATIO'], plot_data['Flux_r_S'], color=colors_freq, s=sizes, alpha=0.6)
    #ax.scatter(plot_data['x'], plot_data['Flux_reconnect'], color=colors['Flux_reconnect'], marker='x', s=60)
    #ax.scatter(plot_data['x'], plot_data['Flux_sb'], color=colors['Flux_sb'], marker='^', s=60)
    # Assign numbers to names
    number_map = {name: i+1 for i, name in enumerate(plot_data['name'])}
    # Annotate points with numbers
    for xi, frs, name in zip(plot_data['STAR_PLANET_MASS_RATIO'], plot_data['Flux_r_S'], plot_data['name']):
        if pd.notna(xi) and pd.notna(frs) and np.isfinite(xi) and np.isfinite(frs):
            if ylim_low < frs < ylim_high:
                ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
                #if 1 < xi < 150:
                #    ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
    # Create a legend mapping numbers to names
    mapping_handles = [Line2D([0], [0], color='none', label=f"{num} = {name}") 
                       for name, num in number_map.items()]
    #ax.legend(handles=mapping_handles, loc='outside upper right', fontsize=6, title='Targets')
    ax.legend(handles=mapping_handles, loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=6, title='Targets')
    # Axis labels & title
    ax.set_xlabel('STAR PLANET MASS RATIO')
    ax.set_ylabel('Flux (mJy)')
    ax.set_title(geom.replace('_', ' ').title())
    #ax.set_xlabel(xlabel,fontsize=20) 
    plt.tight_layout()
    plt.savefig('OUTPUT/alltargets_flux_starplanetmassratio_Alfven_wing_' + geom + '.pdf', bbox_inches='tight')
    plt.close()   

    
    ########### SPECTRAL TYPE ###########


    #### Plotting flux vs spectral type for the Alfvén wing model for all targets     
    fig, ax = plt.subplots(figsize=(10, 12), constrained_layout=True)
    #ax.set_ylim([1.001e-3, 1e2])
    #ax.set_xlim([10, 50])
    ax.set_yscale('log')    
    #ax.set_xscale('log')    
    ylim_low = 1.001e-6
    ylim_high = 1e-1
    ax.set_ylim([ylim_low, ylim_high])
    #ax.set_xlim([1, 150])
    # Scatter plots
    #ax.scatter(plot_data['x'], plot_data['Flux_r_S'], color=colors_freq, s=60)
    ax.scatter(plot_data['Numeric_SpType'], plot_data['Flux_r_S'], color=colors_freq, s=60, alpha=0.6)
    #ax.scatter(plot_data['x'], plot_data['Flux_reconnect'], color=colors['Flux_reconnect'], marker='x', s=60)
    #ax.scatter(plot_data['x'], plot_data['Flux_sb'], color=colors['Flux_sb'], marker='^', s=60)
    # Assign numbers to names
    number_map = {name: i+1 for i, name in enumerate(plot_data['name'])}
    # Annotate points with numbers
    for xi, frs, name in zip(plot_data['Numeric_SpType'], plot_data['Flux_r_S'], plot_data['name']):
        if pd.notna(xi) and pd.notna(frs) and np.isfinite(xi) and np.isfinite(frs):
            if ylim_low < frs < ylim_high:
                ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
                #if 1 < xi < 150:
                #    ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
    # Create a legend mapping numbers to names
    mapping_handles = [Line2D([0], [0], color='none', label=f"{num} = {name}") 
                       for name, num in number_map.items()]
    #ax.legend(handles=mapping_handles, loc='outside upper right', fontsize=6, title='Targets')
    ax.legend(handles=mapping_handles, loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=6, title='Targets')
    # Axis labels & title
    ax.set_xlabel('Spectral type (M)')
    ax.set_ylabel('Flux (mJy)')
    ax.set_title(geom.replace('_', ' ').title())
    #ax.set_xlabel(xlabel,fontsize=20) 
    plt.tight_layout()
    plt.savefig('OUTPUT/alltargets_flux_star_sp_type_Alfven_wing_' + geom + '.pdf', bbox_inches='tight')
    plt.close()   






    #### Plotting P_rot vs star mass for the Alfvén wing model for all targets     
    fig, ax = plt.subplots(figsize=(10, 12), constrained_layout=True)
    #ax.set_ylim([1.001e-3, 1e2])
    #ax.set_xlim([10, 50])
    ax.set_yscale('log')    
    #ax.set_xscale('log')    
    ylim_low = 1
    ylim_high = 200
    ax.set_ylim([ylim_low, ylim_high])
    star_data = plot_data.drop_duplicates(subset="star_name", keep="first")
    #print(star_data["star_name"])
    #ax.set_xlim([1, 150])
    # Scatter plots
    #ax.scatter(plot_data['x'], plot_data['Flux_r_S'], color=colors_freq, s=60)
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
    ax.legend(handles=mapping_handles, loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=6, title='Targets')
    # Axis labels & title
    ax.set_xlabel('Stellar Mass / Sun Mass')
    ax.set_ylabel('Rotation period (days)')
    #ax.set_title(geom.replace('_', ' ').title())
    #ax.set_xlabel(xlabel,fontsize=20) 
    plt.tight_layout()
    plt.savefig('OUTPUT/alltargets_P_rot_star_starmass_Alfven_wing_' + geom + '.pdf', bbox_inches='tight')
    plt.close() 







    #### Plotting P_rot vs spectral type for the Alfvén wing model for all targets     
    fig, ax = plt.subplots(figsize=(10, 12), constrained_layout=True)
    #ax.set_ylim([1.001e-3, 1e2])
    #ax.set_xlim([10, 50])
    ax.set_yscale('log')    
    #ax.set_xscale('log')    
    ylim_low = 1
    ylim_high = 200
    ax.set_ylim([ylim_low, ylim_high])
    star_data = plot_data.drop_duplicates(subset="star_name", keep="first")
    #print(star_data["star_name"])
    #ax.set_xlim([1, 150])
    # Scatter plots
    #ax.scatter(plot_data['x'], plot_data['Flux_r_S'], color=colors_freq, s=60)
    ax.scatter(star_data['Numeric_SpType'], star_data['P_rot'], s=60, alpha=0.6)
    #ax.scatter(plot_data['x'], plot_data['Flux_reconnect'], color=colors['Flux_reconnect'], marker='x', s=60)
    #ax.scatter(plot_data['x'], plot_data['Flux_sb'], color=colors['Flux_sb'], marker='^', s=60)
    # Assign numbers to names

    number_map = {name: i+1 for i, name in enumerate(star_data['star_name'])}
    # Annotate points with numbers
    for xi, frs, name in zip(star_data['Numeric_SpType'], star_data['P_rot'], star_data['star_name']):
        if pd.notna(xi) and pd.notna(frs) and np.isfinite(xi) and np.isfinite(frs):
            if ylim_low < frs < ylim_high:
                ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
                #if 1 < xi < 150:
                #    ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
    # Create a legend mapping numbers to names
    mapping_handles = [Line2D([0], [0], color='none', label=f"{num} = {name}") 
                       for name, num in number_map.items()]
    #ax.legend(handles=mapping_handles, loc='outside upper right', fontsize=6, title='Targets')
    ax.legend(handles=mapping_handles, loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=6, title='Targets')
    # Axis labels & title
    ax.set_xlabel('Spectral type (M)')
    ax.set_ylabel('Rotation period (days)')
    #ax.set_title(geom.replace('_', ' ').title())
    #ax.set_xlabel(xlabel,fontsize=20) 
    plt.tight_layout()
    plt.savefig('OUTPUT/alltargets_P_rot_star_sp_type_Alfven_wing_' + geom + '.pdf', bbox_inches='tight')
    plt.close() 





    
    ########### DISTANCE ###########
    
    #### Plotting flux vs distance to the system for the Alfvén wing model for all targets     
    fig, ax = plt.subplots(figsize=(10, 12), constrained_layout=True)
    #ax.set_ylim([1.001e-3, 1e2])
    #ax.set_xlim([10, 50])
    ax.set_yscale('log')    
    ax.set_xscale('log')    
    ylim_low = 1.001e-6
    ylim_high = 1e-1
    ax.set_ylim([ylim_low, ylim_high])
    #ax.set_xlim([1, 150])
    # Scatter plots
    #ax.scatter(plot_data['x'], plot_data['Flux_r_S'], color=colors_freq, s=60)
    ax.scatter(plot_data['d_star(pc)'], plot_data['Flux_r_S'], color=colors_freq, s=sizes, alpha=0.6)
    #ax.scatter(plot_data['x'], plot_data['Flux_reconnect'], color=colors['Flux_reconnect'], marker='x', s=60)
    #ax.scatter(plot_data['x'], plot_data['Flux_sb'], color=colors['Flux_sb'], marker='^', s=60)
    # Assign numbers to names
    number_map = {name: i+1 for i, name in enumerate(plot_data['name'])}
    # Annotate points with numbers
    for xi, frs, name in zip(plot_data['d_star(pc)'], plot_data['Flux_r_S'], plot_data['name']):
        if pd.notna(xi) and pd.notna(frs) and np.isfinite(xi) and np.isfinite(frs):
            if ylim_low < frs < ylim_high:
                ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
                #if 1 < xi < 150:
                #    ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
    # Create a legend mapping numbers to names
    mapping_handles = [Line2D([0], [0], color='none', label=f"{num} = {name}") 
                       for name, num in number_map.items()]
    #ax.legend(handles=mapping_handles, loc='outside upper right', fontsize=6, title='Targets')
    ax.legend(handles=mapping_handles, loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=6, title='Targets')
    # Axis labels & title
    ax.set_xlabel('Distance (pc)')
    ax.set_ylabel('Flux (mJy)')
    ax.set_title(geom.replace('_', ' ').title())
    #ax.set_xlabel(xlabel,fontsize=20) 
    plt.tight_layout()
    plt.savefig('OUTPUT/alltargets_flux_distance_Alfven_wing_' + geom + '.pdf', bbox_inches='tight')
    plt.close()   
    


    
    
    
    ########### PLANET MASS ###########
    
    #### Plotting flux vs planet mass to the system for the Alfvén wing model for all targets     
    fig, ax = plt.subplots(figsize=(10, 12), constrained_layout=True)
    #ax.set_ylim([1.001e-3, 1e2])
    #ax.set_xlim([10, 50])
    ax.set_yscale('log')    
    ax.set_xscale('log')    
    ylim_low = 1.001e-6
    ylim_high = 1e-1
    ax.set_ylim([ylim_low, ylim_high])
    #ax.set_xlim([1, 150])
    # Scatter plots
    #ax.scatter(plot_data['x'], plot_data['Flux_r_S'], color=colors_freq, s=60)
    ax.scatter(plot_data['mass_planet(m_earth)'], plot_data['Flux_r_S'], color=colors_freq, s=sizes, alpha=0.6)
    #ax.scatter(plot_data['x'], plot_data['Flux_reconnect'], color=colors['Flux_reconnect'], marker='x', s=60)
    #ax.scatter(plot_data['x'], plot_data['Flux_sb'], color=colors['Flux_sb'], marker='^', s=60)
    # Assign numbers to names
    number_map = {name: i+1 for i, name in enumerate(plot_data['name'])}
    # Annotate points with numbers
    for xi, frs, name in zip(plot_data['mass_planet(m_earth)'], plot_data['Flux_r_S'], plot_data['name']):
        if pd.notna(xi) and pd.notna(frs) and np.isfinite(xi) and np.isfinite(frs):
            if ylim_low < frs < ylim_high:
                ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
                #if 1 < xi < 150:
                #    ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
    # Create a legend mapping numbers to names
    mapping_handles = [Line2D([0], [0], color='none', label=f"{num} = {name}") 
                       for name, num in number_map.items()]
    #ax.legend(handles=mapping_handles, loc='outside upper right', fontsize=6, title='Targets')
    ax.legend(handles=mapping_handles, loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=6, title='Targets')
    # Axis labels & title
    ax.set_xlabel('Planet mass')
    ax.set_ylabel('Flux (mJy)')
    ax.set_title(geom.replace('_', ' ').title())
    #ax.set_xlabel(xlabel,fontsize=20) 
    plt.tight_layout()
    plt.savefig('OUTPUT/alltargets_flux_planetmass_Alfven_wing_' + geom + '.pdf', bbox_inches='tight')
    plt.close()       
    



   
    ########### STELLAR ROTATION ###########
    
    #### Plotting flux vs stellar rotation period for all targets    
    fig, ax = plt.subplots(figsize=(8, 12))
    #ax.set_ylim([1.001e-3, 1e2])
    #ax.set_xlim([10, 50])
    ax.set_yscale('log') 
    ax.set_xscale('log') 
    ax.set_ylim([ylim_low, ylim_high])
    #ax.set_xlim([0, 200])
    # Scatter plots
    ax.scatter(plot_data['P_rot'], plot_data['Flux_r_S'], color=colors['Flux_r_S'], s=60)
    #ax.scatter(plot_data['x'], plot_data['Flux_reconnect'], color=colors['Flux_reconnect'], marker='x', s=60)
    #ax.scatter(plot_data['x'], plot_data['Flux_sb'], color=colors['Flux_sb'], marker='^', s=60)
    # Assign numbers to names
    number_map = {name: i+1 for i, name in enumerate(plot_data['name'])}
    # Annotate points with numbers
    for xi, frs, name in zip(plot_data['P_rot'], plot_data['Flux_r_S'], plot_data['name']):
        if pd.notna(xi) and pd.notna(frs) and np.isfinite(xi) and np.isfinite(frs):
            if ylim_low < frs < ylim_high:
                ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
    # Create a legend mapping numbers to names
    mapping_handles = [Line2D([0], [0], color='none', label=f"{num} = {name}") 
                       for name, num in number_map.items()]
    ax.legend(handles=mapping_handles, loc='upper right', fontsize=6, title='Targets')
    # Axis labels & title
    ax.set_xlabel('Rotation period (days)')
    ax.set_ylabel('Flux (mJy)')
    ax.set_title(geom.replace('_', ' ').title())
    #ax.set_xlabel(xlabel,fontsize=20)
    
    plt.tight_layout()
    plt.savefig('OUTPUT/alltargets_flux_P_rot_' + geom + '.pdf', bbox_inches='tight')
    plt.close()   
        
    
    
    
    
    
    
    ########### OBSERVING FREQUENCY ###########
    
    
    #### Plotting flux vs observing frequency for all targets    
    fig, ax = plt.subplots(figsize=(8, 12))
    #ax.set_ylim([1.001e-3, 1e2])
    #ax.set_xlim([10, 50])
    ax.set_yscale('log') 
    ax.set_xscale('log') 
    ax.set_ylim([ylim_low, ylim_high])
    #ax.set_xlim([400, 4000])
    # Scatter plots
    ax.scatter(plot_data['obs_freq'], plot_data['Flux_r_S'], color=colors['Flux_r_S'], s=60)
    #ax.scatter(plot_data['x'], plot_data['Flux_reconnect'], color=colors['Flux_reconnect'], marker='x', s=60)
    #ax.scatter(plot_data['x'], plot_data['Flux_sb'], color=colors['Flux_sb'], marker='^', s=60)
    # Assign numbers to names
    number_map = {name: i+1 for i, name in enumerate(plot_data['name'])}
    # Annotate points with numbers
    for xi, frs, name in zip(plot_data['obs_freq'], plot_data['Flux_r_S'], plot_data['name']):
        if pd.notna(xi) and pd.notna(frs) and np.isfinite(xi) and np.isfinite(frs):
            if ylim_low < frs < ylim_high:
                ax.text(xi, frs, str(number_map[name]), fontsize=8, ha='right', va='bottom', rotation=0)
    # Create a legend mapping numbers to names
    mapping_handles = [Line2D([0], [0], color='none', label=f"{num} = {name}") 
                       for name, num in number_map.items()]
    ax.legend(handles=mapping_handles, loc='upper right', fontsize=6, title='Targets')
    # Axis labels & title
    ax.set_xlabel('Observing frequency (MHz)')
    ax.set_ylabel('Flux (mJy)')
    ax.set_title(geom.replace('_', ' ').title())
    #ax.set_xlabel(xlabel,fontsize=20)
    
    plt.tight_layout()
    plt.savefig('OUTPUT/alltargets_flux_freq_' + geom + '.pdf', bbox_inches='tight')
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
    ax.legend(handles=mapping_handles, loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=6, title='Targets')
    # Axis labels & title
    ax.set_xlabel('Spectral type (M)')
    ax.set_ylabel('Stellar Mass / Sun Mass')
    #ax.set_title(geom.replace('_', ' ').title())
    #ax.set_xlabel(xlabel,fontsize=20) 
    plt.tight_layout()
    plt.savefig('OUTPUT/alltargets_starmass_spectraltype_' + geom + '.pdf', bbox_inches='tight')
    plt.close() 
    
    
    
    
        
    
        
    '''    
    #### Plotting flux vs orbital separation for all targets     
    fig, ax = plt.subplots(figsize=(8, 12))
    #ax.set_ylim([1.001e-3, 1e2])
    #ax.set_xlim([10, 50])
    ax.set_yscale('log')    
    ax.set_ylim([1.001e-6, 1e-1])
    #ax.set_xlim([1, 150])
    # Scatter plots
    ax.scatter(plot_data['x'], plot_data['Flux_r_S'], color=colors['Flux_r_S'], s=60)
    #ax.scatter(plot_data['x'], plot_data['Flux_reconnect'], color=colors['Flux_reconnect'], marker='x', s=60)
    #ax.scatter(plot_data['x'], plot_data['Flux_sb'], color=colors['Flux_sb'], marker='^', s=60)
    # Assign numbers to names
    number_map = {name: i+1 for i, name in enumerate(plot_data['name'])}
    # Annotate points with numbers
    for xi, frs, name in zip(plot_data['x'], plot_data['Flux_r_S'], plot_data['name']):
        if pd.notna(xi) and pd.notna(frs) and np.isfinite(xi) and np.isfinite(frs):
            if 1e-6 < frs < 1e-1:
                if 1 < xi < 150:
                    ax.text(xi, frs, str(number_map[name]), fontsize=6, ha='right', va='bottom', rotation=0)
                    
    ax.scatter(plot_data['x'], plot_data['Flux_reconnect'], color=colors['Flux_reconnect'], marker='x', s=60)
    for xi, frs, name in zip(plot_data['x'], plot_data['Flux_reconnect'], plot_data['name']):
        if pd.notna(xi) and pd.notna(frs) and np.isfinite(xi) and np.isfinite(frs):
            if 1e-6 < frs < 1e-1:
                if 1 < xi < 150:
                    ax.text(xi, frs, str(number_map[name]), fontsize=6, ha='right', va='bottom', rotation=0)
    
    ax.scatter(plot_data['x'], plot_data['Flux_sb'], color=colors['Flux_sb'], marker='^', s=60)
    for xi, frs, name in zip(plot_data['x'], plot_data['Flux_sb'], plot_data['name']):
        if pd.notna(xi) and pd.notna(frs) and np.isfinite(xi) and np.isfinite(frs):
            if 1e-6 < frs < 1e-1:
                if 1 < xi < 150:
                    ax.text(xi, frs, str(number_map[name]), fontsize=6, ha='right', va='bottom', rotation=0)                    
    # Create a legend mapping numbers to names
    mapping_handles = [Line2D([0], [0], color='none', label=f"{num} = {name}") 
                       for name, num in number_map.items()]
    ax.legend(handles=mapping_handles, loc='upper right', fontsize=6, title='Targets')
    # Axis labels & title
    ax.set_xlabel('Orbital separation / Stellar radius')
    ax.set_ylabel('Flux (mJy)')
    ax.set_title(geom.replace('_', ' ').title())
    #ax.set_xlabel(xlabel,fontsize=20) 
    plt.tight_layout()
    plt.savefig('OUTPUT/test_alltargets_flux_separation_' + geom + '.pdf', bbox_inches='tight')
    plt.close()   
    '''
        
    
    

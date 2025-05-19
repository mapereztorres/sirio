

# SIRIO

**S**tar-planet **I**nteraction and **R**adio-**I**nduced **O**bservations 

`SIRIO` (pronounced `See-ree-oh`) is a public Python code to predict the radio emission from
Sub-Alfvénic star-planet interaction (SPI), via the electron-cyclotron maser
instability mechanism. `SIRIO` is the Spanish name of the famous star `Sirius`.  


SIRIO omputes the expected flux density from SPI (1) as a function of the orbital
distance (i.e., the semi-major axis) to the planet, or (2) as a function of the
exoplanetary magnetic field, or (3) as a function of the stellar wind mass loss rate.
For cases (2) and (3), the orbital distance of the planet to its host star is kept
fixed. 

The current version of the code assumes an isothermal Parker wind, and considers up to three
different geometries of the stellar wind magnetic field: a closed dipolar field, an open Parker 
spiral, and a Potential Field Source Surface (PFSS), which is a hybrid model, having a dipolar topology up to a certain radial separation from the host star, and adopting an open Parker spiral beyond that radius. The effective radius of the planet, i.e., the magnetosphere radius, is obtained as an equilibrium of pressures between the wind and the planet.  The code can also take into account the free-free extinction of the radio emission within the stellar wind, which is particularly relevant for relatively low
stellar wind temperatures, and/or high-values of the stellar wind mass-loss rate, and/or low observing frequencies. 

The code can be run for a single target, or for a whole table of targets, provided by
the user.  


##  Developers

* Miguel Pérez-Torres
* Luis Peña-Moñino

## Running SIRIO

The main script is `sirio.py`. To be sure the code runs without any issues, run it within the `sirio` environment. To create the `sirio` environment, run the following commad:

``` 
conda env create -f sirio.yml 
```

To activate this environment, use 

``` 
$ conda activate sirio 
```

At this point, you simply run 

```
python sirio.py
```

## File structure

The following file structure of spirou should be self-explanatory. 

```
├── INPUT          - folder containing the input target/targets
│   ├── SINGLE-TARGETS - folder to host single targets of your choice
│   ├── table.csv  - A table containing multiple targets
│   └── target.py  - A file containing a single target 
├── LICENSE
├── OUTPUT         - folder containing the output. Creates a sub-folder per target
├── output.py      - module handling the output
├── pics           - folder containing figures and the logo of the SPIROU code
│   ├── earth.jpg
│   ├── earth.png
│   └── spirou-logo.png
├── plotting       - folder containing modules to handle the plots
│   ├── plot_diagnostic_plots.py
│   ├── plot_effective_radius.py
│   └── plot_flux_density.py
├── README.md
├── setup.py       - file to set up the parameters for the runs of spirou.py
├── sirio.py      - main script
├── sirio.yml     - yml file to create a python environment to run the code
└── SPIworkflow
    ├── constants.py  - file containing useful constants
    ├── freefree.py   - Module computing the free-free extinction
    ├── load_data.py  - Module to handle reading the data
    ├── spi.mplstyle  - Style file for the plots
    └── SPIutils.py   - Main module containing many useful functions
```

When you are done, it is good to deactivate the environment:

``` 
$ conda deactivate 
```

In case there are modifications of some dependencies, the file `sirio.yml` will eventually need to be modified, and you will need to update the environment, by running the followin commands.

```
conda activate sirio 
conda env update --name sirio --file sirio.yml --prune
```

If you encounter any issues, please submit an issue and/or let us know via e-mail: 
Miguel Pérez-Torres (torres@iaa.es) and Luis Peña-Moñino (lpm@iaa.es).


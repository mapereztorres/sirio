<img src="pics/spirou-logo.png" alt="Alt text" style="width: 40%; transform: rotate(90deg);">

# SPIROU

**S**tar-**P**lanet **I**nteraction and **R**adio **O**bservations **U**nited

`SPIROU` is a public Python code to predict the radio emission from
Sub-Alfvénic star-planet interaction (SPI), via the electron-cyclotron maser
instability mechanism.  


The code omputes the expected flux density from SPI (1) as a function of the orbital
distance (i.e., the semi-major axis) to the planet, or (2) as a function of the
exoplanetary magnetic field, or (3) as a function of the stellar wind mass loss rate.
For cases (2) and (3), the orbital distance of the planet to its host star is kept
fixed. 

The current version of the code assumes an isothermal Parker wind, and considers  one
(or both) of the following two geometries of the stellar wind magnetic field: a closed
dipolar one and/or an open Parker spiral. The effective radius of the planet, i.e., the
magnetosphere radius, is obtained as an equilibrium of pressures between the wind and
the planet.  The code can also take into account the free-free extinction of the radio
emission within the stellar wind, which is particularly relevant for relatively low
stellar wind temperatures and/or low observing frequencies. 

The code can be run for a single target, or for a whole table of targets, provided by
the user.  


##  Developers

* Miguel Pérez-Torres
* Luis Peña-Moñino

## Running SPIROU

The main script is `spirou.py`. To be sure the code runs without any issues, run it within the `spirou` environment. To create the `spirou` environment, run the following commad:

``` 
conda env create -f spirou.yml 
```

To activate this environment, use 

``` 
$ conda activate spirou 
```

At this point, you simply run 

```
python spirou.py
```
And SPIROU will run.

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
├── spirou.py      - main script
├── spirou.yml     - yml file to create a python environment to run the code
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

In case there are modifications of some dependencies, the file `spirou.yml` will eventually need to be modified, and you will need to update the environment, by running the followin commands.

```
conda activate spirou 
conda env update --name spirou --file spirou.yml --prune
```

If you encounter any issues, please submit an issue and/or let us know via e-mail: 
Miguel Pérez-Torres (torres@iaa.es) and Luis Peña-Moñino (lpm@iaa.es).


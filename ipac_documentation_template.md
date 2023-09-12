About the FDL PyATMOS Dataset
Overview

The Frontier Development Lab (FDL) is a partnership between NASA, ESA, the SETI Institute, and commercial AI partners, formed with the purpose of applying ML technologies to space science and pushing the frontiers of research across a broad range of space-related topics. The program consists of an annual intensive eight-week research sprint during which small teams, with members drawn from the fields of AI, data science, and the space sector, work on a range of challenge areas.

The PyATMOS dataset is the result of the work of one of the teams from the FDL 2018 Astrobiology challenge. It comprises a set of model earth-like exoplanet atmospheres, spanning a range of plausible parameter space.

Software Used to Create the Dataset

ATMOS

ATMOS is a software package produced and maintained by the Virtual Planetary Laboratory which is publicly available online. The ATMOS package is a coupled photochemistry-climate model, designed to simulate "stable" atmospheres given thermodynamic and chemical inputs. The inputs to ATMOS can be either a set of gas concentrations or a set of gas fluxes at the surface of the planet, as well as planetary parameters and stellar parameters (e.g. the gravitational field strength of the planet, and the stellar spectrum). If the inputs are gas concentrations, then the software will find what gas fluxes are required to maintain the given concentrations in a stable atmosphere. Conversely, if the inputs are fluxes, then the software will calculate the gas concentrations of the resulting stable atmosphere. A stable atmosphere is defined as one where the gas concentrations do not change over a sufficiently long timestep—the length of which varies depending on context but can be as long as the age of the universe. Sometimes a "stable" atmosphere cannot be maintained, for example if certain chemical species are broken down or produced by photochemical reactions too quickly. The output of ATMOS is a 1D column of the resultant atmosphere's temperature, pressure, gas concentrations and gas fluxes as a function of altitude.

PyATMOS

PyATMOS (available on GitHub) is the software used to actually run ATMOS hundreds of thousands of times and create a database of exoplanet atmospheres. It is a dockerized wrapper for ATMOS controllable via python, and therefore should be usable on any machine with docker and python, regardless of the operating system.

The Exoplanet Archive PyATMOS Dataset

The NASA Exoplanet Archive hosts around 125,000 simulated exoplanet atmospheres. The astrobiology team of the 2018 NASA FDL programme used the PyATMOS software to create a database of possible stable exoplanet atmospheres. All of these are based around an Earth-like planet that orbits a star similar to the Sun, but with different gas mixtures in their atmospheres. A parameter space of possible atmospheres was scanned by varying the concentrations of the following gasses:

Carbon dioxide (CO2)
Oxygen (O2)
Water vapour (H2O)
Methane (CH4)
Hydrogen (H2)
Nitrogen (N2)
A number of other gasses remain in the atmosphere, such as ozone (O3), but these have small concentrations which are fixed for all planets. The planetary and stellar parameters were not varied.
Using the Interactive Interface to Browse the Dataset

Clicking on the 'PyATMOS Dataset' tab opens an interface for interactively browsing the dataset. The interface is composed of two panels, each comprising a table and an interactive plot.

Summary Table

The upper panel displays a summary table of all the available atmospheric models, with one model per row, and with columns representing the basic input and output parameters for each model. By default, the plot to the right of the table shows a two-dimensional histogram of the distribution of the models across parameter space—initially against surface temperature and pressure.

Model Preview

The bottom panel shows a preview of the model indicated by the currently highlighted row in the summary table (top panel). The table in the bottom panel shows the content of one of the main output files as returned by PyATMOS, parsed_clima_final.csv (see 'Data Products Description' for more information). This table records the variation of several quantities (pressure, temperature, etc.) through each layer of the atmosphere, with one row per layer in the atmosphere. The plot to its right shows a visualisation of the data in the table, by default plotting the curve of temperature against altitude.

Navigating the Dataset

Clicking on any row in the summary table (top panel) highlights that row, and updates the model preview panel (bottom panel) accordingly. The highlighted row (atmospheric layer) in the bottom table is indicated by a similarly-colored point on the curve in the model preview plot to the right. Clicking any row in the lower table udpates the marked point in that plot.

Both of the tables can be sorted by clicking on the column headings, and filtered by entering filter criteria into the text boxes underneath the column headings (e.g., you might enter '<300' in the 'Temperature' column of the summary table). The respective plots will be updated accordingly.

The table and plot interfaces are inherited from a toolkit developed at IPAC to support the IRSA Viewer and Catalog Search tools, and provide a lot more functionality via the icons in the gray bar at the top of each. For more information, refer to the following parts of the IRSA Catalog Search Tool documentation:

Interacting with the Tables
Interacting with the Plots
Downloading data

Downloading a single model
To download the full data for a single model, simply highlight the row corresponding to the model you wish to download in the summary table in the top panel, and then click 'Download This Model' in the bottom panel. This will immediately download a tarred set of all the files associated with the chosen model.

Downloading a selection of models
To download data for a selection of models, indicate each of the models you wish to download by checking their respective checboxes at the left side of the summary table in the top panel. If you have a large number of models to download, judicious application of filters on the columns, followed by use of the 'check-all' check box at the very top-left of the table should save considerable effort. Once you have checked all the required models, click 'Download All Checked Models'. The system will generate a single (possibly large!) wget script which you can run to retrieve all the required models. The script will retrieve a .csv file containing selected rows of the summary table (pyatmos_models.csv), and a single gzipped tar bundle for every model that you have requested. Note that this can result in a large number of files, and a long download time - you may wish to split the wget script in order to manage the download on your filesystem. If you are trying to download a significant fraction of the dataset, you may find it more efficient to consider instead downloading the entire dataset in compressed form (see below).

Downloading the entire dataset
To download the entire dataset comprising all the atmospheric models and the complete summary table, download and run the following wget script:

Download Entire FDL PyATMOS Dataset (compressed: ~40 GB; uncompressed: ~110GB)

Note: that the dataset is large, and download may take considerable time!

For a description of the files delivered with each model, see 'Data Products Description' below.

Data Products Description

There is a set of 20 data files associated with each model atmosphere, organized into a single directory labelled with a unique hash for each atmosphere (identified in the 'Directory Name' column of the summary table).

In each "atmosphere" directory we see the following files:

Clima_log.txt
Photo_log.txt
TempIn.dat
TempOut.dat
clima_allout.tab
in.dist
input_clima.dat *
mixing_ratios.dat *
out.dist
out.out
parsed_clima_final.csv
parsed_clima_final.npy.npz
parsed_clima_initial.csv
parsed_clima_initial.npy.npz
parsed_clima_iterations.csv
parsed_clima_iterations.npy.npz
parsed_photochem_fluxes.csv
parsed_photochem_fluxes.npy.npz
parsed_photochem_mixing_ratios.csv
parsed_photochem_mixing_ratios.npy.npz
run_metadata.json
species.dat
(* Files marked with an asterisk may or may not be present for a given model - see below.)

In addition, if you downloaded multiple models, a single summary table file is included corresponding to the table in the top panel of the interactive tool:

pyatmos_models.csv
Below is a description of each of these files. The products for each model can broadly be separated into two classes, those produced and used directly by the ATMOS software, and those produced by PyATMOS.

pyatmos_models.csv
Single summary table for all the downloaded data, with one model per row, and with columns representing the basic input and output parameters for each model, as well as the hashes used to label the directories for each model. See 'Column Definitions' for a complete description of the columns. Provided only if more than one model is downloaded.

species.dat
This file controls the chemical species that were present during the atmosphere simulation. One of the important fields is "FIXEDMR", which corresponds to the mixing ratio of this particular gas at the surface of the planet. The other important field is "SGFLUX" which corresponds to the gas flux at the surface of the planet. One can only have either a constant mixing ratio or a constant flux, and this is toggled by the LBOUND field. An LBOUND of 1 indicates a constant mixing ratio is being set, a value of 2 indicates a constant flux is being set. Other fields are described in more detail inside each species.dat file.

Clima_log.txt and Photo_log.tex
These two files are simply a recording of the text output from running the "clima" and "photochem" parts of the ATMOS software, respectively.

TempIn.dat and TempOut.dat
These two files are used internally by the climate model. They contain the temperature profile (as a function of altitude, or atmospheric layer number), and also the "FH20" values in the atmosphere at that altitude. This information is repeated in parsed_clima_initial.csv and parsed_clima_final.csv

out.out
Contains lists and tables of many of the results of ATMOS, in particular the photochemical model. These results include the calculated fluxes and mixing ratios of gasses of the stable planetary atmosphere. This file has been parsed for the user by the PyATMOS software to extract the relevant information. The two resulting files are "parsed_photochem_fluxes.csv" and "parsed_photochem_mixing_ratios.csv".

clima_allout.tab
This file stores tabulated results from the climate model within ATMOS, as well as some of the inputs to the climate model. The inputs include the initial state of the atmosphere (temperature, pressure as a function of altitude) as well as some stellar parameters and planetary albedo. The atmosphere's final state (again temperature and pressure as a function of altitude) is also stored. This file has been parsed for the user by the PyATMOS software to extract the relevant information. The three resulting files are "parsed_clima_initial.csv", "parsed_clima_final.csv" and "parsed_clima_iterations.csv".

parsed_clima_initial.csv and parsed_clima_final.csv
These contain the pressure and temperature of the atmosphere as a function of altitude for the initial state of the atmosphere as well as the final state, respectively. The relevant columns are "ALT" (altitude of the atmosphere layer in km), "P" (pressure of the atmosphere layer in bar), "T" (temperature of the atmosphere layer in Kelvin). The other columns are not relevant.

parsed_clima_iterations.csv
A record of the temperature and temperature change between each iteration of the photochemical model computation.

parsed_photochem_fluxes.csv and parsed_photochem_mixing_ratios.csv
This file contains information parsed from the "out.out" file into a more computer-friendly format. The information is the fluxes and mixing ratios of each gas species (labelled by the column headings) at each altitude ("Z") in cm, respectively.

in.dist
The set of parameters used to configure the photochemical mode.

input_clima.dat
This file contains the inputs to the climate model of ATMOS, including the stellar parameters and some of the planetary parameters—for example, the gravitational field strength of the planet (at the planet's surface), or the fraction of solar radiance received (relative to what the earth receives). Mostly the default settings in this file are used, except in cases where the input methane mixing ratio was > 1e-4, in which case, a special flag was set to allow the climate model to consider hazes, which are relevant for atmospheres with large methane concentrations. If this file is not present in the download, then one can assume the default settings were used.

mixing_ratios.dat
Mixing ratios of the relevant gases passed from the photochemical part of the model to the climate model. This file is handled internally by ATMOS and so is not always included in the downloadable files.

out.dist
The set of the parameters produced by the photochemical model. This can be used as a replacement for "in.dist" in further iterations of the photochemical model if a slightly different set of parameters were used.

run_metadata.json
This file stores metadata about the PyATMOS run that was used to generate this atmosphere, including the initial gas concentrations, duration of the calculations, number of iterations required by the photochemical model to find a stable atmosphere.

*npy.npz files
These are duplicates of the csv files that were created by PyATMOS in a compressed numpy format.

Column Definitions

Summary Table Column Definitions

Column Name	Table Label	Units	Description
concentration_CH4	CH4 Concentration	fractional	CH4 concentration at planet surface calculated by model
concentration_CO2	CO2 Concentration	fractional	CO2 concentration at planet surface calculated by model
concentration_H2	H2 Concentration	fractional	H2 concentration at planet surface calculated by model
concentration_H2O	H2O Concentration	fractional	H2O concentration at planet surface calculated by model
concentration_O2	O2 Concentration	fractional	O2 concentration at planet surface calculated by model
flux_CH4	CH4 Flux	molecules/s/cm2	CH4 flux required to maintain gas concentrations at planet surface, calculated by model
flux_CO2	CO2 Flux	molecules/s/cm2	CO2 flux required to maintain gas concentrations at planet surface, calculated by model
flux_H2	H2 Flux	molecules/s/cm2	H2 flux required to maintain gas concentrations at planet surface, calculated by model
flux_H2O	H2O Flux	molecules/s/cm2	H2O flux required to maintain gas concentrations at planet surface, calculated by model
flux_O2	O2 Flux	molecules/s/cm2	O2 flux required to maintain gas concentrations at planet surface, calculated by model
hash	Directory Name	none	Identifier for the directory containing full model data
input_CH4	Input CH4 concentration	fractional	CH4 concentration input to model (at planet surface)
input_CO2	Input CO2 concentration	fractional	CO2 concentration input to model (at planet surface)
input_H2	Input H2 concentration	fractional	H2 concentration input to model (at planet surface)
input_H2O	Input H2O concentration	fractional	H2O concentration input to model (at planet surface)
input_O2	Input O2 concentration	fractional	O2 concentration input to model (at planet surface)
pressure_bar	Pressure	bar	Pressure at planet surface
temperature_kelvin	Temperature	K	Temperature at planet surface
Model Preview Column Definitions

Column Name	Table Label	Units	Description
J	Layer Number	none	Atmospheric layer number
P	Pressure	bar	Atmospheric pressure at this layer
Alt	Altitude	km	Altitude at this layer
T	Temperature	K	Temperature at this layer
CONVEC	Is Layer Convective?	none	Is the atmosphere convective in this layer? (0 = no; 1 = yes)
FH2O	H2O Fraction	none	Fraction of the atmosphere composed of H2O in this layer
FO3	O3 Fraction	none	Fraction of the atmosphere composed of O3 in this layer

The FDL PyATMOS Team

FDL team members
Aaron Bell (University of Tokyo), Aditya Chopra (Australian National University), William Fawcett (University of Cambridge), Rodd Tabeli (Georgia Institute of Technology).

FDL mentors
Daniel Angerhausen (Universität Bern), Anamaria Bera (University of Central Florida), Natalie Cabrol (SETI Institute), Chris Kempes (Santa Fe Institute), Massimo Mascaro (Applied AI, Google).

Last updated: 02 August 2019
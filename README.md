tsetse_climate_change
=================

# <span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">tsetse_climate_change code repository</span>

This repository contains all files needed to run the analyses in:

> <span xmlns:cc="http://creativecommons.org/ns#" property="cc:attributionName">Lord J.S., Hargrove J.W., Torr S.J. and Vale, G.A. (2018). Climate change and African try-
panosomiasis vector populations in Zimbabwe's Zambezi Valley: a mathematical modelling study.
PLoS Medicine. https://doi.org/10.1371/journal.pmed.1002675.


The code is made available under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>. You are free to reuse this code provided that you give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use. Giving appropriate credit includes citation of the above publication and providing a link to this repository:

<a xmlns:dct="http://purl.org/dc/terms/" href="https://github.com/jenniesuz/tsetse_climate_change" rel="dct:source">https://github.com/jenniesuzL/tsetse_climate_change</a>

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />

I suggest you open [this project file](tsetse_climate_change.Rproj) file in [R Studio](rstudio.org). The files that form part of this project file are detailed below. I suggest that you view the files in the following order:
- [**r_1_data_bioassay.R**](r_1_data_bioassay.R)
- [**r_2_subfuncs_adult_mortality.R**](r_2_subfuncs_adult_mortality.R) 
- [**r_2_subfuncs_larviposition.R**](r_2_subfuncs_larviposition.R)
- [**r_2_subfuncs_pupal_duration.R**](r_2_subfuncs_pupal_duration.R)
- [**r_2_subfuncs_pupal_mortality.R**](r_2_subfuncs_pupal_mortality.R)
- [**r_3_fit_funcs.R**](r_3_fit_funcs.R)
- [**r_4_tsetsemod.R**](r_4_tsetsemod.R)


These files are explained briefly below, the data used in each file are also provided here in .csv files.

## Analyses
- **r_1_data_bioassay.R** - reads in tsetse count data and temperature data, calculates monthly average temperatureture
- **r_2_subfuncs_adult_mortality.R** - adult mortality as a function of temperature
- **r_2_subfuncs_larviposition.R** - larviposition as a function of temperature
- **r_2_subfuncs_pupal_duration.R** - pupal duration as a function of temperature
- **r_2_subfuncs_pupal_mortality.R** - pupal mortality as a function of temperature
- **r_3_fit_funcs.R** - functions for fitting the model
- **r_4_tsetsemod.R** - tsetse population dynamics model

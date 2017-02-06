# Impacts of representing sea-level rise uncertainty on future flood risks: An example from San Francisco Bay

## Overview
This code uses "*a case study to quantify and illustrate how neglecting sea-level rise uncertainties can bias risk projections.*" The case study focuses on the future 1-in-100 year storm surge in the year 2100 in the San Francisco Bay area. Accounting for the uncertainties associated with future sea-level rise is important because it "*can have considerable impacts on the design of flood risk management strategies.*"

This R code and the [GIS tutorial](https://download.scrim.psu.edu/Ruckert_etal_San_Francisco/) is intended to help users who wish to work with the storm surge and sea-level projections/ hindcasts in greater detail than provided in the text. Key functionality of these scripts include:

1. Project changes in sea level from 1880 to 2300 with associated probabilities
2. Estimate the current 100-year storm surge for the San Francisco Bay area
3. Estimate potential future storm surge in 2100 from 3 methods of accounting for sea level change:
  * no change in sea-level
  * adding the mean sea-level rise estimate
  * accounting for sea-level uncertainty (i.e., the entire distribution)
4. Determine the area at risk to potential future coastal flooding in the San Francisco Bay area
5. Produce figures from the paper 


Additionally, a "shiny" application (https://clima.shinyapps.io/ROK\_SFB\_2016/) was developed by Kelsey Ruckert to showcase the results in Ruckert et al. (in review). The objective of the Shiny app is to provide an interactive tutorial illustrating how changes in sea level affect flooding events and their associated probability of occurrence, especially over time.

### Citation:
>Ruckert, KL, Oddo, PC, and Keller, K. Impacts of representing sea-level rise uncertainty on future flood risks: An example from San Francisco Bay, *in Review*.

## Requirements
### R
The scripts are written in R (tested under R v3.2.1;  https://www.r-project.org/) using the following packages:
>mcmc  
coda  
RColorBrewer  
DEoptim  
compiler  
extRemes  
fExtremes  
ismev  
lubridate  
zoo  
shiny  
shinyRGL  
plotrix  
ggvis  

You can install and open the packages in R as shown in the example below, which installs the mcmc package.

```R
install.packages("mcmc")
library(mcmc)
```
##### Main contents
* `San_francisco_tides.R`: Runs GEV analysis on tide gauge data from San Francisco
* `Project_global_sealevel.R`: Calibrates and projects global sea-level rise using the Rahmstorf (2007) model and Markov chain Monte Carlo calibration
* `Converge_test_MCMC_SLR.R`: Tests MCMC convergence using the potential scale reduction factor
* `SLR_StormSurge_100yrFlood.R`: Estimates potential future 100-yr flood for the San Francisco Bay area based on the current 100-yr storm surge plus various SLR projections in the year 2100
* `Distribution_test.R`: Estimates potential future 100-yr flood for the San Francisco Bay area based on the current 100-yr storm surge plus various distribution shapes of SLR projections in the year 2100
* `Plot_Ruckertetal_SFB.R`, `Plot_distributions.R`, & `Plot_hypospectic.R`: Generates figures in the paper

### ArcGIS
The GIS analysis is conducted using ArcGIS ArcMAP (tested under version 10.3.1 with the Spatial Analyst extension). **All files needed to conduct the analysis can be found at: https://download.scrim.psu.edu/Ruckert\_etal\_San\_Francisco/.**

##### Main contents
* `Flood_map_documentation.pdf`: Instructions on creating a flood map for the San Francisco Bay area
* `Python.py`: Python codes that accompany the instructions to be used in ArcGIS ArcMAP


## Instructions
* Download the files included
* Make a `Figures`, `Trace`, and `Workspace` directory for R output
* Open and source the `.R` files in R or RStudio in the following order; 1) `San_francisco_tides.R`, 2) `Project_global_sealevel.R`, 3) `Converge_test_MCMC_SLR.R`, 4) `SLR_StormSurge_100yrFlood.R`, 5) `Distribution_test.R`, 6) `Plot_Ruckertetal_SFB.R`, 7) `Plot_distributions.R`, & 8) `Plot_hypospectic.R`.
* Download the GIS files from: https://download.scrim.psu.edu/Ruckert\_etal\_San\_Francisco/
* Follow the ArcGIS instructions in `Flood_map_documentation.pdf` and use the `Python.py` script for additional help

## Contact
Kelsey Ruckert  
E-mail: <klr324@psu.edu>

Corresponding author:  
Klaus Keller  
E-mail: <klaus@psu.edu>
 

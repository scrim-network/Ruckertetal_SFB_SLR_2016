#Current and future flood risk analysis of the San Francisco Bay area from Ruckert et al. (in prep)

README file last updated by Kelsey Ruckert, klr324-at-psu-dot-edu, Fri July 8 12:19:51 EST 2016

## Citation

This code is intended to accompany the results of 

>Ruckert, KL, Oddo, PC, and Keller, K. Accounting for sea-level rise uncertainty increases flood risk area: An example from San Francisco Bay, (in prep).

### Other credits:
>Ruckert, KL, Guan, Y, Forest, FE, and Keller, K. Improving the statistical method can raise the upper tail of sea-level projections, (in prep.)

##Overview

This code requires R and ArcGIS ArcMAP with the following libraries:
- coda
- mcmc
- RColorBrewer
- DEoptim
- compiler

This R code and the GIS tutorial is intended to help users who wish to work with the storm surge and sea-level projections or hindcasts or flood risk analysis in greater detail than provided in the text. Key functionality of these scripts include:

1. Project changes in sea level from 1880 to 2300 with associated probabilities
2. Estimate the current 100-year storm surge for the San Francisco Bay area
3. Estimate potential future storm surge from 3 methods of accounting for sea level change:
  1. no change in sea-level
  2. adding the mean sea-level rise estimate
  3. accounting for the sea-level rise probability distribution
4. Determine potential future flood risk in the San Francisco Bay area
5. Produce figures from the paper

The RFILES directory contains all the scripts and data necessary to run the analysis along with a README file. The prerun analysis output used to generate the Ruckert et al. (in prep.) figures exceeds 100 MB. For access to the prerun analysis please contact the corresponding author. _(Note that the folder directory MUST be in the same format as when downloaded otherwise the scripts will not locate the files/scripts needed to run. Additionally, create a 'Workspace' folder that is empty before running the analysis. Output will be saved to this folder.)_

The most important functions are **Project_global_sealevel** for MCMC calibration and projections, **San_francisco_tides** for Generalized Extreme Value analysis of the San Francisco hourly tide gauge data, and **SLR_StormSurge_100yrFlood** for generating potential future flood frequency curves based on storm surge and sea-level rise.

To run the R analysis, simply open a terminal and run **Project_global_sealevel** (~30 mins.). Then source **San_francisco_tides** (~5 mins), **SLR_StormSurge_100yrFlood** (~8 hrs), and **PlotRuckertetal_SanFranStormSurSLR** (~1hr) to generate plots. (Specific details can be found in the README.txt file.

To run the GIS analysis, follow the instructions provided in the Flood_Analysis_tutorial.pdf files.

##Contact
Kelsey Ruckert: <klr324@psu.edu>  
Corresponding author: Klaus Keller at <klaus@psu.edu>
 

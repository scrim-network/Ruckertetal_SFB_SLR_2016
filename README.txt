  Current and Future 100-yr Storm Surge in San Francisco Bay Codes

  Reference:
  Ruckert, KL, Oddo, PC, Keller, K. Accounting for sea-level rise
  uncertainty increases flood risk area: An example from San Francisco
  Bay, (in prep).

  This program is distributed in the hope that it will be useful,
  but WITH NO WARRANTY (NEITHER EXPLICITNOR IMPLICIT). We are not liable
  for the behavior of these codes in your own application. You are free
  to share this code so long as the authors(s) and version history remain
  intact.

    Kelsey Ruckert, klr324@psu.edu
    Perry Oddo, poddo@psu.edu

 ============================================================
| For the impatient:(NOTE:This program runs for 5+ hours)    |
|  1. Start R                                                |
|  2. Type source(“Project_global_sealevel.R”)               |
|  3. Type source(“San_francisco_tides.R”)                   |
|  4. Type source(“SLR_StormSurge_100yrFlood.R”)             |
|  5. Type source(“PlotRuckertetal_SanFranStormSurSLR.R”)    |
|  6. Open one of the .mxd files in the ROK_analysis folder  |
 ============================================================

Required software:
    R
    ArcGIS ArcMAP

 =============================================================================
| Short description of main R scripts:
|  1. Project_global_sealevel.R: Markov Chain Monte Carlo calibration of the 
|     Rahmstorf SLR model. 
|  2. San_francisco_tides.R: Generalized Extreme Value analysis of the hourly 
|     tide gauge data.
|  3. SLR_StormSurge_100yrFlood.R: Generates flood frequency curves from storm 
|     surge and SLR.
|  4. PlotRuckertetal_SanFranStormSurSLR.R: Plots figures and outputs return 
|     level & period values from the paper.
 =============================================================================


 =============================================================================
| Important NOTE:
|  ~  To run the codes without error the folder must have the following structure: 
|  Main folder
|    Project_global_sealevel.R
|    San_francisco_tides.R
|    SLR_StormSurge_100yrFlood.R
|    PlotRuckertetal_SanFranStormSurSLR.R
|  
|    1. Data
|    2. SLR_scripts
|    3. Figures
|      3.1  SuppFigures
|    4. Workspace
 =============================================================================


 =============================================================================
| Short description of GIS folders:
|  1. Analysis: Contains most if not all of the shapefiles used for flood 
|     risk analysis. NOTE: Rasters and DEMs not included because of size.
|     The .mxd files are ArcMap files set to create the figures in the paper. 
|  2. Tutorial: Contains a step by step tutorial on how flood risk was analzyed
|      using Digital Elevation models.
|     The data needed for this tutorial are included in the folder.
 =============================================================================


Credits:
SLR_scripts: 
    Codes used in:
    Ruckert, KL, Guan, Y, Forest, FE, and Keller, K. Improving the statistical
    method can raise the upper tail of sea-level projections, (submitted to ERL).


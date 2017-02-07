########################################################################################
#
#  -file = "Project_global_sealevel.R"   Code written April 2014, updated Dec 2015
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This program loads in temperature and global sea-level data for use in 
#       projecting a global sea-level model and Markov chain Monte Carlo calibration described
#       in Ruckert et al. (2016; 2017). For further description and references, please read the paper.
#
#  -TIME: This script takes 1 hour to run
#
#  -NOTE: This program applies the Markov chain Monte Carlo method assuming
#       heteroskedastic errors to the Rahmstorf (2007) semi-empirical sea-level model:
#
#       -preserving the AR(1) structure of the observations
#       -creating random realizations of "natural variability"
#       -projecting new realizations to 2300
#
#  -NOTE: This file contains data that is sourced into the other programs. Information
#       regarding this data can be found below:
#
#       -RCP8.5 is used to create temperature simulations to 2300
#       -RCP8.5 simulates "Business as usual" and is similar to the A2 scenario
#       -The RCP8.5 temperatures are from the CNRM-CM5 model
#           (Centre National de Recherches Meteorologiques)
#       -These RCP8.5 temperatures closly resemble historical temperatures
#           from the National Oceanic and Atmospheric Administration (NOAA)
#
#       -Annual global land and ocean temperature anomalies (C)
#       -Anomalies with respect to the 20th century average
#       -http://www.ncdc.noaa.gov/cag/time-series/global
#
#       -Tide guage data from Church & White_GRL_2006
#       -http://www.psmsl.org/products/reconstructions/church.php
#       -SLR values are in relation to the 1990 mean value
#
##==============================================================================
## Copyright 2016 Kelsey Ruckert
## This file is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This file is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this file.  If not, see <http://www.gnu.org/licenses/>.
##==============================================================================
########################################################################################
# Clear global environment
rm(list =ls())

# Open packages:
library(DEoptim)
library(compiler)
enableJIT(3)
enableJIT(3)

# Set the seed
set.seed(111)

#------------------------------ Step 1: Read in Temperature and Sea Level Data ------------------------
# Temperature data
data = read.csv("Data/NOAA_IPCC_RCPtempsscenarios.csv")

## Historical time frame and temperatures from NOAA
# hist.temp, historical global mean temperatures from 1880 to 2300 in C
# alltime, time in years from 1880 to 2300 (1yr increments)
hist.temp = data[1:122,2] 
alltime = data[,1] 

## IPCC temperature scenarios added to NOAA temperatures from 1880-2100
# scenario_time, , time in years from 1880 to 2100 (5yr increments)
# max to b2, merged historical + IPCC temperature scenarios from 1880-2100 in C (5yr increments)
scenario_time = data[1:45,5] 
max = data[1:45,6]  #6.195 C in 2100
min = data[1:45,7]  #1.773 C in 2100
a1fi = data[1:45,8] #4.892 C in 2100
a1b = data[1:45,9]  #3.347 C in 2100
a1t = data[1:45,10] #2.938 C in 2100
a2 = data[1:45,11]  #4.194 C in 2100
b1 = data[1:45,12]  #2.382 C in 2100
b2 = data[1:45,13]  #3.087 C in 2100

## RCP 8.5 temperatures from 2014-2300 and NOAA historical temps from 1880-2013
# rcp85, merged historical + rcp 8.5 temperatures from 1880 to 2300 in C
rcp85 = data[,4]    #4.284 C in 2100 and 9.849284058 C in 2300

# Historical global mean sea levels from tide gauges & estimated errors
# from the Church and White (2006) sea-level file.
# year, time in years from 1880 to 2001, 1yr increments
# slr, global mean sea level in mm
# err.obs, global mean sea level measurement error in mm
# Divide vectors by 10 to convert to cm
church = read.table("Data/church_13221.txt")
year = church[11:132, 1] #timeframe 1880 to 2002 in 1yr incriments
slr = church[11:132, 2]/10 #mm to cm
err.obs = church[11:132,3]/10 #mm to cm

## To match Heberger et al. (2009), sea-level values are set in reference to the 2000 mean SLR value
SLR2000 = slr[120]  # 120 equals the year 2000 
slr = slr - SLR2000

# Set the observational errors by adding and subtracting the errors to the sea-level values
err_pos = slr+err.obs
err_neg = slr-err.obs

# Set the hindcast and projection length
hindcast_length = 122 # there are 122 years from 1880 to 2002
projection_length = 421 # from 1880 to 2300

#----------------------------- Step 2: Calibrate and Project Sea-level Rise ------------------------
#------------------------------ Find Initial Parameter & Initial Hindcast ------------------------
# Physical model parameters [Original parameter estimates in Rahmstorf 2007]:
# [1] alpha =  .34      sensitivity of SLR to temperature change (cm/year/C)
# [2] T_0   = -0.5      baseline temp at which the sea level anomaly is zero (C)
# [3] H_0   = -15       initial sea-level anomaly (cm)

timestep=1 # timesteps are 1 year
from=2 # start from the second year since the first year is an uncertain parameter
to=hindcast_length #122

# Run differential evolution optimization to find initial starting values for
# parameters to use in optimization of the likelihood function.
source("SLR_scripts/Deoptim_rahm_model.R") # physical model
source("SLR_scripts/minimize_residuals.R") # function to minimize the residuals

lower=c(0, -3, err_neg[1])
upper=c(2,  2, err_pos[1])
iter=1000  # specify number of iterations
outDEoptim <- DEoptim(min_res, lower, upper, 
                      DEoptim.control(itermax=iter,
                                      trace=FALSE))
print(outDEoptim$optim$bestmem)# find best initial parameters
parms = c(outDEoptim$optim$bestmem[1], outDEoptim$optim$bestmem[2], outDEoptim$optim$bestmem[3])

# Run the model with the initial parameters to create a simulation of the observations
# Load the physical sea-level model converted to R from the equations in Rahmstorf (2007).
source("SLR_scripts/sealevel_rahm_model.R")

# Use the optimized parameters to generate a fit to the data.
slr.est = rahmfunction(parms, hist.temp)
to = projection_length  # number of years in the projection
proj.est = rahmfunction(parms, rcp85)
to = hindcast_length

#------------------------ Calculate the Residuals & AR(1) Coefficient --------------------------
# Calculate residuals from the fit to the data. Modification of Equation (2).
res = slr-slr.est$sle
nyears.obs = length(year) # number of years in observational time series

# Apply the auto-correlation function to determine a starting value for rho,
# the correlation coefficient.
#pdf(file="ARheter1.pdf", family="Helvetica", pointsize=11, height=4.5, width=4.5)
rho=rep(NA,3)
ac=acf(res, lag.max=5, plot=FALSE, main="")
rho[1]=ac$acf[1]
rho[2]=ac$acf[2]
rho[3]=ac$acf[3]
rho[4]=ac$acf[4]
rho[5]=ac$acf[5]
#dev.off()

#--------------------------------- Run MCMC Calibration ----------------------------------------
# Set up priors.
bound.lower = c(0, -3, err_neg[1], 0, -0.99)
bound.upper = c(2,  2, err_pos[1], 1,  0.99)

# Set the measurement errors for heteroskedastic assumption.
y.meas.err = err.obs 

# Name the model parameters and specify the number of model parameters.
# Sigma and rho are statistical parameters and are not counted in the number.
model.p=3
parnames=c("alpha","base temp","initialvalue", "sigma.y", “rho.y”)

# Load the likelihood model assuming correlated residuals.
source("SLR_scripts/Obs_likelihood_AR.R")

# Optimize the likelihood function to estimate initial starting values.
p = c(outDEoptim$optim$bestmem[1], outDEoptim$optim$bestmem[2], 
      outDEoptim$optim$bestmem[3], sd(res), rho[2]) 
p0 = c(0.34, -0.5, slr[1], 0.6, 0.5) # Rahmstorf estimated best guess parameters
p0 = optim(p0, function(p) -log.post(p))$par
print(round(p0,4))
library(mcmc)

# Set up the step size, burnin, and number of iterations.
step = c(0.02, 0.02, 0.1, 0.01, 0.01)
NI = 2.5e7 # number of iterations
burnin = seq(1,0.01*NI,1) # 1% burnin

# Run the MCMC calibration
mcmc.out1 = metrop(log.post, p0, nbatch=NI, scale=step)
prechain1 = mcmc.out1$batch

# Print the acceptance rate as a percent. Should be ~ 25%
acceptrate = mcmc.out1$accept * 100
cat("Accept rate =", acceptrate, "%\n")

#-------------------------- Estimate sea-level rise projection & Median Estimates ----------------------------
# Thin the chain to a subset; ~20,000 is sufficient.
# Take a subset of every 1237th number in the data set from the burnin to the 25 millionth run
ssprechain1 = prechain1[seq(length(burnin),NI,1237),]
h = length(ssprechain1[,1]) #length of the subset ~20,000

# Find median of the estimated parameters.
median.slr.par = c(median(prechain1[-burnin,1]), median(prechain1[-burnin,2]), median(prechain1[-burnin,3]))

# Estimate model projection from parameter medians.
to = projection_length #421 
median.slr.proj = rahmfunction(median.slr.par, rcp85) 

#------------- Hindcast Sea-level Rise with Uncertainty & Find Plausible Parameters ------------------
# Calculate all possible hindcasts from the subset parameter estimates.
to = hindcast_length

# Set up empty matrices for sea level rate and sea level output.
new.rate = mat.or.vec(h, nyears.obs) # hindcast sea-level rate
mcmc.fit = mat.or.vec(h, nyears.obs) # hindcast SLR simulation
par = mat.or.vec(h, 2) # MCMC chains of parameters 1 and 2

# Load model that estimates the sea level rate of change: Equation (1).
source("SLR_scripts/searate_func.R") #searate function finds the hindcast rates for SLR
for(i in 1:h) {
  par[i,1] = ssprechain1[i,1] # alpha parameter
  par[i,2] = ssprechain1[i,2] # T0 parameter
  
  # Estimate the sea level rate of change: Equation (1).
  new.rate[i,] = thermal(par[i,], hist.temp)[[1]]
  mcmc.fit[i,1] = ssprechain1[i,3]  # Initial value
  
  # Use Forward Euler to estimate sea level over time.
  for (n in from:to){
    mcmc.fit[i,n] = mcmc.fit[i,n-1]+new.rate[i,n-1]*timestep
  }
}

# Estimate the residuals with the AR(1) coefficient (ssprechain1[n,5]) and sigma (ssprechain1[n,4]).
res.mcmc_hind = mat.or.vec(h, nyears.obs) #(nr,nc)
for(n in 1:h) {
  for(i in 2:nyears.obs) {
    # Equation (3-4).
    res.mcmc_hind[n,i] = ssprechain1[n,5]*res.mcmc_hind[n,i-1] + 
      rnorm(1, mean = 0, sd = ssprechain1[n,4]) # add in the AR(1) noise
  }
}
# Estimate the hindcasts: add the residuals onto the model simulations, Equation (2) & (3).
slr.mcmc_hind = res.mcmc_hind
for(i in 1:h) {
  slr.mcmc_hind[i,] = mcmc.fit[i,]+res.mcmc_hind[i,]
}

#----------------------------- Project Sea-level Rise with Uncertainty --------------------------------
### Project SLR with uncertainty using the parameters generated from MCMC assuming heteroskedastic
### errors and RCP8.5 temp. emission
years.mod = (alltime) # all time represent the years from 1880 to 2300
nyears.mod = length(years.mod)
to = projection_length #421

# Set up empty matrices for sea level rate and sea level output.
pred.rate = mat.or.vec(h, nyears.mod) # project sea-level rate
fit.mcmc_proj = mat.or.vec(h, nyears.mod) # project SLR simulation

# Loop over the sea level model to generate a distribution of sea level rates
# and sea level simulations.
for(n in 1:h) {
  # Estimate the sea level rate of change: Equation (1).
  pred.rate[n,] = thermal(par[n,], rcp85)[[1]]
  fit.mcmc_proj[n,1] = ssprechain1[n,3] # initial value in 1880
  
  # Use Forward Euler to estimate sea level over time.
  for (i in from:to){
    fit.mcmc_proj[n,i] = fit.mcmc_proj[n,i-1]+pred.rate[n,i-1]*timestep
  }
}

# Estimate the residuals with the AR(1) coefficient (ssprechain1[n,5]) and sigma (ssprechain1[n,4]).
res.mcmc_proj = mat.or.vec(h, nyears.mod) #(nr,nc)
for(n in 1:h) {
  for(i in 2:nyears.mod) {
    # Equation (3-4).
    res.mcmc_proj[n,i] = ssprechain1[n,5]*res.mcmc_proj[n,i-1] + 
      rnorm(1, mean = 0, sd = ssprechain1[n,4])
  }
}
# Estimate the hindcasts: add the residuals onto the model simulations, Equation (2) & (3).
slr.mcmc_proj = res.mcmc_proj
for(i in 1:h) {
  slr.mcmc_proj[i,] = fit.mcmc_proj[i,]+res.mcmc_proj[i,]
}

#--------------------- Step 3: Recreate the Rahmstorf predictions: ----------------
# Set up the parameters specified in Rahmstorf 2007:
# [1] alpha =  .34      sensitivity of SLR to temperature change (cm/year/C)
# [2] T_0   = -0.5      baseline temp at which the sea level anomaly is zero (C)
# [3] H_0   = -15       initial sea-level anomaly (cm)

original = c(0.34, -0.5, slr[1]) # slr[1] is equal to -0.1464 meters

# Run the hindcast in one year increments from 1880 to 2002
from = 2 # Start from the 2nd year since the first value is defined as an uncertain parameter
to = 122 # Number of observations from 1880 to 2002
timestep = 1 # Data is in increments of one year
hindcast = rahmfunction(original, hist.temp)  # hist.temp are the historical temperatures

# Run the projections using the emission scenarios used in Rahmstorf 2007
to = 45     # There are 45 points in this projection from 1880 to 2100
timestep = 5 # Projections are in intervals of 5 years

# The names of the projections correspond to the names of the emission scenario
max_p = rahmfunction(original, max)
min_p = rahmfunction(original, min)
a1fi_p = rahmfunction(original, a1fi)
a1b_p = rahmfunction(original, a1b)
a1t_p = rahmfunction(original, a1t)
a2_p = rahmfunction(original, a2)
b1_p = rahmfunction(original, b1)
b2_p = rahmfunction(original, b2)

# Make with respect to the year 2000:
a2new = (a2_p$sle/100) - (a2_p$sle[25]/100)

#----------- Step 4: Calculate additional runoff from dams & reservoirs: ----------------
# Heberger 2009: 1.4m (1.38m) from the A2 estimate using rahmstorf 2007 model plus 
# additional accounting for changing runoff from dams & reservoirs
DamsReservoirs = 1.38 - a2new[45]
DamsReservoirs <- round(DamsReservoirs,2) 

# Set up a vector for the sea-level anomaly distribution in 2100.
# Estimate the probability density function of sea-level anomalies in 2100
# WITHOUT accounting for dams and reservoirs.
# The year 2100 is the 221 number in the sequence.
prob_proj2100=mat.or.vec(h,1)
prob_proj2100=slr.mcmc_proj[,221]
pdf2100 <- density(prob_proj2100/100)

# Set up a vector for the sea-level anomaly distribution in 2100.
# Estimate the probability density function of sea-level anomalies in 2100
# ACCOUNTING for dams and reservoirs.
# The year 2100 is the 221 number in the sequence.
proj2100DamsRes <- mat.or.vec(h,1)
proj2100DamsRes <- (slr.mcmc_proj[,221]/100) + DamsReservoirs 
pdf2100DamsRes <- density(proj2100DamsRes)

# plot(pdf2100DamsRes, col="red")
# lines(pdf2100)
# abline(v=1.4)

save.image(file = "Workspace/SFB_globalMeanSeaLevel2.RData")
######################################## END ###############################################
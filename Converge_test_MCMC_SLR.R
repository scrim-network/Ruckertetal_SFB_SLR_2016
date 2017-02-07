########################################################################################
#
#  -file = "Converge_test_MCMC_SLR.R"   Code written April 2014, updated Dec 2015
#  -Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This program runs MCMC calibration on the Rahmstorf (2007) to global sea-level data and
#       tests for markov chain convergence using the potential scale reduction factor.
#
#  -TIME: This script takes 1 hour to run
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
library(coda)
library(DEoptim)
library(compiler)
enableJIT(3)
enableJIT(3)

# Set the seed
set.seed(1780)

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

# Historical global mean sea-levels from tide gauges & estimated errors
# from the Church and White (2006) sea-level file.
# year, time in years from 1880 to 2001, 1yr increments
# slr, global mean sea level in mm
# err.obs, global mean sea level measurement error in mm
# Divide vectors by 10 to convert to cm
church = read.table("Data/church_13221.txt")
year = church[11:132, 1] #timeframe 1880 to 2002 in 1yre incriments
slr = church[11:132, 2]/10 #mm to cm
err.obs = church[11:132,3]/10 #mm to cm

## To match Heberger et al. (2009), sea-level values are set in reference to the 2000 mean SLR value
SLR2000 = slr[120]  # 120 equals the year 2000
slr = slr - SLR2000

# Set the observational errors by adding and subtracting the errors to the sea-level values
err_pos=slr+err.obs
err_neg=slr-err.obs

# Set the hindcast and projection length
hindcast_length=122 # there are 122 years from 1880 to 2002
projection_length=421 # from 1880 to 2300

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
print(outDEoptim$optim$bestmem) # print best initial parameters
parms = c(outDEoptim$optim$bestmem[1], outDEoptim$optim$bestmem[2], outDEoptim$optim$bestmem[3])

# Run the model with the initial parameters to create a simulation of the observations
# Load the physical sea-level model converted to R from the equaitons in Rahmstorf (2007).
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
parnames=c("alpha","base temp","initialvalue", "sigma.y", "rho.y")

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

# Run MCMC calibration.
mcmc.out1780 = metrop(log.post, p0, nbatch=NI, scale=step)
prechain1780 = mcmc.out1780$batch

# Print the acceptance rate as a percent. Should be ~ 25%
acceptrate = mcmc.out1780$accept * 100
cat("Accept rate =", acceptrate, "%\n")

#-------------------------- Convergence Test ----------------------------
# Read in MCMC data with different seed (seed 111)
load("Workspace/SFB_globalMeanSeaLevel2.RData")

# Potential Scale Reduction Factor
# MCMC is converged when the potental scale reduction factor is less than 1.1
heter = as.mcmc(prechain1)
heter2 = as.mcmc(prechain1780)
heterlist = mcmc.list(list(heter, heter2))
gelman.diag(heterlist)

# Create trace plots:
# Remove burnin
results1 = prechain1[length(burnin):NI, ]
results1780 = prechain1780[length(burnin):NI, ]

plot_names = c("01", "02", "03", "04", "05")

# Trace plots for MCMC with seed 111
for(i in 1:5){
    jpeg(file = paste("Trace/MCMC1_", plot_names[i], ".png", sep=""), family="Helvetica", width=1200, height=700, units="px", pointsize=20)
    par(mfrow=c(2,1), mar=c(3, 7, 1, 7), mgp=c(1.5,.5,0))
    plot(results1[ ,i], type="l",
    ylab = parnames[i],
    xlab = "Number of Runs", main="")
    
    par(mar = c(3,7,1,7))
    hist(results1[ ,i], freq=FALSE, col="gray",
    xlab = parnames[i],
    ylab = "Density [Dimensionless]", main="")
    dev.off()
}

# Trace plots for MCMC with seed 1780
for(i in 1:5){
    jpeg(file = paste("Trace/MCMC2_", plot_names[i], ".png", sep=""), family="Helvetica", width=1200, height=700, units="px", pointsize=20)
    par(mfrow=c(2,1), mar=c(3, 7, 1, 7), mgp=c(1.5,.5,0))
    plot(results1780[ ,i], type="l",
    ylab = parnames[i],
    xlab = "Number of Runs", main="")
    
    par(mar = c(3,7,1,7))
    hist(results1780[ ,i], freq=FALSE, col="gray",
    xlab = parnames[i],
    ylab = "Density [Dimensionless]", main="")
    dev.off()
}
######################################## END ###############################################
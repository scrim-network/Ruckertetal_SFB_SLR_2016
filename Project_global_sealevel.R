########################################################################################
#
#  -file = "Project_global_sealevel.R"   Code written April 2014, updated Dec 2015
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This program loads in temperature and global sea-level data for use in 
#       projecting a global sea-level model and calibration described in Ruckert et al. (in prep). For
#       further description and references, please read the paper.
#
#  -This program calibrates Rahmstorf 2007 sea-level model with Markov Chain Monte Carlo
#		accounting for autocorrelated residuals and heteroskedastic errors.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOR IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
#   -TIME: This script takes 1 hour to run
#   -NOTE: This file contains data that is sourced into the other programs. Information
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
########################################################################################
#rm(list =ls()) #Clear all previous variables
library(DEoptim)
library(compiler)
enableJIT(3)
enableJIT(3)

# Set the seed
set.seed(111)  #seed #1

#------------------------------ Step 1: Read in Temperature and Sea Level Data ------------------------
# Temperature data
data = read.csv("Data/NOAA_IPCC_RCPtempsscenarios.csv")

#Historical time frame and temperatures from NOAA
hist.temp = data[1:122,2] #temperature data
alltime = data[,1] #1880-2300 (1 yr incriments)

#IPCC temperature scenarios added to NOAA temperatures from 1880-2100
scenario_time = data[1:45,5] #1880-2100 (5yr incriments)
max = data[1:45,6]  #6.195 C in 2100
min = data[1:45,7]  #1.773 C in 2100
a1fi = data[1:45,8] #4.892 C in 2100
a1b = data[1:45,9]  #3.347 C in 2100
a1t = data[1:45,10] #2.938 C in 2100
a2 = data[1:45,11]  #4.194 C in 2100
b1 = data[1:45,12]  #2.382 C in 2100
b2 = data[1:45,13]  #3.087 C in 2100

#RCP 8.5 temperatures from 2014-2300 and NOAA historical temps from 1880-2013
rcp85 = data[,4]    #4.284 C in 2100 and 9.849284058 C in 2300

# Historical global mean sea-levels from tide gauges & estimated errors
church = read.table("Data/church_13221.txt")
year = church[11:132, 1] #timeframe 1880 to 2002 in 1yre incriments
slr = church[11:132, 2]/10 #mm to cm
err.obs = church[11:132,3]/10 #mm to cm

## To match Heberger et al. (2009), sea-level values are set in reference to the 2000 mean SLR value
SLR2000 = slr[120]  #120 equals the year 2000 
slr = slr - SLR2000

# Set the observational errors by adding and subtracting the errors to the sea-level values
err_pos=slr+err.obs
err_neg=slr-err.obs

# Set the hindcast and projection length
hindcast_length=122 # there are 122 years from 1880 to 2002
projection_length=421 # from 1880 to 2300

#----------------------------- Step 2: Calibrate and Project Sea-level Rise ------------------------
#------------------------------ Find Initial Parameter & Initial Hindcast ------------------------
#parm = c(.34, -0.5) [Original parameter estimates in Rahmstorf 2007]
#[1] 0.34 cm/year per C the sensitivity of SLR to temperature change
#[2] -0.5 #baseline temp (C) at which sea level anomaly is zero
timestep=1 # timesteps are 1 year
from=2 # start from the second year since the first year is an uncertain parameter
to=hindcast_length #122

#Run DEoptim in R to find good initial parameters
source("SLR_scripts/Deoptim_rahm_model.R") # source the model
source("SLR_scripts/minimize_residuals.R") # find the minimum residuals
lower=c(0,-3,err_neg[1])
upper=c(2,2,err_pos[1])
iter=1000  # specify number of iterations
outDEoptim <- DEoptim(min_res, lower, upper, 
                      DEoptim.control(itermax=iter,
                                      trace=FALSE))
print(outDEoptim$optim$bestmem)# find best initial parameters
parms = c(outDEoptim$optim$bestmem[1], outDEoptim$optim$bestmem[2], outDEoptim$optim$bestmem[3])

#Run the model with the initial parameters to create a simulation of the observations
source("SLR_scripts/sealevel_rahm_model.R") #sealevel_rahm_model.R is the model equation
slr.est = rahmfunction(parms, hist.temp)
to=projection_length  #number of years in the projection
proj.est = rahmfunction(parms, rcp85)
to=hindcast_length

#------------------------ Calculate the Residuals & AR(1) Coefficient --------------------------
#Calculate Residuals
res=slr-slr.est$sle
nyears.obs=length(year) #number of years in observational time series

### Estimate and save the lag-1 autocorrelation coefficient (rho[2])
#pdf(file="ARheter1.pdf", family="Helvetica", pointsize=11, height=4.5, width=4.5)
rho=rep(NA,3)
ac=acf(res, lag.max=5, plot=TRUE, main="")# apply  auto-correlation to determine correlation coefficients
rho[1]=ac$acf[1]
rho[2]=ac$acf[2]
rho[3]=ac$acf[3]
rho[4]=ac$acf[4]
rho[5]=ac$acf[5]
#dev.off()

#--------------------------------- Run MCMC Calibration ----------------------------------------
# Set up the prior ranges for the MCMC
bound.lower = c(0,-3,err_neg[1],0,-0.99)
bound.upper = c(2,2,err_pos[1],1,0.99)
y.meas.err = err.obs # measurement error, this changes over time making it heteroskedastic

# Define number of model parameters
model.p=3
parnames=c("alpha","base temp","initialvalue", "sigma.y", "phi11")

# Source the physical model and statistical model
source("SLR_scripts/sealevel_rahm_model.R")
source("SLR_scripts/Obs_likelihood_AR.R")

# Set up the initial parameters from the DEoptim best guess parameters
p = c(outDEoptim$optim$bestmem[1], outDEoptim$optim$bestmem[2], 
      outDEoptim$optim$bestmem[3], sd(res), rho[2]) 
p0 = c(0.34,-0.5,slr[1],0.6, 0.5) # Rahmstorf estimated best guess parameters
p0 = optim(p0, function(p) -log.post(p))$par
print(round(p0,4))
library(mcmc)

# Set up the step size, burnin, and number of iterations to run
step = c(0.02,0.02,0.1,0.01,0.01)
NI = 2.5e7 #number of iterations
burnin = seq(1,0.01*NI,1) # 1% burnin

#Run the MCMC chain
mcmc.out1 = metrop(log.post, p0, nbatch=NI, scale=step)
prechain1 = mcmc.out1$batch
mcmc.out1$accept
# Calculate the parameter acceptance rate
acceptrate = mcmc.out1$accept * 100
#Print the acceptance rate as a percent. Should be ~ 25%
cat("Accept rate =", acceptrate, "%\n")

# Estimating with all 25 million runs is not neccasary if the chains have
# converged. So we will take a subset of every 1237th number in the data set from
# the burnin to the 25 millionth run
ssprechain1 = prechain1[seq(length(burnin),NI,1237),]
h = length(ssprechain1[,1]) #length of the subset ~20,000

#------------- Estimate the median parameters and sea-level rise projection based on the posterior distribution
# Estimate median parameters
median.slr.par = c(median(prechain1[-burnin,1]), median(prechain1[-burnin,2]), median(prechain1[-burnin,3]))
to = projection_length #421 
median.slr.proj = rahmfunction(median.slr.par, rcp85) 

#------------- Hindcast Sea-level Rise with Uncertainty & Find Plausible Parameters ------------------
# Calculate all possible hindcasts from the subset parameter estimates.
to=hindcast_length
new.rate=mat.or.vec(h, nyears.obs)
mcmc.fit=mat.or.vec(h, nyears.obs) # mcmc.fit is the hindcast SLR simulation without noise
par=mat.or.vec(h, 2)
source("SLR_scripts/searate_func.R") #searate function finds the hindcast rates for SLR
for(i in 1:h) {
  to=hindcast_length
  par[i,1]=ssprechain1[i,1] # alpha parameter
  par[i,2]=ssprechain1[i,2] # T0 parameter
  new.rate[i,] = thermal(par[i,], hist.temp)[[1]]
  mcmc.fit[i,1] = ssprechain1[i,3]  # Initial value
  for (n in from:to){
    mcmc.fit[i,n]=mcmc.fit[i,n-1]+new.rate[i,n-1]*timestep
  }
}

### Calculate hindcast residuals with the lag-1 autocorrelation coefficient estimates: ssprechain1[n,5]
### and the standard deviation (sigma) estimates: ssprechain1[n,4]
res.mcmc_hind=mat.or.vec(h, nyears.obs) #(nr,nc)
slr.mcmc_hind=res.mcmc_hind
for(n in 1:h) {
  for(i in 2:nyears.obs) {
    res.mcmc_hind[n,i] = ssprechain1[n,5]*res.mcmc_hind[n,i-1] + 
      rnorm(1,mean=0,sd=ssprechain1[n,4]) # add in the AR(1) noise
  }
}
### superimpose residuals on the hindcasts from the SLR model ###
for(i in 1:h) {
  slr.mcmc_hind[i,]=mcmc.fit[i,]+res.mcmc_hind[i,]
}

#----------------------------- Project Sea-level Rise with Uncertainty --------------------------------
### Project SLR with uncertainty using the parameters generated from MCMC assuming heteroskedastic
### errors and RCP8.5 temp. emission
years.mod=(alltime) # all time represent the years from 1880 to 2300
nyears.mod=length(years.mod)
pred.rate=mat.or.vec(h, nyears.mod) #(nr,nc)
fit.mcmc_proj=mat.or.vec(h, nyears.mod) #(nr,nc) # fit.mcmc_proj is SLR projections without noise
for(n in 1:h) {
  to=projection_length #421
  source("SLR_scripts/searate_func.R") #searate function finds the projected rates for SLR to 2300
  pred.rate[n,] = thermal(par[n,], rcp85)[[1]]
  fit.mcmc_proj[n,1] = ssprechain1[n,3] # initial value in 1880
  for (i in from:to){
    fit.mcmc_proj[n,i]=fit.mcmc_proj[n,i-1]+pred.rate[n,i-1]*timestep
  }
}

###Calculate projection residuals with the lag-1 autocorrelation coefficient estimates: ssprechain1[n,5]
### and the standard deviation (sigma) estimates: ssprechain1[n,4]
res.mcmc_proj=mat.or.vec(h, nyears.mod) #(nr,nc)
slr.mcmc_proj=res.mcmc_proj
for(n in 1:h) {
  for(i in 2:nyears.mod) {
    res.mcmc_proj[n,i] = ssprechain1[n,5]*res.mcmc_proj[n,i-1] + 
      rnorm(1,mean=0,sd=ssprechain1[n,4])
  }
}
### superimpose residuals on the projections from the SLR model ###
for(i in 1:h) {
  slr.mcmc_proj[i,]=fit.mcmc_proj[i,]+res.mcmc_proj[i,]
}

#--------------------- Step 3: Recreate the Rahmstorf predictions: ----------------
# Set up the parameters specified in Rahmstorf 2007:
# The parameters in the model include (in order of appreciance in the vector)
# a = sensitivity of sea-level to changes in temperature. Given as meters/Celsius/year
# T0 = base temperature or temperature when sea-level is zero. Given as Celius
# Initial = the value of sea-level in the year 1880. 1880 is the starting year for this
#       analysis. Given as meters.

original = c(0.34, -0.5, slr[1]) #slr[1] is equal to -0.1464 meters

# Run the hindcast in one year incriments from 1880 to 2002
from = 2 # Start from the 2nd year since the first value is defined as a uncertain parameter
to = 122 # Number of observations from 1880 to 2002
timestep = 1 # Data is in criments of one year
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

# Find the probability density function of sea-level estimates in 2100
# WithOUT accounting for dams and reservoirs
prob_proj2100=mat.or.vec(h,1)
prob_proj2100=slr.mcmc_proj[,221] #The year 2100 is the 221 number in the sequence
pdf2100 <- density(prob_proj2100/100)

# Accounting for dams and reservoirs
proj2100DamsRes <- mat.or.vec(h,1)
proj2100DamsRes <- (slr.mcmc_proj[,221]/100) + DamsReservoirs #The year 2100 is the 221st number
pdf2100DamsRes <- density(proj2100DamsRes)

# plot(pdf2100DamsRes, col="red")
# lines(pdf2100)
# abline(v=1.4)

save.image(file = "Workspace/SFB_globalMeanSeaLevel.RData")
######################################## END ###############################################
#########################################################################
# file: SLR_StormSurge_100yrFlood.R
#------------------------------------------------------------------------
# Author and copyright: Kelsey Ruckert
# Pennsylvania State University
# klr324@psu.edu
# Code written Dec. 2015
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOR IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# NOTE: This script takes 5+ hours to run
#------------------------------------------------------------------------
# Estimates potential future 100-yr flood for the San Francisco Bay area
# based on the current 100-yr storm surge (San_francisco_tides.R) plus
# various SLR projections in the year 2100 (Project_global_sealevel.R)
#
# The various SLR projections include our mean estimate, the estimate
# in Heberger et al. 2009, and accounting for the SLR uncertainty in 2100.
#
# San Francisco, CA
#########################################################################

rm(list =ls()) #Clear global environment
library(extRemes)
library(fExtremes)
library(ismev)
library(lubridate)
library(zoo)
library(compiler)
enableJIT(3)
enableJIT(3)

#------------------- Source in scripts -------------------------------------
source("Scripts/inv_sf.r")

#------------------- Load in the San Francisco tide & sea-level rise workspaces 
# San Francisco tides
load("Workspace/SFB_stormSurge.RData")

# Global mean sea-level
load("Workspace/SFB_globalMeanSeaLevel.RData")

############################# STORM SURGE ESTIMATION ##################################
#------------------- Estimate current storm surge -------------------------------------
storm_surgeL <- 1e5 # desired storm surge level
q = seq(0,1,length.out= storm_surgeL +1)  # quantile array

# Find closed-form solution of GEV fit
fit_q_year = qgev(q, year.res.max.fit2@fit$par.ests[1], year.res.max.fit2@fit$par.ests[2],
                  year.res.max.fit2@fit$par.ests[3])
fit_q_year = fit_q_year[fit_q_year< max(fit_q_year)]

# Find which q is the 100-yr value 
num = which(q <= 0.99)
num.max = which.max(num)
year100prob <- num.max +1
# check to make sure the year100prob value is the 100-yr value (10^-2)
print(1-q[year100prob])

# We know the first value in the GEV fit [107cm] has a 100% probability when we account for any
# amount of future sea-level rise. For finding the inverse of the survival function any values
# that due not exist will be returned as an NA. To return a value we wil add 107cm will a 100%
# probability to the vectors accounting for sea-level rise.
# sea_levelq <- rep(NA, length(q)+1)
# sea_levelq[1] <- 0
# sea_levelq[2:length(sea_levelq)] <- q
# 
# num.sl = which(sea_levelq <= 0.99)
# num.max.sl = which.max(num.sl)
# year100prob.sl <- num.max.sl # -1
# # check to make sure the year100prob.sl value is the 100-yr value (10^-2)
# print(1-q[year100prob.sl])

#------------------- Storm surge plus Heberger et al. 2009 --------------------
#Estimate the storm surge + SLR in heberger et al. 2009 (1.38 m)
#heberger09.2100 <- fit_q_year/100 + 1.38

# If the first the value is -inf, set the first value to 107 (the second value)
# which has a probability of 100%.
print(1-q[1]); print(fit_q_year[1])

# Add 107 cm with a probability of 100% to the vector
heberger09.2100 <- rep(NA, length(fit_q_year))
heberger09.2100[1] = fit_q_year[2]/100
heberger09.2100[2:length(heberger09.2100)] = (fit_q_year[2:length(fit_q_year)]/100) + 1.38

#------------------- Storm surge plus mean sea-level rise ---------------------
#Estimate the mean storm surge + SLR based on the mean of the SLR uncertainty from MCMC
# heteroskedastic calibration
# Mean (expected) return level 
#mean.stormSLR.2100HET <- (fit_q_year/100) + mean(proj2100DamsRes)
#mean.stormSLR.2100HET <- apply(slr.storm.DR2100, 1, mean)

# Add 107 cm with a probability of 100% to the vector
mean.stormSLR.2100HET <- rep(NA, length(fit_q_year))
mean.stormSLR.2100HET[1] = fit_q_year[2]/100
mean.stormSLR.2100HET[2:length(mean.stormSLR.2100HET)] = (fit_q_year[2:length(fit_q_year)]/100) + mean(proj2100DamsRes)

#------------------- Storm surge plus potential sea-level rises ---------------
#Estimate all the potential storm surge + SLRs based on MCMC heteroskedastic calibration
#slr.storm.DR2100 <- mat.or.vec(storm_surgeL,length(proj2100DamsRes))
#for(i in 1:length(proj2100DamsRes)){
#  slr.storm.DR2100[,i] <- (fit_q_year/100) + proj2100DamsRes[i]
#}

slr.storm.DR2100 <- mat.or.vec(length(fit_q_year),length(proj2100DamsRes))
for(i in 1:length(proj2100DamsRes)){
    slr.storm.DR2100[1,i] <- fit_q_year[2]/100
    slr.storm.DR2100[2:length(fit_q_year), i] <- (fit_q_year[2:length(fit_q_year)]/100) + proj2100DamsRes[i]
}

##################### ESTIMATE NEW PROBABILITY OF OCCURRENCE #####################
#------------------- Estimating new Heberger et al. 2009 probabilities --------
# Storm surge plus Heberger et al. 2009
# What is the probability of the 100-yr storm surge associated with Heberger et al. 2009
# when you account for all the potential [states of the world] rises in sea-level by the year 2100?
heb.probs <- apply(slr.storm.DR2100, 2, inv.sf, val=heberger09.2100[year100prob])

# Check if probabilities are estimated for each state of the world
# If so then the minimum SLR plus the maximum storm surge should be equal to or larger than
# the 100-yr Heberger 09 value
print(heberger09.2100[year100prob])
print(fit_q_year[100000]/100 + min(proj2100DamsRes))

# If False then NAs will be returned for values with probabilities smaller than 1: 100,000
# Set all NAs to 0 since they are smaller than 1:100,000. This gives us a conservative estimate,
# but the difference/impact of doing so is minor. This effect was tested in the plotting script.
# See for example supplementary figure 2.
write.csv(heb.probs,"Data/NewHebergerProbabilities.csv",na="0")
heb.probs.table <- read.csv("Data/NewHebergerProbabilities.csv")
new.heb.probs = heb.probs.table[,2]

#Find the new probability for the Heberger 09 value:
new.average.prob = mean(new.heb.probs)

remove(heb.probs)
#------------------- Estimating new Mean probabilities -----------------------
# Storm surge plus mean SLR
# What is the probability of the 100-yr storm surge associated with mean SLR
# when you account for all the potential [states of the world] rises in sea-level by the year 2100?
mean.slr.probs <- apply(slr.storm.DR2100, 2, inv.sf, val=mean.stormSLR.2100HET[year100prob])

# Check if probabilities are estimated for each state of the world
# If so then the minimum SLR plus the maximum storm surge should be equal to or larger than
# the 100-yr mean SLR + storm surge value
print(mean.stormSLR.2100HET[year100prob])
print(fit_q_year[storm_surgeL]/100 + min(proj2100DamsRes))

# If False then NAs will be returned for values with probabilities smaller than 1: 100,000
# Set all NAs to 0 since they are smaller than 1:100,000. This gives us a conservative estimate,
# but the difference/impact of doing so is minor. This effect was tested in the plotting script.
# See for example supplementary figure 2.
write.csv(mean.slr.probs,"Data/NewMeanProbabilities.csv",na="0")
mean.probs.table <- read.csv("Data/NewMeanProbabilities.csv")
new.mean.probs = mean.probs.table[,2]

#Find the new probability for the mean value:
new.meanMSLR.prob = mean(new.mean.probs)

remove(mean.slr.probs)
#------------------- Estimating new current storm surge probabilities ----------
# Storm surge
# What is the probability of the 100-yr storm surge associated when you account for all 
# the potential [states of the world] rises in sea-level by the year 2100?
nw.storm.slr.probs <- apply(slr.storm.DR2100, 2, inv.sf, val= fit_q_year[year100prob]/100)

# Check if probabilities are estimated for each state of the world
# If so then the minimum SLR plus the maximum storm surge should be equal to or larger than
# the 100-yr storm surge value
print(fit_q_year[year100prob]/100)
print(fit_q_year[storm_surgeL]/100 + min(proj2100DamsRes))

which(nw.storm.slr.probs == 'is.na') # integer(0) when no NAs

# If False then NAs will be returned for values with probabilities smaller than 1: 100,000
# Set all NAs to 0 since they are smaller than 1:100,000. This gives us a conservative estimate.
write.csv(nw.storm.slr.probs,"Data/NewStormProbabilities.csv") #,na="0")
storm.probs.table <- read.csv("Data/NewStormProbabilities.csv")
new.storm.probs = storm.probs.table[,2]

#Find the new probability for the Heberger 09 value:
new.ssurge.prob = mean(new.storm.probs)

remove(nw.storm.slr.probs)
#------------------- Estimate the 100-yr storm surge accounting for all potential SLR anomalies 
# Create a range of values from:
# The smallest storm surge value plus the smallest SLR anomaly
small.num = fit_q_year[2]/100 + min(proj2100DamsRes)
small.num = round(small.num, 1)
end.num = round(mean.stormSLR.2100HET[year100prob], 1)
num.range1 = seq(small.num, end.num, length.out=20)
# To 1 and 1/2 the 100-yr storm surge plus mean sea-level anomaly
max.pot.num = mean.stormSLR.2100HET[year100prob] + (mean.stormSLR.2100HET[year100prob])/2
max.pot.num = round(max.pot.num, 1)
num.range2 = seq(end.num+0.1, max.pot.num, length.out=80)
num.range = c(num.range1, num.range2)

# num.range = seq(small.num, max.pot.num, length.out=200)

# Find the probabilities of the values in the range using all the potential storm surge plus SLR anomalies
new.probs <- mat.or.vec(length(slr.storm.DR2100[1,]), length(num.range))
for(i in 1:length(num.range)){
  new.probs[,i] <- apply(slr.storm.DR2100, 2, inv.sf, val= num.range[i])
}

# Set all NAs to 0 since they are smaller than 1:100,000. This gives us a conservative estimate.
write.csv(new.probs,"Data/NewUncertaintyProbabilities.csv",na="0")
uncertainty.probs.table <- read.csv("Data/NewUncertaintyProbabilities.csv")
new.uncertainty.probs <- uncertainty.probs.table[,2:101]

average.uncertainty.probs <- rep(NA, length(num.range))
for(i in 1:length(num.range)){
  average.uncertainty.probs[i] <- mean(new.uncertainty.probs[,i])
}

#Find the new 100-yr probability
new.returnLevel <- which(average.uncertainty.probs >= 0.01)
max.returnL <- which.max(new.returnLevel)

print(average.uncertainty.probs[max.returnL]) #The probability should be 0.01
#If not try max.returnL +1
print(average.uncertainty.probs[max.returnL+1])

remove(new.probs)
##################### PRINT THE NEW PROBABILITY OF OCCURRENCE #####################
# Current 100-yr storm surge values:
round(fit_q_year[year100prob]/100,2) #Current storm surge for San Fran.
round(mean.stormSLR.2100HET[year100prob],2) #Current storm surge + Mean SLR for San Fran.
round(heberger09.2100[year100prob],2) #Current storm surge + SLR for San Fran.
round(num.range[max.returnL],2)  # New storm surge + SLR for San Fran. [Accounting for uncertainty]

#Factor of increase:
# From Heberger et al 2009 to new estimate:
# 0.01 to:
round(new.average.prob, 2)

# From mean slr estimate to new estimate:
# 0.01 to:
round(new.meanMSLR.prob, 2)

# From current storm surge (w/o slr) estimate to new estimate:
# 0.01 to:
round(new.ssurge.prob, 2)

# Save the workspace for plotting
save.image(file = "Workspace/SanFranSLR_StormSurge_Analysis.RData")

################################ PLOTTING ##########################################
# source("PlotRuckertetal_SanFranStormSurSLR.R")
################################## END #############################################

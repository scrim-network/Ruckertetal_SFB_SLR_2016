#########################################################################
# file: SLR_StormSurge_100yrFlood.R
#------------------------------------------------------------------------
# Author and copyright: Kelsey Ruckert
# Pennsylvania State University
# klr324@psu.edu
# Code written Dec. 2015, updated Nov. 2016
#
##==============================================================================
## Copyright 2015 Kelsey Ruckert
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
#
# NOTE: This script takes 8+ hours to run
#------------------------------------------------------------------------
# Estimates potential future 100-yr flood for the San Francisco Bay area
# based on the current 100-yr storm surge (San_francisco_tides.R) plus
# various SLR projections in the year 2100 (Project_global_sealevel.R)
#
# The various SLR projections include our mean estimate, our 90% CI estimate,
# the estimate in Heberger et al. 2009, and accounting for the SLR uncertainty
# in 2100.
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
load("Workspace/SFB_globalMeanSeaLevel2.RData")

############################# STORM SURGE ESTIMATION ##################################
#------------------- Estimate current storm surge -------------------------------------
storm_surgeL <- 1e5 # desired storm surge level
q = seq(0,1,length.out= storm_surgeL +1)  # quantile array

# Find closed-form solution of GEV fit, Equation (5).
fit_q_year = qgev(q, year.res.max.fit2@fit$par.ests[1], year.res.max.fit2@fit$par.ests[2],
                  year.res.max.fit2@fit$par.ests[3])
fit_q_year = fit_q_year[fit_q_year< max(fit_q_year)]

# Find which q is the 100-yr value 
num = which(q <= 0.99)
num.max = which.max(num)
year100prob <- num.max +1
# check to make sure the year100prob value is the 100-yr value (10^-2 or 0.01)
print(1-q[year100prob])

# We know the first value in the GEV fit [107cm] has a 100% probability when we account for any
# amount of future sea-level rise. For finding the inverse of the survival function any values
# that due not exist will be returned as an NA. To return a value we will add 107cm will a 100%
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
# If the first the value is -inf, set the first value to 107 (the second value)
# which has a probability of 100%.
print(1-q[1]); print(fit_q_year[1])

#Estimate the storm surge + SLR in heberger et al. 2009
# Add 107 cm with a probability of 100% to the vector
heberger09.2100 <- rep(NA, length(fit_q_year))
heberger09.2100[1] = fit_q_year[2]/100
heberger09.2100[2:length(heberger09.2100)] = (fit_q_year[2:length(fit_q_year)]/100) + a2new[45]

# Uncomment to account for water behind and in dams and reserviors (not done in main analysis)
#heberger09.2100[2:length(heberger09.2100)] = (fit_q_year[2:length(fit_q_year)]/100) + 1.38

#------------------- Storm surge plus mean sea-level rise ---------------------
#Estimate the mean storm surge + SLR based on the mean of the SLR uncertainty from MCMC
# heteroskedastic calibration
# Add 107 cm with a probability of 100% to the vector
mean.stormSLR.2100HET <- rep(NA, length(fit_q_year))
mean.stormSLR.2100HET[1] = fit_q_year[2]/100
mean.stormSLR.2100HET[2:length(mean.stormSLR.2100HET)] = (fit_q_year[2:length(fit_q_year)]/100) + mean(prob_proj2100/100)

# Uncomment to account for water behind and in dams and reserviors (not done in main analysis)
#mean.stormSLR.2100HET[2:length(mean.stormSLR.2100HET)] = (fit_q_year[2:length(fit_q_year)]/100) + mean(proj2100DamsRes)

#------------------- Storm surge plus 90% quantile sea-level rise ---------------------
#Estimate the storm surge + 90% SLR based on the mean of the SLR uncertainty from MCMC
# heteroskedastic calibration

# Add 107 cm with a probability of 100% to the vector
quant90.stormSLR.2100HET <- rep(NA, length(fit_q_year))
quant90.stormSLR.2100HET[1] = fit_q_year[2]/100
quant90.stormSLR.2100HET[2:length(quant90.stormSLR.2100HET)] = (fit_q_year[2:length(fit_q_year)]/100) + quantile(prob_proj2100/100, 0.90)

# Uncomment to account for water behind and in dams and reserviors (not done in main analysis)
#quant90.stormSLR.2100HET[2:length(quant90.stormSLR.2100HET)] = (fit_q_year[2:length(fit_q_year)]/100) + quantile(proj2100DamsRes, 0.90)

#------------------- Storm surge plus potential sea-level rises ---------------
#Estimate all the potential storm surge + SLRs based on MCMC heteroskedastic calibration
slr.storm.DR2100 <- mat.or.vec(length(fit_q_year),length(prob_proj2100))
for(i in 1:length(prob_proj2100)){
    slr.storm.DR2100[1,i] <- fit_q_year[2]/100
    slr.storm.DR2100[2:length(fit_q_year), i] <- (fit_q_year[2:length(fit_q_year)]/100) + (prob_proj2100[i]/100)
    
    # Uncomment to account for water behind and in dams and reserviors (not done in main analysis)
    #slr.storm.DR2100[2:length(fit_q_year), i] <- (fit_q_year[2:length(fit_q_year)]/100) + proj2100DamsRes[i]
}

##################### ESTIMATE NEW PROBABILITY OF OCCURRENCE #####################
#------------------- Estimating new Heberger et al. 2009 probabilities --------
# Storm surge plus Heberger et al. 2009
# What is the probability of the 100-yr storm surge associated with Heberger et al. 2009
# when you account for all the potential [states of the world] rises in sea-level by the year 2100?
heb.probs <- apply(slr.storm.DR2100, 2, inv.sf, val=heberger09.2100[year100prob])

# Check if probabilities are estimated for each state of the world
# If so, then the minimum SLR plus the maximum storm surge should be equal to or larger than
# the 100-yr Heberger 09 value
print(heberger09.2100[year100prob])
print(fit_q_year[100000]/100 + min(prob_proj2100/100))

# If not, all NAs produced have a probability smaller than 1:100,000.
# Set all NAs to 0 since they are smaller than 1:100,000. This gives us a conservative estimate,
# but the difference/impact of doing so is minor. This effect was tested in the plotting script.
# See for example supplementary figure 3.
write.csv(heb.probs,"Workspace/NewHebergerProbabilities_norm.csv",na="0")
heb.probs.table <- read.csv("Workspace/NewHebergerProbabilities_norm.csv")
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
# If so, then the minimum SLR plus the maximum storm surge should be equal to or larger than
# the 100-yr mean SLR + storm surge value
print(mean.stormSLR.2100HET[year100prob])
print(fit_q_year[storm_surgeL]/100 + min(prob_proj2100/100))

# If not, all NAs produced have a probability smaller than 1:100,000.
# Set all NAs to 0 since they are smaller than 1:100,000. This gives us a conservative estimate,
# but the difference/impact of doing so is minor. This effect was tested in the plotting script.
# See for example supplementary figure 3.
write.csv(mean.slr.probs,"Workspace/NewMeanProbabilities_norm.csv",na="0")
mean.probs.table <- read.csv("Workspace/NewMeanProbabilities_norm.csv")
new.mean.probs = mean.probs.table[,2]

#Find the new probability for the mean value:
new.meanMSLR.prob = mean(new.mean.probs)

remove(mean.slr.probs)
#------------------- Estimating new 90th quantile probabilities -----------------------
# Storm surge plus 90% SLR estimate
# What is the probability of the 100-yr storm surge associated with 90% SLR
# when you account for all the potential [states of the world] rises in sea-level by the year 2100?
quant90.slr.probs <- apply(slr.storm.DR2100, 2, inv.sf, val=quant90.stormSLR.2100HET[year100prob])

# Check if probabilities are estimated for each state of the world
# If so, then the minimum SLR plus the maximum storm surge should be equal to or larger than
# the 100-yr 90%  SLR + storm surge value
print(quant90.stormSLR.2100HET[year100prob])
print(fit_q_year[storm_surgeL]/100 + min(prob_proj2100/100))

# If not, all NAs produced have a probability smaller than 1:100,000.
# Set all NAs to 0 since they are smaller than 1:100,000. This gives us a conservative estimate,
# but the difference/impact of doing so is minor. This effect was tested in the plotting script.
# See for example supplementary figure 2.
write.csv(quant90.slr.probs,"Workspace/Newquant90Probabilities_norm.csv",na="0")
quant90.probs.table <- read.csv("Workspace/Newquant90Probabilities_norm.csv")
new.quant90.probs = quant90.probs.table[,2]

#Find the new probability for the mean value:
new.quant90MSLR.prob = mean(new.quant90.probs)

remove(quant90.slr.probs)
#------------------- Estimating new current storm surge probabilities ----------
# Storm surge
# What is the probability of the 100-yr storm surge associated when you account for all 
# the potential [states of the world] rises in sea-level by the year 2100?
nw.storm.slr.probs <- apply(slr.storm.DR2100, 2, inv.sf, val= fit_q_year[year100prob]/100)

# Check if probabilities are estimated for each state of the world
# If so, then the minimum SLR plus the maximum storm surge should be equal to or larger than
# the 100-yr storm surge value
print(fit_q_year[year100prob]/100)
print(fit_q_year[storm_surgeL]/100 + min(prob_proj2100/100))

which(nw.storm.slr.probs == 'is.na') # integer(0) when no NAs

# If not, all NAs produced have a probability smaller than 1:100,000.
# Set all NAs to 0 since they are smaller than 1:100,000. This gives us a conservative estimate.
write.csv(nw.storm.slr.probs,"Workspace/NewStormProbabilities_norm.csv") #,na="0")
storm.probs.table <- read.csv("Workspace/NewStormProbabilities_norm.csv")
new.storm.probs = storm.probs.table[,2]

#Find the new probability for the Heberger 09 value:
new.ssurge.prob = mean(new.storm.probs)

remove(nw.storm.slr.probs)
#------------------- Estimate the 100-yr storm surge accounting for all potential SLR anomalies 
# Create a range of values from:
# The smallest storm surge value plus the smallest SLR anomaly
small.num = fit_q_year[2]/100 + min(prob_proj2100/100)
small.num = round(small.num, 1)
end.num = round(mean.stormSLR.2100HET[year100prob], 1)
num.range1 = seq(small.num, end.num, length.out=20)
# To 1 and 1/2 the 100-yr storm surge plus mean sea-level anomaly
max.pot.num = mean.stormSLR.2100HET[year100prob] + (mean.stormSLR.2100HET[year100prob])/2
max.pot.num = round(max.pot.num, 1)
num.range2 = seq(end.num+0.1, max.pot.num, length.out=80)
num.range = c(num.range1, num.range2)

# num.range = seq(small.num, max.pot.num, length.out=200)

# Test: The first number in num.range must be larger than the minimum value in each flood frequency curve
# That way all NAs produced have a frequency smaller than 1:100,000.
print(num.range[1])
print(all(slr.storm.DR2100[1,] == min(slr.storm.DR2100)))
print(min(slr.storm.DR2100))

# Find the probabilities of the values in the range using all the potential storm surge plus SLR anomalies
new.probs <- mat.or.vec(length(slr.storm.DR2100[1,]), length(num.range))
for(i in 1:length(num.range)){
  new.probs[,i] <- apply(slr.storm.DR2100, 2, inv.sf, val= num.range[i])
}

# Set all NAs to 0 since they are smaller than 1:100,000. This gives us a conservative estimate.
write.csv(new.probs,"Workspace/NewUncertaintyProbabilities_norm.csv",na="0")
uncertainty.probs.table <- read.csv("Workspace/NewUncertaintyProbabilities_norm.csv")
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
round(fit_q_year[year100prob]/100,2) # Current storm surge for San Fran.
round(mean.stormSLR.2100HET[year100prob],2) # Current storm surge + Mean SLR for San Fran.
round(quant90.stormSLR.2100HET[year100prob],2) # Current storm surge + 90% SLR for San Fran.
round(heberger09.2100[year100prob],2) # Current storm surge + SLR for San Fran.
round(num.range[max.returnL+1],2)  # New storm surge + SLR for San Fran. [Accounting for uncertainty]

#Factor of increase:
# From Heberger et al 2009 to new estimate:
# 0.01 to:
round(new.average.prob, 2)

# From mean slr estimate to new estimate:
# 0.01 to:
round(new.meanMSLR.prob, 2)

# From 90% slr estimate to new estimate:
# 0.01 to:
round(new.quant90MSLR.prob, 2)

# From current storm surge (w/o slr) estimate to new estimate:
# 0.01 to:
round(new.ssurge.prob, 2)

# Save the workspace for plotting
save.image(file = "Workspace/SanFranSLR_StormSurge_Analysis_norm.RData")
################################## END #############################################

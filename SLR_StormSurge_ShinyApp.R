#########################################################################
# file: SLR_StormSurge_ShinyApp.R
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
# NOTE: This script takes 5+ hours to run
#------------------------------------------------------------------------
# Estimates potential future 100-yr flood for the San Francisco Bay area
# based on the current 100-yr storm surge (San_francisco_tides.R) plus
# various SLR projections in multiple years (Project_global_sealevel.R)
#
# The various SLR projections include our mean estimate and accounting for
# the SLR uncertainty in the year 2020, 2050, 2080, and 2100.
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
source("../Scripts/inv_sf.r")

#------------------- Load in the San Francisco tide & sea-level rise workspaces 
# San Francisco tides
load("Workspace/SFB_stormSurge.RData")

# Global mean sea-level
load("Workspace/SFB_globalMeanSeaLevel2.RData")

############################# STORM SURGE ESTIMATION ##################################
#------------------- Estimate current storm surge -------------------------------------
storm_surgeL <- 1e4 # desired storm surge level
q = seq(0,1,length.out= storm_surgeL +1)  # quantile array

# Find closed-form solution of GEV fit, Equation (5).
fit_q_year = qgev(q, year.res.max.fit2@fit$par.ests[1], year.res.max.fit2@fit$par.ests[2],
                  year.res.max.fit2@fit$par.ests[3])
fit_q_year = fit_q_year[fit_q_year< max(fit_q_year)]

# Find which q is the 100-yr value 
num = which(q <= 0.99)
num.max = which.max(num)
year100prob <- num.max +1
# check to make sure the year100prob value is the 100-yr value (10^-2)
print(1-q[year100prob])

freq_mat = matrix(c(1-q[1:(storm_surgeL)], fit_q_year/100), nrow = storm_surgeL, ncol=2)

############################ SLR PROJECTION #############################################
#--------------------- Step 3: Recreate the Rahmstorf predictions: ----------------
# Set up the parameters specified in Rahmstorf 2007:
# The parameters in the model include (in order of appreciance in the vector)
# a = sensitivity of sea-level to changes in temperature. Given as meters/Celsius/year
# T0 = base temperature or temperature when sea-level is zero. Given as Celius
# Initial = the value of sea-level in the year 1880. 1880 is the starting year for this
#       analysis. Given as meters.

# original = c(0.34, -0.5, slr[1]) #slr[1] is equal to -0.1464 meters
# 
# # Run the projections using the emission scenarios used in Rahmstorf 2007
# from = 2 # Start from the 2nd year since the first value is defined as a uncertain parameter
# to = 45     # There are 45 points in this projection from 1880 to 2100
# timestep = 5 # Projections are in intervals of 5 years
# 
# # The names of the projections correspond to the names of the emission scenario
# a2_p = rahmfunction(original, a2)
# 
# # Make with respect to the year 2000:
# a2new = (a2_p$sle/100) - (a2_p$sle[25]/100)
# 
# #----------- Step 4: Calculate additional runoff from dams & reservoirs: ----------------
# # Heberger 2009: 1.4m (1.38m) from the A2 estimate using rahmstorf 2007 model plus
# # additional accounting for changing runoff from dams & reservoirs
# print(scenario_time)
# 
# DamsReservoirs2100 = 1.38 - a2new[45]
# DamsReservoirs2100 <- round(DamsReservoirs2100,2)
# 
# DamsReservoirs2080 = 0.95 - a2new[41]
# DamsReservoirs2080 <- round(DamsReservoirs2080,2)
# 
# DamsReservoirs2050 = 0.47 - a2new[35]
# DamsReservoirs2050 <- round(DamsReservoirs2050,2)
# 
# DamsReservoirs2020 = 0.12 - a2new[29]
# DamsReservoirs2020 <- round(DamsReservoirs2020,2)

# Find sea-level estimates in 2100
# NOTE: Uncomment the above an below code if you wish to account for water behind and in dams and 
# reserviors. We do not account for this in the paper or in the Shiny App.
projSLRDamsRes <- mat.or.vec(length(slr.mcmc_proj[,1]),4)
projSLRDamsRes[,1] <- (slr.mcmc_proj[,141]/100) #+ DamsReservoirs2020 #The year 2020 is the 221st number
projSLRDamsRes[,2] <- (slr.mcmc_proj[,171]/100) #+ DamsReservoirs2050 #The year 2050 is the 221st number
projSLRDamsRes[,3] <- (slr.mcmc_proj[,201]/100) #+ DamsReservoirs2080 #The year 2080 is the 221st number
projSLRDamsRes[,4] <- (slr.mcmc_proj[,221]/100) #+ DamsReservoirs2100 #The year 2100 is the 221st number

colnames(projSLRDamsRes, do.NULL = FALSE)
colnames(projSLRDamsRes) = c("2020", "2050", "2080", "2100")

write.csv(projSLRDamsRes,"ShinyApp/Data/Sea_level_proj_anomalies.csv")
write.csv(freq_mat,"ShinyApp/Data/Storm_surge_freq.csv")

############################# STORM SURGE ESTIMATION ##################################
#------------------- Storm surge plus potential sea-level rises ---------------
#Estimate all the potential storm surge + SLRs based on MCMC heteroskedastic calibration
#slr.storm.DR2100 <- mat.or.vec(storm_surgeL,length(proj2100DamsRes))
#for(i in 1:length(proj2100DamsRes)){
#  slr.storm.DR2100[,i] <- (fit_q_year/100) + proj2100DamsRes[i]
#}
# If the first the value is -inf, set the first value to 107 (the second value)
# which has a probability of 100%.
print(1-q[1]); print(fit_q_year[1])

#2020
slr.storm.2020 <- mat.or.vec(length(fit_q_year),length(projSLRDamsRes[,1]))
for(i in 1:length(projSLRDamsRes[,1])){
  slr.storm.2020[1,i] <- fit_q_year[2]/100
  slr.storm.2020[2:length(fit_q_year), i] <- (fit_q_year[2:length(fit_q_year)]/100) + projSLRDamsRes[i,1]
}
#2050
slr.storm.2050 <- mat.or.vec(length(fit_q_year),length(projSLRDamsRes[,2]))
for(i in 1:length(projSLRDamsRes[,2])){
  slr.storm.2050[1,i] <- fit_q_year[2]/100
  slr.storm.2050[2:length(fit_q_year), i] <- (fit_q_year[2:length(fit_q_year)]/100) + projSLRDamsRes[i,2]
}
#2080
slr.storm.2080 <- mat.or.vec(length(fit_q_year),length(projSLRDamsRes[,3]))
for(i in 1:length(projSLRDamsRes[,3])){
  slr.storm.2080[1,i] <- fit_q_year[2]/100
  slr.storm.2080[2:length(fit_q_year), i] <- (fit_q_year[2:length(fit_q_year)]/100) + projSLRDamsRes[i,3]
}
#2100
slr.storm.2100 <- mat.or.vec(length(fit_q_year),length(projSLRDamsRes[,4]))
for(i in 1:length(projSLRDamsRes[,4])){
  slr.storm.2100[1,i] <- fit_q_year[2]/100
  slr.storm.2100[2:length(fit_q_year), i] <- (fit_q_year[2:length(fit_q_year)]/100) + projSLRDamsRes[i,4]
}

# Year 2020
#------------------- Estimate the 100-yr storm surge accounting for all potential SLR anomalies 
# Create a range of values from:
# The smallest storm surge value plus the smallest SLR anomaly
small.num2020 = fit_q_year[2]/100
end.num2020 = round(fit_q_year[year100prob]/100, 1)
num.range12020 = seq(small.num2020, end.num2020, length.out=20)
# To 1 and 1/2 the 100-yr storm surge plus mean sea-level anomaly
max.pot.num2020 = fit_q_year[year100prob]/100 + (max(projSLRDamsRes[,1]))
max.pot.num2020 = round(max.pot.num2020, 1)
num.range22020 = seq(end.num2020+0.02, max.pot.num2020, length.out=80)
num.range2020 = c(num.range12020, num.range22020)

# Find the probabilities of the values in the range using all the potential storm surge plus SLR anomalies
new.probs2020 <- mat.or.vec(length(slr.storm.2020[1,]), length(num.range2020))
for(i in 1:length(num.range2020)){
  new.probs2020[,i] <- apply(slr.storm.2020, 2, inv.sf, val= num.range2020[i])
}

# Set all NAs to 0 since they are smaller than 1:100,000. This gives us a conservative estimate.
write.csv(new.probs2020,"ShinyApp/Data/NewUncertaintyProbs2020.csv",na="0")
uncertainty.probs.table2020 <- read.csv("ShinyApp/Data/NewUncertaintyProbs2020.csv")
new.uncertainty.probs2020 <- uncertainty.probs.table2020[,2:101]

average.uncertainty.probs2020 <- rep(NA, length(num.range2020))
for(i in 1:length(num.range2020)){
  average.uncertainty.probs2020[i] <- mean(new.uncertainty.probs2020[,i])
}

#Find the new 100-yr probability
new.returnLevel2020 <- which(average.uncertainty.probs2020 >= 0.01)
max.returnL2020 <- which.max(new.returnLevel2020)

print(average.uncertainty.probs2020[max.returnL2020]) #The probability should be 0.01
#If not try max.returnL +1
print(average.uncertainty.probs2020[max.returnL2020+1])

remove(new.probs2020)


accountUNC_2020 <- mat.or.vec(length(average.uncertainty.probs2020),3)
accountUNC_2020[,1] <- average.uncertainty.probs2020   # + DamsReservoirs #The year 2020 is the 221st number
accountUNC_2020[,2] <- num.range2020 # + DamsReservoirs #The year 2050 is the 221st number
accountUNC_2020[1,3] <- max.returnL2020 # + DamsReservoirs #The year 2080 is the 221st number

colnames(accountUNC_2020, do.NULL = FALSE)
colnames(accountUNC_2020) = c("ReturnPeriod", "ReturnLevel", "100LevelNum")

write.csv(accountUNC_2020,"ShinyApp/Data/StormFreqAccountUNC2020.csv")

#2050
#------------------- Estimate the 100-yr storm surge accounting for all potential SLR anomalies 
# Create a range of values from:
# The smallest storm surge value plus the smallest SLR anomaly
meanyear2050 = fit_q_year[year100prob]/100 + mean(projSLRDamsRes[,2])
small.num2050 = fit_q_year[2]/100
end.num2050 = round(meanyear2050/1.1, 1)
num.range12050 = seq(small.num2050, end.num2050, length.out=20)
# To 1 and 1/2 the 100-yr storm surge plus mean sea-level anomaly
max.pot.num2050 = meanyear2050 + mean(projSLRDamsRes[,2]) + 0.15
max.pot.num2050 = round(max.pot.num2050, 1)
num.range22050 = seq(end.num2050+0.05, max.pot.num2050, length.out=80)
num.range2050 = c(num.range12050, num.range22050)

# Find the probabilities of the values in the range using all the potential storm surge plus SLR anomalies
new.probs2050 <- mat.or.vec(length(slr.storm.2050[1,]), length(num.range2050))
for(i in 1:length(num.range2050)){
  new.probs2050[,i] <- apply(slr.storm.2050, 2, inv.sf, val= num.range2050[i])
}

# Set all NAs to 0 since they are smaller than 1:100,000. This gives us a conservative estimate.
write.csv(new.probs2050,"ShinyApp/Data/NewUncertaintyProbs2050.csv",na="0")
uncertainty.probs.table2050 <- read.csv("ShinyApp/Data/NewUncertaintyProbs2050.csv")
new.uncertainty.probs2050 <- uncertainty.probs.table2050[,2:101]

average.uncertainty.probs2050 <- rep(NA, length(num.range2050))
for(i in 1:length(num.range2050)){
  average.uncertainty.probs2050[i] <- mean(new.uncertainty.probs2050[,i])
}

#Find the new 100-yr probability
new.returnLevel2050 <- which(average.uncertainty.probs2050 >= 0.010)
max.returnL2050 <- which.max(new.returnLevel2050)

print(average.uncertainty.probs2050[max.returnL2050]) #The probability should be 0.01
#If not try max.returnL +1
print(average.uncertainty.probs2050[max.returnL2050+1])

accountUNC_2050 <- mat.or.vec(length(average.uncertainty.probs2050),3)
accountUNC_2050[,1] <- average.uncertainty.probs2050   # + DamsReservoirs #The year 2020 is the 221st number
accountUNC_2050[,2] <- num.range2050 # + DamsReservoirs #The year 2050 is the 221st number
accountUNC_2050[1,3] <- max.returnL2050+1 # + DamsReservoirs #The year 2080 is the 221st number

colnames(accountUNC_2050, do.NULL = FALSE)
colnames(accountUNC_2050) = c("ReturnPeriod", "ReturnLevel", "100LevelNum")

write.csv(accountUNC_2050,"ShinyApp/Data/StormFreqAccountUNC2050.csv")

remove(new.probs2050)

# 2080
#------------------- Estimate the 100-yr storm surge accounting for all potential SLR anomalies 
# Create a range of values from:
# The smallest storm surge value plus the smallest SLR anomaly
meanyear2080 = fit_q_year[year100prob]/100 + mean(projSLRDamsRes[,3])
small.num2080 = fit_q_year[2]/100
end.num2080 =  round(meanyear2080, 1)
num.range12080 = seq(small.num2080, end.num2080, length.out=20)
# To 1 and 1/2 the 100-yr storm surge plus mean sea-level anomaly
max.pot.num2080 =  meanyear2080 + mean(projSLRDamsRes[,3]) + 0.1
max.pot.num2080 = round(max.pot.num2080, 1)
num.range22080 = seq(end.num2080+0.05, max.pot.num2080, length.out=80)
num.range2080 = c(num.range12080, num.range22080)

# Find the probabilities of the values in the range using all the potential storm surge plus SLR anomalies
new.probs2080 <- mat.or.vec(length(slr.storm.2080[1,]), length(num.range2080))
for(i in 1:length(num.range2080)){
  new.probs2080[,i] <- apply(slr.storm.2080, 2, inv.sf, val= num.range2080[i])
}

# Set all NAs to 0 since they are smaller than 1:100,000. This gives us a conservative estimate.
write.csv(new.probs2080,"ShinyApp/Data/NewUncertaintyProbs2080.csv",na="0")
uncertainty.probs.table2080 <- read.csv("ShinyApp/Data/NewUncertaintyProbs2080.csv")
new.uncertainty.probs2080 <- uncertainty.probs.table2080[,2:101]

average.uncertainty.probs2080 <- rep(NA, length(num.range2080))
for(i in 1:length(num.range2080)){
  average.uncertainty.probs2080[i] <- mean(new.uncertainty.probs2080[,i])
}

#Find the new 100-yr probability
new.returnLevel2080 <- which(average.uncertainty.probs2080 >= 0.01)
max.returnL2080 <- which.max(new.returnLevel2080)

print(average.uncertainty.probs2080[max.returnL2080]) #The probability should be 0.01
#If not try max.returnL +1
print(average.uncertainty.probs2080[max.returnL2080+1])

accountUNC_2080 <- mat.or.vec(length(average.uncertainty.probs2080),3)
accountUNC_2080[,1] <- average.uncertainty.probs2080   # + DamsReservoirs #The year 2020 is the 221st number
accountUNC_2080[,2] <- num.range2080 # + DamsReservoirs #The year 2050 is the 221st number
accountUNC_2080[1,3] <- max.returnL2080 # + DamsReservoirs #The year 2080 is the 221st number

colnames(accountUNC_2080, do.NULL = FALSE)
colnames(accountUNC_2080) = c("ReturnPeriod", "ReturnLevel", "100LevelNum")

write.csv(accountUNC_2080,"ShinyApp/Data/StormFreqAccountUNC2080.csv")

# #------------------- Estimate the 100-yr storm surge accounting for all potential SLR anomalies 
# # Create a range of values from:
# # The smallest storm surge value plus the smallest SLR anomaly
# small.num = fit_q_year[2]/100 + min(proj2100DamsRes)
# small.num = round(small.num, 1)
# end.num = round(mean.stormSLR.2100HET[year100prob], 1)
# num.range1 = seq(small.num, end.num, length.out=20)
# # To 1 and 1/2 the 100-yr storm surge plus mean sea-level anomaly
# max.pot.num = mean.stormSLR.2100HET[year100prob] + (mean.stormSLR.2100HET[year100prob])/2
# max.pot.num = round(max.pot.num, 1)
# num.range2 = seq(end.num+0.1, max.pot.num, length.out=80)
# num.range = c(num.range1, num.range2)
# 
# # Find the probabilities of the values in the range using all the potential storm surge plus SLR anomalies
# new.probs <- mat.or.vec(length(slr.storm.2100[1,]), length(num.range))
# for(i in 1:length(num.range)){
#   new.probs[,i] <- apply(slr.storm.2100, 2, inv.sf, val= num.range[i])
# }
# 
# # Set all NAs to 0 since they are smaller than 1:100,000. This gives us a conservative estimate.
# write.csv(new.probs,”Workspace/NewUncertaintyProbabilities.csv",na="0")
# uncertainty.probs.table <- read.csv(“Workspace/NewUncertaintyProbabilities.csv")
# new.uncertainty.probs <- uncertainty.probs.table[,2:101]
# 
# average.uncertainty.probs <- rep(NA, length(num.range))
# for(i in 1:length(num.range)){
#   average.uncertainty.probs[i] <- mean(new.uncertainty.probs[,i])
# }

#Find the new 100-yr probability
# new.returnLevel <- which(average.uncertainty.probs >= 0.01)
# max.returnL <- which.max(new.returnLevel)

# Load in the San Francisco tide & sea-level rise workspaces 
load("Workspace/SanFranSLR_StormSurge_Analysis_norm.RData")

print(average.uncertainty.probs[max.returnL]) #The probability should be 0.01
#If not try max.returnL +1
print(average.uncertainty.probs[max.returnL+1])

accountUNC_2100 <- mat.or.vec(length(num.range),3)
accountUNC_2100[,1] <- average.uncertainty.probs   # + DamsReservoirs #The year 2020 is the 221st number
accountUNC_2100[,2] <- num.range # + DamsReservoirs #The year 2050 is the 221st number
accountUNC_2100[1,3] <- max.returnL # + DamsReservoirs #The year 2080 is the 221st number

colnames(accountUNC_2100, do.NULL = FALSE)
colnames(accountUNC_2100) = c("ReturnPeriod", "ReturnLevel", "100LevelNum")

write.csv(accountUNC_2100,”ShinyApp/Data/StormFreqAccountUNC2100.csv")
##################### PRINT THE NEW PROBABILITY OF OCCURRENCE #####################
# Current 100-yr storm surge values:
# round(num.range[max.returnL+1],2)  # New storm surge + SLR for San Fran. [Accounting for uncertainty]

################################## END #############################################

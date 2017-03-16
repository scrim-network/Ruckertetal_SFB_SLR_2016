#########################################################################
# file: Distribution_test.R
#------------------------------------------------------------------------
# Author and copyright: Kelsey Ruckert
# Pennsylvania State University
# klr324@psu.edu
# Code written Nov. 2016
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
#
# NOTE: This script takes 5+ hours to run
#------------------------------------------------------------------------
# Estimates potential future 100-yr flood for the San Francisco Bay area
# based on the current 100-yr storm surge (San_francisco_tides.R) plus
# various SLR projections in the year 2100 (Project_global_sealevel.R)
# The various SLR projections include our distribution, a normal distribution,
# a log-normal distribution, and a pareto distribution.
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

# Global mean sea level
load("Workspace/SFB_globalMeanSeaLevel2.RData")

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

#------------------- Storm surge plus potential sea-level rises ---------------
#Estimate all the potential storm surge + SLRs based on MCMC heteroskedastic calibration

slr.storm.DR2100 <- mat.or.vec(length(fit_q_year),length(prob_proj2100))
for(i in 1:length(prob_proj2100)){
    slr.storm.DR2100[1,i] <- fit_q_year[2]/100
    slr.storm.DR2100[2:length(fit_q_year), i] <- (fit_q_year[2:length(fit_q_year)]/100) + (prob_proj2100[i]/100)
    #slr.storm.DR2100[2:length(fit_q_year), i] <- (fit_q_year[2:length(fit_q_year)]/100) + proj2100DamsRes[i]
}

#------------------- Normal distribution ---------------
# Approximate normal distribution from SLR mean and standard deviation
normal_dist_2100 = rnorm(1e4, mean = mean(prob_proj2100/100), sd = sd(prob_proj2100/100))

slr.storm.NormalD <- mat.or.vec(length(fit_q_year),length(normal_dist_2100))
for(i in 1:length(normal_dist_2100)){
    slr.storm.NormalD[1,i] <- fit_q_year[2]/100 + min(normal_dist_2100) # the minimum value is negative
    slr.storm.NormalD[2:length(fit_q_year), i] <- (fit_q_year[2:length(fit_q_year)]/100) + normal_dist_2100[i]
}

#------------------- Log Normal distribution ---------------
# Approximate log normal distribution from SLR mean and standard deviation on a log scale
meanlog <- log(mean(prob_proj2100/100)) - 0.5 * log(1 + (sd(prob_proj2100/100)^2)/(mean(prob_proj2100/100)^2))
sdlog <- sqrt(log(1 + (sd(prob_proj2100/100)^2)/(mean(prob_proj2100/100)^2)))

lognorm_dist_2100 = rlnorm(1e4, mean = meanlog, sd = sdlog)

slr.storm.log_norm <- mat.or.vec(length(fit_q_year),length(lognorm_dist_2100))
for(i in 1:length(lognorm_dist_2100)){
    slr.storm.log_norm[1,i] <- fit_q_year[2]/100
    slr.storm.log_norm[2:length(fit_q_year), i] <- (fit_q_year[2:length(fit_q_year)]/100) + lognorm_dist_2100[i]
}

#------------------- Pareto distribution ---------------
# pareto r script (written by Yawen Guan)

# b is the lower cutoff of the distribution. saying smaller values than b
# is not allowed. From optim, the cutoff is set to 0.4743678.
# a is the shape parameter.

library(DEoptim)
opt = function(x){
    alpha  =x[1]
    xm = x[2]
    (mean(prob_proj2100/100) - alpha*xm/(alpha-1))^2 +((sd(prob_proj2100/100)^2) - alpha*xm^2/((alpha-2)*(alpha-1)^2 ))^2
}
out = DEoptim(opt,lower = c(2.00001,0.00001),upper=c(100,10))$optim$bestmem
print(out)
# par1 is the parameter for pareto distribution
# par1      par2
# 3.7516450 0.4537994

xm = out[2]
alpha = out[1]
cat("mean",alpha*xm/(alpha-1))
cat("sd",sqrt(alpha*xm^2/((alpha-2)*(alpha-1)^2 )))

# generate pareto
# Create funciton as http://sites.stat.psu.edu/~dhunter/astrostatistics/PennState2010/simboot.html
qpareto <- function(u, a=alpha, b=xm) b/(1-u)^(1/a)
rpareto <- function(n, a=alpha, b=xm) qpareto(runif(n),a,b)

pareto_dist_2100 = rpareto(n=1e4, b=xm, a=alpha)

# double check the mean and sd
mean(pareto_dist_2100)
sd(pareto_dist_2100)

slr.storm.pareto <- mat.or.vec(length(fit_q_year),length(pareto_dist_2100))
for(i in 1:length(pareto_dist_2100)){
    slr.storm.pareto[1,i] <- fit_q_year[2]/100
    slr.storm.pareto[2:length(fit_q_year), i] <- (fit_q_year[2:length(fit_q_year)]/100) + pareto_dist_2100[i]
}

##################### ESTIMATE NEW PROBABILITY OF OCCURRENCE #####################
#------------------- Estimate the 100-yr storm surge accounting for all potential SLR anomalies
# Create a range of values from:
# The smallest storm surge value plus the smallest SLR anomaly
small.num = fit_q_year[2]/100 + min(prob_proj2100/100)
small.num = round(small.num, 1)
end.num = round(2.21, 1)
num.range1 = seq(small.num, end.num, length.out=20)
# To 1 and 1/2 the 100-yr storm surge plus mean sea-level anomaly
max.pot.num = 2.21 + (2.21)/2
max.pot.num = round(max.pot.num, 1)
num.range2 = seq(end.num+0.1, max.pot.num, length.out=80)
num.range = c(num.range1, num.range2)

# Test: The first number in num.range must be larger than the minimum value in each flood frequency curve
# That way all NAs produced have a frequency smaller than 1:100,000.
print(num.range[1])
print(all(slr.storm.DR2100[1,] == min(slr.storm.DR2100)))
print(min(slr.storm.DR2100))

# Set all NAs to 0 since they are smaller than 1:100,000. This gives us a conservative estimate.
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

#------------------- Estimate the 100-yr storm surge accounting for Normal approximation
# Create a range of values from:
# The smallest storm surge value plus the smallest SLR anomaly
norm.small.num = fit_q_year[2]/100 + min(normal_dist_2100)
norm.small.num = round(norm.small.num, 1)
norm.end.num = round(2.21, 1)
norm.num.range1 = seq(norm.small.num, norm.end.num, length.out=20)
# To 1 and 1/2 the 100-yr storm surge plus mean sea-level anomaly (2.21)
norm.max.pot.num = 2.21 + (2.21)/2
norm.max.pot.num = round(norm.max.pot.num, 1)
norm.num.range2 = seq(norm.end.num+0.1, norm.max.pot.num, length.out=80)
norm.num.range = c(norm.num.range1, norm.num.range2)

# Test: The first number in num.range must be larger than the minimum value in each flood frequency curve
# That way all NAs produced have a frequency smaller than 1:100,000.
print(norm.num.range[1])
print(all(slr.storm.NormalD[1,] == min(slr.storm.NormalD)))
print(min(slr.storm.NormalD))

# num.range = seq(small.num, max.pot.num, length.out=200)

# Find the probabilities of the values in the range using all the potential storm surge plus SLR anomalies
norm.new.probs <- mat.or.vec(length(slr.storm.NormalD[1,]), length(norm.num.range))
for(i in 1:length(norm.num.range)){
    norm.new.probs[,i] <- apply(slr.storm.NormalD, 2, inv.sf, val= norm.num.range[i])
}

# Set all NAs to 0 since they are smaller than 1:100,000. This gives us a conservative estimate.
write.csv(norm.new.probs,"Workspace/NewUncertaintyProbabilities_NormalDist.csv",na="0")
norm.uncertainty.probs.table <- read.csv("Workspace/NewUncertaintyProbabilities_NormalDist.csv")
norm.new.uncertainty.probs <- norm.uncertainty.probs.table[,2:101]

save.image(file = "Workspace/SanFranSLR_StormSurge_Analysis_distribution_test.RData")

norm.average.uncertainty.probs <- rep(NA, length(norm.num.range))
for(i in 1:length(norm.num.range)){
    norm.average.uncertainty.probs[i] <- mean(norm.new.uncertainty.probs[,i])
}

#Find the new 100-yr probability
norm.new.returnLevel <- which(norm.average.uncertainty.probs >= 0.01)
norm.max.returnL <- which.max(norm.new.returnLevel)

print(norm.average.uncertainty.probs[norm.max.returnL]) #The probability should be 0.01
#If not try max.returnL +1
print(norm.average.uncertainty.probs[norm.max.returnL+1])

remove(norm.new.probs)

#------------------- Estimate the 100-yr storm surge accounting for Log Normal approximation
# Create a range of values from:
# The smallest storm surge value plus the smallest SLR anomaly
Lnorm.small.num = fit_q_year[2]/100 + min(lognorm_dist_2100)
Lnorm.small.num = round(Lnorm.small.num, 1)
Lnorm.end.num = round(2.21, 1)
Lnorm.num.range1 = seq(Lnorm.small.num, Lnorm.end.num, length.out=20)
# To 1 and 1/2 the 100-yr storm surge plus mean sea-level anomaly (2.21)
Lnorm.max.pot.num = 2.21 + (2.21)/2
Lnorm.max.pot.num = round(Lnorm.max.pot.num, 1)
Lnorm.num.range2 = seq(Lnorm.end.num+0.1, Lnorm.max.pot.num, length.out=80)
Lnorm.num.range = c(Lnorm.num.range1, Lnorm.num.range2)

# Test: The first number in num.range must be larger than the minimum value in each flood frequency curve
# That way all NAs produced have a frequency smaller than 1:100,000.
print(Lnorm.num.range[1])
print(all(slr.storm.log_norm[1,] == min(slr.storm.log_norm)))
print(min(slr.storm.log_norm))

# num.range = seq(small.num, max.pot.num, length.out=200)

# Find the probabilities of the values in the range using all the potential storm surge plus SLR anomalies
Lnorm.new.probs <- mat.or.vec(length(slr.storm.log_norm[1,]), length(Lnorm.num.range))
for(i in 1:length(Lnorm.num.range)){
    Lnorm.new.probs[,i] <- apply(slr.storm.log_norm, 2, inv.sf, val= Lnorm.num.range[i])
}

# Set all NAs to 0 since they are smaller than 1:100,000. This gives us a conservative estimate.
write.csv(Lnorm.new.probs,"Workspace/NewUncertaintyProbabilities_LogNormalDist.csv",na="0")
Lnorm.uncertainty.probs.table <- read.csv("Workspace/NewUncertaintyProbabilities_LogNormalDist.csv")
Lnorm.new.uncertainty.probs <- Lnorm.uncertainty.probs.table[,2:101]

save.image(file = "Workspace/SanFranSLR_StormSurge_Analysis_distribution_test.RData")

Lnorm.average.uncertainty.probs <- rep(NA, length(Lnorm.num.range))
for(i in 1:length(Lnorm.num.range)){
    Lnorm.average.uncertainty.probs[i] <- mean(Lnorm.new.uncertainty.probs[,i])
}

#Find the new 100-yr probability
Lnorm.new.returnLevel <- which(Lnorm.average.uncertainty.probs >= 0.01)
Lnorm.max.returnL <- which.max(Lnorm.new.returnLevel)

print(Lnorm.average.uncertainty.probs[Lnorm.max.returnL]) #The probability should be 0.01
#If not try max.returnL +1
print(Lnorm.average.uncertainty.probs[Lnorm.max.returnL+1])

remove(Lnorm.new.probs)

#------------------- Estimate the 100-yr storm surge accounting for Pareto approximation
# Create a range of values from:
# The smallest storm surge value plus the smallest SLR anomaly
pareto.small.num = fit_q_year[2]/100 + min(pareto_dist_2100)
pareto.small.num = round(pareto.small.num, 1)
pareto.end.num = round(2.21, 1)
pareto.num.range1 = seq(pareto.small.num, pareto.end.num, length.out=20)
# To 1 and 1/2 the 100-yr storm surge plus mean sea-level anomaly (2.21)
pareto.max.pot.num = 2.21 + (2.21)/2
pareto.max.pot.num = round(pareto.max.pot.num, 1)
pareto.num.range2 = seq(pareto.end.num+0.1, pareto.max.pot.num, length.out=80)
pareto.num.range = c(pareto.num.range1, pareto.num.range2)

# Test: The first number in num.range must be larger than the minimum value in each flood frequency curve
# That way all NAs produced have a frequency smaller than 1:100,000.
print(pareto.num.range[1])
print(all(slr.storm.pareto[1,] == min(slr.storm.pareto)))
print(min(slr.storm.pareto))

# num.range = seq(small.num, max.pot.num, length.out=200)

# Find the probabilities of the values in the range using all the potential storm surge plus SLR anomalies
pareto.new.probs <- mat.or.vec(length(slr.storm.pareto[1,]), length(pareto.num.range))
for(i in 1:length(pareto.num.range)){
    pareto.new.probs[,i] <- apply(slr.storm.pareto, 2, inv.sf, val= pareto.num.range[i])
}

# Set all NAs to 0 since they are smaller than 1:100,000. This gives us a conservative estimate.
write.csv(pareto.new.probs,"Workspace/NewUncertaintyProbabilities_paretoDist.csv",na="0")
pareto.uncertainty.probs.table <- read.csv("Workspace/NewUncertaintyProbabilities_paretoDist.csv")
pareto.new.uncertainty.probs <- pareto.uncertainty.probs.table[,2:101]

save.image(file = "Workspace/SanFranSLR_StormSurge_Analysis_distribution_test.RData")

pareto.average.uncertainty.probs <- rep(NA, length(pareto.num.range))
for(i in 1:length(pareto.num.range)){
    pareto.average.uncertainty.probs[i] <- mean(pareto.new.uncertainty.probs[,i])
}

#Find the new 100-yr probability
pareto.new.returnLevel <- which(pareto.average.uncertainty.probs >= 0.01)
pareto.max.returnL <- which.max(pareto.new.returnLevel)

print(pareto.average.uncertainty.probs[pareto.max.returnL]) #The probability should be 0.01
#If not try max.returnL +1
print(pareto.average.uncertainty.probs[pareto.max.returnL+1])

remove(pareto.new.probs)

##################### PRINT THE NEW PROBABILITY OF OCCURRENCE #####################
# Current 100-yr storm surge values:
round(num.range[max.returnL+1],2)  # New storm surge + SLR for San Fran. [Accounting for uncertainty]
round(norm.num.range[norm.max.returnL+1],2)  # New storm surge + SLR for San Fran. [Normal approximation]
round(Lnorm.num.range[Lnorm.max.returnL],2)  # New storm surge + SLR for San Fran. [Log normal approximation]
round(Lnorm.num.range[pareto.max.returnL],2)  # New storm surge + SLR for San Fran. [Log normal approximation]

# Save the workspace for plotting
save.image(file = "Workspace/SanFranSLR_StormSurge_Analysis_distribution_test.RData")
################################## END #############################################

#########################################################################
# file: Input_preanalysed_data.R
#------------------------------------------------------------------------
# Author and copyright: Kelsey Ruckert
# Pennsylvania State University
# klr324@psu.edu
# Code written May 2016, updated Oct. 2016 & Jan. 2017
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
# Read in pre-analyzed data and set variables for the Shiny App.
#########################################################################

# Input block maxima data
load("Data/year_res_max.RData")
obs.x = coredata(year.res.max)/100
num.x <- length(obs.x)
sf <- seq(1,1/num.x,by=-1/num.x)
order.x <- order(obs.x)

# Input sea-level data
sea_data = read.csv("Data/Sea_level_proj_anomalies.csv")
med_x = mean(sea_data$X2050)
max_x = max(sea_data$X2100)

# Input storm surge data
storm_data = read.csv("Data/Storm_surge_freq.csv")
frequency = storm_data$V1
flood_meters = storm_data$V2
remove(storm_data)

MinMaxsf <- mat.or.vec(length(flood_meters), 2)
MinMaxsf[1, 1:2] = flood_meters[2]

new_line <- rep(NA, length(flood_meters))
new_line[1] = flood_meters[2]

slr.plotted <- rep(NA, length(flood_meters))
slr.plotted[1] = flood_meters[2]

# Find which q is the 100-yr value
my_num = which(frequency >= 0.01)
my_num.max = which.max(my_num)
year100prob <- my_num.max +1

# Input flood probability data which account for uncertainty in
# sea-level projections
#------------------------------------------------------------------------
stormFreqAcc2100 = read.csv("Data/StormFreqAccountUNC2100.csv")
returnperiods2100 <- stormFreqAcc2100[,2]
returnlevels2100 <- stormFreqAcc2100[,3]
return100_yr2100 <- stormFreqAcc2100[1,4] + 1

stormFreqAcc2080 = read.csv("Data/StormFreqAccountUNC2080.csv")
returnperiods2080 <- stormFreqAcc2080[,2]
returnperiods2080[1] <-returnperiods2100[1]
returnlevels2080 <- stormFreqAcc2080[,3]
return100_yr2080 <- stormFreqAcc2080[1,4] + 1

stormFreqAcc2050 = read.csv("Data/StormFreqAccountUNC2050.csv")
returnperiods2050 <- stormFreqAcc2050[,2]
returnperiods2050[1:4] <-returnperiods2100[1]
returnlevels2050 <- stormFreqAcc2050[,3]
return100_yr2050 <- stormFreqAcc2050[1,4] - 1

stormFreqAcc2020 = read.csv("Data/StormFreqAccountUNC2020.csv")
returnperiods2020 <- stormFreqAcc2020[,2]
returnlevels2020 <- stormFreqAcc2020[,3]
return100_yr2020 <- stormFreqAcc2020[1,4]

################################## END #############################################

#########################################################################
# file: Plot_hypospectic.R
#------------------------------------------------------------------------
# Author and copyright: Kelsey Ruckert
# Pennsylvania State University
# klr324@psu.edu
# Code written Aug. 2016
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
# Plot hypospectic figure for Ruckert et al. (in review)
#
#########################################################################
# Clear global environment
rm(list =ls())

# open packages
library(RColorBrewer)

# Create unique colors
myheatcolors <- brewer.pal(9,"YlOrRd")
mygreen <- brewer.pal(6,"Greens")

# Read in elevation data from subset of San Francisco county
sanfran_dem_sub = read.table("Data/sf_dem.txt", skip = 6, na.strings = "-9999")
matrix_dem = as.matrix(sanfran_dem_sub)
old = matrix_dem

# Analyze elevations between -2 and 8 meters above mean sea level
matrix_dem[matrix_dem < -2] <- NA
matrix_dem[matrix_dem > 8] <- NA

# Plot the elevation cdf as a hypospectic plot. To create the figure in the paper
# this figure will have to be rotated, text and colors need to be added, and the subset
# figure will have to be added in the corner. All of this can be done in powerpoint, gimp, and etc.
png(file="Figures/sanfran_subset2.tiff", family="Times", units="in", width=8, height=8, pointsize=11, res=300)
plot.ecdf(matrix_dem, ylab="", xlab="", main="")
abline(v = 0, lty=2)
abline(v = 1.59, lty=2, col=myheatcolors[3])
abline(v = 2.21, lty=2, col=myheatcolors[5])
abline(v = 2.72, lty=2, col=myheatcolors[9])
axis(side=3, labels=TRUE)
dev.off()

# Estimate and print the fraction of land suspetible to the following water levels.
func = ecdf(matrix_dem)
print(paste("Fraction of land suspetible to 1.59 m: ", round(func(1.59),2)))
print(paste("Fraction of land suspetible to 2.21 m: ", round(func(2.21),2)))
print(paste("Fraction of land suspetible to 2.72 m: ", round(func(2.72),2)))

#################################### END #######################################

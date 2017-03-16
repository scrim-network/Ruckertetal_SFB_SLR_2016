#########################################################################
# file: Plot_DamsReservoirs.R
#------------------------------------------------------------------------
# Author and copyright: Kelsey Ruckert
# Pennsylvania State University
# klr324@psu.edu
# Code written May 2016; edits Feb. 2017
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
# Plot figures for Ruckert et al. (in review)
#
#########################################################################
# Clear global environment
rm(list =ls())

# open packages
library(extRemes)
library(fExtremes)
library(ismev)
library(lubridate)
library(zoo)
library(gridExtra)
library(compiler)
enableJIT(3)
enableJIT(3)

#------------------- Source in scripts -------------------------------------
source("Scripts/plot_sf.r")
source("Scripts/put_fig_letter.r")

#------------------- Create unique colors ----------------------------------
makeTransparent<- function(someColor, alpha=100){
  someColor = someColor
  newColor<-col2rgb(someColor)
  apply(newColor,2 ,
        function(curcoldata)
        {rgb(red=curcoldata[1],
             green=curcoldata[2],
             blue=curcoldata[3], alpha=alpha,
             maxColorValue=255)})
}
someColor <-"gray"
mygray <- makeTransparent(someColor)

#install.packages('RColorBrewer')
library(RColorBrewer)
myheatcolors <- brewer.pal(9,"YlOrRd")
mygreen <- brewer.pal(6,"Greens")

#------------------- Load in the San Francisco tide & sea-level rise workspaces
load("Workspace/SanFranSLR_StormSurge_Analysis_DamRes.RData")

#-------------------------- Set widths and heights ------------------------------
inches_to_dpi = function(inch){ inch * 300 }

text_column_width   = 5.2
minimum_width       = 2.63
full_page_width     = 7.5
full_page_height    = 8.75
single_panel_height = 4

############################# Supp. Figure #8 ##################################
# Plot figure #3 accounting for water stored behind dams and in reservoirs

png(file="Figures/S8_Fig.tiff", family="Times", units="in", width=text_column_width, height=single_panel_height*2, pointsize=11, res=300)
par(mfrow=c(2,1), mgp=c(1.5,.5,0),mar=c(4, 3, 1, 2.5))

pdf2100DamsRes <- density(proj2100DamsRes)

plot(pdf2100DamsRes, col=myheatcolors[9], main="", xlab="Sea-level anomaly in 2100 (m)",
     sub="With respect to the year 2000",ylim=c(0, 2), xlim=c(0,4.5), lwd=1.5, ylab="Probability density",
     yaxt="n") # SLR (MCMC heteroskedastic calib.)

lines(c(mean(proj2100DamsRes),mean(proj2100DamsRes)), c(-1,1.73), col=myheatcolors[5], lty=3,lwd=1.5)
points(mean(proj2100DamsRes), 1.73, pch=21, bg=myheatcolors[5])

lines(c(1.38, 1.38), c(-1,0.95), col=myheatcolors[7], lty=3,lwd=1.5) # Heberger 09 (in paper they round up) Here we will stick to two decimal points
points(1.38, .94, pch=21, bg=myheatcolors[7])

text(mean(proj2100DamsRes), 1.9, "mean SLR", col=myheatcolors[5], cex=0.9)
text((1.4+0.5), 0.9, "Heberger et al.
     (2009) SLR", col=myheatcolors[7], cex=0.9)
text(3.3, 0.2, "full probability density\nusing the simple model", col=myheatcolors[9], cex=0.9)
put.fig.letter("a.",font=2, cex=1)
#----------------------- Return level/ survival function plot ------------------------
par(mgp=c(1.5,.5,0),mar=c(4, 3, 1, 2.5))
#par(mgp=c(1.5,.5,0),mar=c(4, 4, 1, 1))
plot.sf(coredata(year.res.max)/100, pch = 21, bg = "black",
        ylab = "Survival function [1 - cumulative frequency]",
        xlab = "Return level in 2100 (m)", sub="San Francisco Bay",
        yaxt = "n", yaxs = 'i',
        ylim = c(10^-2.5, 10^0+0.25), 
        xlim = c(1, 5.5))

for(i in 1:length(proj2100DamsRes)){
  lines(slr.storm.DR2100[,i], 1-q[1:storm_surgeL], col=mygray, lwd=1)
}

lines(heberger09.2100, 1-q[1:(storm_surgeL)], type="l",col = myheatcolors[7],lwd=1.5)
lines(mean.stormSLR.2100HET, 1-q[1:(storm_surgeL)], type="l",col = myheatcolors[5],lwd=1.5)
lines(fit_q_year/100, 1-q[1:storm_surgeL], type="l",lwd=1.5, col = myheatcolors[3])

#Add in 100-yr value lines
abline(h=0.01, lty=2) # add in the 1 in 100
axis(2, at=10^(-4:-2), label=parse(text=paste("10^", -4:-2, sep="")))
text(4, 0.0115, "100-yr return period", cex=0.85)

axis(4, at=c(10^0, 10^-1, 10^-2), label=c("1", "10", "100"), cex=0.9)
mtext(4, text="Return period (years)", line=1)

# Plot the new probabilities
lines(num.range, average.uncertainty.probs, type="l", col=myheatcolors[9], lwd=1.5)
points(fit_q_year[year100prob]/100, 1-q[year100prob], pch=21, bg= myheatcolors[3])
points(fit_q_year[year100prob]/100, new.ssurge.prob, pch=21, bg= myheatcolors[3])
points(mean.stormSLR.2100HET[year100prob], 1-q[year100prob], pch=21, bg=myheatcolors[5])
points(mean.stormSLR.2100HET[year100prob], new.meanMSLR.prob, pch=21, bg=myheatcolors[5])
points(heberger09.2100[year100prob], 1-q[year100prob], pch=21, bg=myheatcolors[7])
points(heberger09.2100[year100prob], new.average.prob, pch=21, bg=myheatcolors[7])
points(num.range[max.returnL], average.uncertainty.probs[max.returnL], pch=21, bg=myheatcolors[9])

legend("topright",
       c("Observations", "Baseline storm surge +\npotential SLR", "Baseline storm surge",
         "Flood height accounting\nfor mean SLR", "Flood height accounting\nfor Heberger et al. (2009) SLR", "Flood height accounting\nfor SLR uncertainty"),
       col = c("black", "snow2", myheatcolors[3], myheatcolors[5], myheatcolors[7], myheatcolors[9]),
       pch = c(21, NA, 21, 21, 21, 21),
       pt.bg = c("black", NA, myheatcolors[3], myheatcolors[5], myheatcolors[7], myheatcolors[9]),
       lty = c(NA, 1, 1, 1, 1, 1), y.intersp=c(0.9,0.9,0.9,0.9,1,1.1),
       lwd = c(NA, 1.5, 1.5, 1.5, 1.5, 1.5),
       cex=0.8)

put.fig.letter("b.",font=2, cex=1)

dev.off()

#################################### END #######################################
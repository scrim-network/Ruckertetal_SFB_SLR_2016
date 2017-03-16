#########################################################################
# file: Plot_distributions.R
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
load("Workspace/SanFranSLR_StormSurge_Analysis_distribution_test.RData")

#-------------------------- Set widths and heights ------------------------------
inches_to_dpi = function(inch){ inch * 300 }

text_column_width   = 5.2
minimum_width       = 2.63
full_page_width     = 7.5
full_page_height    = 8.75
single_panel_height = 4

############################# Supp. Figure #7 ##################################
# Test different distribution shapes

png(file="Figures/S4_Fig.tif", family="Times", units="in", width=text_column_width, height=single_panel_height*2, pointsize=11, res=300)
par(mfrow=c(2,1), mgp=c(1.5,.5,0),mar=c(4, 3, 1, 2.5))

#pdf2100DamsRes <- density(proj2100DamsRes)
pdf2100 <- density(prob_proj2100/100)

plot(pdf2100, col=myheatcolors[9], main="", xlab="Sea-level anomaly in 2100 (m)",
sub="With respect to the year 2000",ylim=c(0, 6), xlim=c(-0.3,5.6), lwd=2, ylab="Probability density",
yaxt="n") # SLR (MCMC heteroskedastic calib.)

lines(density(pareto_dist_2100), lwd=2, col="royalblue4")
lines(density(lognorm_dist_2100), lwd=2, col="steelblue2")
lines(pdf2100, col=myheatcolors[9], lwd=2)
lines(density(normal_dist_2100), lwd=2, col="powderblue")

legend("topright",c("This study empirical pdf", "Normal dist. approximation", "Log normal dist. approximation", "Pareto dist. approximation"),
col = c(myheatcolors[9], "powderblue", "steelblue2", "royalblue4"),
lty = 1, lwd = 1.5, cex=0.8)

put.fig.letter("a.",font=2, cex=1)

#--Return level/ survival function plot ------------------------
par(mgp=c(1.5,.5,0),mar=c(4, 3, 1, 2.5))
plot.sf(coredata(year.res.max)/100, pch = 21, bg = "black",
ylab = "Survival function [1 - cumulative frequency]",
xlab = "Return level in 2100 (m)",sub="San Francisco Bay",
yaxt = "n", yaxs = 'i',
ylim = c(10^-2.4, 10^0+0.25),
xlim = c(norm.num.range[1], 5))

lines(fit_q_year/100, 1-q[1:storm_surgeL], type="l",lwd=2, col = myheatcolors[3])

abline(h=0.01, lty=2) # add in the 1 in 100
axis(2, at=10^(-4:-2), label=parse(text=paste("10^", -4:-2, sep="")))
text(4, 0.0115, "100-yr return period", cex=0.85)

axis(4, at=c(10^0, 10^-1, 10^-2), label=c("1", "10", "100"), cex=0.9)
mtext(4, text="Return period (years)", line=1)

# Plot the new probabilities
lines(pareto.num.range, pareto.average.uncertainty.probs, type="l", col="royalblue4", lwd=2)
lines(Lnorm.num.range, Lnorm.average.uncertainty.probs, type="l", col="steelblue2", lwd=2)
lines(num.range, average.uncertainty.probs, type="l", col=myheatcolors[9], lwd=2)
lines(norm.num.range, norm.average.uncertainty.probs, type="l", col="powderblue", lwd=2)

#points(fit_q_year[year100prob]/100, 1-q[year100prob], pch=21, bg=myheatcolors[3])
#points(fit_q_year[year100prob]/100, new.ssurge.prob, pch=21, bg=myheatcolors[3])
#points(num.range[max.returnL+1], average.uncertainty.probs[max.returnL+1], pch=21, bg=myheatcolors[9])

legend("topright", c("Observations", "Baseline storm surge", "Flood height accounting for\nthis study empirical pdf",
                     "Flood height accounting for\nnormal dist. approximation", 
                     "Flood height accounting for\nlog normal dist. approximation", 
                     "Flood height accounting for\npareto dist. approximation"),
       col = c("black", myheatcolors[3], myheatcolors[9], "powderblue", "steelblue2", "royalblue4"),
       pch = c(21, NA, NA, NA, NA, NA),
       pt.bg = c("black", NA, NA, NA, NA, NA),
       lty = c(NA, 1, 1, 1, 1, 1), y.intersp=c(0.9,0.9,0.9,1,1.1,1.2),
       lwd = c(NA, 1.5, 1.5, 1.5, 1.5, 1.5),
       cex=0.8)

put.fig.letter("b.",font=2, cex=1)

dev.off()

#################################### END #######################################

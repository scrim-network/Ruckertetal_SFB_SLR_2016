#########################################################################
# file: Plot_Ruckertetal_SFB.R
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
load("Workspace/SanFranSLR_StormSurge_Analysis_norm.RData")

############################# Figure #1 ##################################
# SLR and Tide gauge plot 
datasanfrantide <- read.csv("Data/9414290_MeanSeaLevelTrends.csv", na.strings = '-999')

# Calculate the 90% confidence interval for the heteroskedastic method: 
# Without consideration of Dams and Res.
het_5 <-
  het_50 <-
  het_mean <-
  het_95 <- rep(NA,nyears.mod) # nyears.mod is the total number of data (421 years)
for(i in 1:nyears.mod){
  het_5[i] = quantile(slr.mcmc_proj[,i]/100,0.05) #heteroskedastic SLR values 90%
  het_50[i] = quantile(slr.mcmc_proj[,i]/100,0.5) #median
  het_mean[i] = mean(slr.mcmc_proj[,i]/100) #mean
  het_95[i] = quantile(slr.mcmc_proj[,i]/100,0.95)
}
het_x_90=c(het_5, rev(het_95)); het_y_90=c(years.mod, rev(years.mod))

# #Print the 90% confidence interval, median, and mean for the year 2100 
# (without consideration of Dams and Res:
print(round(median.slr.proj$sle[221]/100,2)) #median based on median parameters
print(round(het_50[221],2)) #median
print(round(mean(prob_proj2100/100),2)) #mean
print(round(het_mean[221],2)) #mean
print(round(het_95[221],2)) #95% quartile
print(round(het_5[221],2)) #5% quartile

#-------------------------- Set widths and heights ------------------------------
inches_to_dpi = function(inch){ inch * 300 }

text_column_width   = 5.2
minimum_width       = 2.63
full_page_width     = 7.5
full_page_height    = 8.75
single_panel_height = 4

#-------------------------- Save plot ------------------------------
png(file="Figures/Fig1.tif", family="Times", units="in", width=text_column_width, height=single_panel_height, pointsize=11, res=300)
par(mfrow=c(1,1), mgp=c(1.5,.5,0), mar=c(3.5,4,1,1))

plot(datasanfrantide$Year[55], datasanfrantide$Monthly_MSL[55], typ="l", 
     ylab="Sea-level anomalies (m)", xlab="Year", xlim=c(1860, 2100),
     ylim=c(-0.32, 1.25))

#polygon(het_y_90, het_x_90, col=mygray, border=NA)
polygon(het_y_90, het_x_90, col="snow2", border=NA)
points(datasanfrantide[55:1991,9], datasanfrantide$Monthly_MSL[55:1991], pch=20, cex=1/72)
lines(years.mod[122:421], het_mean[122:421], col="deepskyblue4", lwd=2)
lines(year, slr/100, col=mygreen[4], lwd=2)
# lines(years.mod, median.slr.proj$sle/100, col="blue")

abline(v=2002, lty=2, col=mygreen[4])
text(1970, 1.2, "Observations used\nfor global mean sea-\nlevel anomalies", cex=0.7, col=mygreen[4])
arrows(2001, 0.8, 1981, 0.8, length = 0.1, lty=1, col = mygreen[4])

text(2040, 1.1, cex=0.7, "Projection mean and 90%\ncredible interval\nfor the global mean\nsea-level anomaly\nusing the simple\nillustrative model", col="deepskyblue4")
arrows(2003, 0.8, 2021, 0.8, length = 0.1, lty=1, col = "deepskyblue4")

axis(side=4, labels=FALSE)
text(1950, 0.3, cex=0.7, "Global sea-level anomalies\n(Church and White 2006)", col=mygreen[4])
text(1880, 0.15, cex=0.7, "San Francisco Bay\nTide gauge\n9414290", col="black")
box(lwd = 1)
dev.off()

############################# Figure #2 ##################################
### Block Maxima plots
png(file="Figures/Fig2.tif", family="Times", units="in", width=text_column_width, height=single_panel_height*2, pointsize=11, res=300)
par(mfrow=c(2,1), mgp=c(1.5,.5,0),mar=c(4, 3, 1, 2))
plot(index(year.res.max), coredata(year.res.max)/100, pch=20, type = "b", col = "black",#, xaxs = 'i',
ylab = "Annual block maxima (m)", xlab = "Year", lwd=1.5)
box(lwd = 1)
put.fig.letter("a.",font=2, cex=1)

# Return level plot function (written by Perry Oddo)
return_level_plot <- function(block_maxima, max_return_period, legend)
{
    # Data (block_maxima) should be a numeric vector of block maxima tide observations
    # max_return_period should be the maximum number of years returned on plot
    # e.g. 10^4 for a 1/10,000 year return level)
    # legend = TRUE to plot legend
    
    require(ismev)
    fit.obj <- gev.fit(block_maxima, show = FALSE)
    
    a <- fit.obj$mle
    mat <- fit.obj$cov
    dat <- fit.obj$data
    
    eps <- 1e-06; a1 <- a; a2 <- a; a3 <- a
    a1[1] <- a[1] + eps; a2[2] <- a[2] + eps; a3[3] <- a[3] + eps
    f <- c(seq(0.01, 0.09, by = 0.01), 0.1, 0.2, 0.3, 0.4, 0.5,
    0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.995, 0.999, 0.9999, 0.99999, 0.999999)
    q <- gevq(a, 1 - f[which((1/(1-f)) <= (max_return_period+0.01))])
    d <- t(gev.rl.gradient(a = a, p = 1 - f))
    v <- apply(d, 1, q.form, m = mat)
    plot(1/(1-(f[1:length(q)])), q, log = "x", type = "n", xlim = c(0.8, 10^3),
         ylim = c(min(dat, q), 1.7), #c(min(dat, q), max(dat, q + 1.96 * sqrt(v[1:length(q)]))), 
         xaxt = 'n', cex=1,
         xlab = "Return period (years)", 
         ylab = "Return level (m)")
    axis(1, lwd = 1, at=10^(seq(-1,log10(10^3), by = 1)), label=c(0.1, 1, 10, 100, 1000))
    axis(2, lwd = 1)
    lines(1/(1-(f[1:length(q)])), q, lty = 1, lwd = 2)
    # lines(1/(1-(f[1:length(q)])), q + 1.96 * sqrt(v[1:length(q)]), col = "#0080FFFF", lwd = 1.5)
    # lines(1/(1-(f[1:length(q)])), q - 1.96 * sqrt(v[1:length(q)]), col = "#0080FFFF", lwd = 1.5)
    points(1/(1-(1:length(dat)/(length(dat) + 1))), sort(dat), lwd = 1, cex = 0.75, pch = 21, bg = "white")
    box(lwd = 1)
    
    if(legend == TRUE | legend == T)
    {
      legend("topleft",
             c("Annual block maxima observations", "Best estimate"),
             col = c("black", "black"),
             pt.bg = c("white", NA),
             pch = c(21, NA),
             lty = c(NA, 1), cex=1,
             lwd = c(1.5, 1.5),
             bty = 'n',
             inset = c(0.01, -0.01))
    }
    
}

par(mgp=c(1.5,.5,0),mar=c(4, 3, 1, 2))
return_level_plot(coredata(year.res.max)/100, 10^4, legend = T)
put.fig.letter("b.",font=2, cex=1)
dev.off()

############################# Figure #3 ##################################
#Plot pdf of sea-level

#-- Probability density plot ------------------------
png(file="Figures/Fig3.tif", family="Times", units="in", width=text_column_width, height=single_panel_height*2, pointsize=11, res=300)
par(mfrow=c(2,1), mgp=c(1.5,.5,0),mar=c(4, 3, 1, 2.5))

#pdf2100DamsRes <- density(proj2100DamsRes)
pdf2100 <- density(prob_proj2100/100)

plot(pdf2100, col=myheatcolors[9], main="", xlab="Sea-level anomaly in 2100 (m)",
sub="With respect to the year 2000",ylim=c(0, 2), xlim=c(0,3.1), lwd=1.5, ylab="Probability density",
     yaxt="n") # SLR (MCMC heteroskedastic calib.)

lines(c(mean(prob_proj2100/100),mean(prob_proj2100/100)), c(-1,1.73), col=myheatcolors[5], lty=3,lwd=1.5)
points(mean(prob_proj2100/100), 1.73, pch=21, bg=myheatcolors[5])

lines(c(a2new[45],a2new[45]), c(-1,0.95), col=myheatcolors[7], lty=3,lwd=1.5)
points(a2new[45], .94, pch=21, bg=myheatcolors[7])

text(mean(prob_proj2100/100), 1.9, "mean SLR", col=myheatcolors[5], cex=0.9)
text((0.8+0.5), 0.9, "Heberger et al.
(2009) SLR", col=myheatcolors[7], cex=0.9)
text(2.3, 0.2, "full probability density\nusing the simple model", col=myheatcolors[9], cex=0.9)
put.fig.letter("a.",font=2, cex=1)

#--Return level/ survival function plot ------------------------
par(mgp=c(1.5,.5,0),mar=c(4, 3, 1, 2.5))
plot.sf(coredata(year.res.max)/100, pch = 21, bg = "black",
        ylab = "Survival function [1 - cumulative frequency]",
        xlab = "Return level in 2100 (m)", sub="San Francisco Bay",
        yaxt = "n", yaxs = 'i',
        ylim = c(10^-2.5, 10^0+0.25), 
        xlim = c((fit_q_year[2]/100), 5))

for(i in 1:length(prob_proj2100)){
  lines(slr.storm.DR2100[,i], 1-q[1:storm_surgeL], col=mygray, lwd=1)
}

lines(heberger09.2100, 1-q[1:(storm_surgeL)], type="l",col = myheatcolors[7],lwd=1.5)
lines(mean.stormSLR.2100HET, 1-q[1:(storm_surgeL)], type="l",col = myheatcolors[5],lwd=1.5)
lines(fit_q_year/100, 1-q[1:storm_surgeL], type="l",lwd=1.5, col = myheatcolors[3])

abline(h=0.01, lty=2) # add in the 1 in 100
axis(2, at=10^(-4:-2), label=parse(text=paste("10^", -4:-2, sep="")))
text(4, 0.0115, "100-yr return period", cex=0.85)

axis(4, at=c(10^0, 10^-1, 10^-2), label=c("1", "10", "100"), cex=0.9)
mtext(4, text="Return period (years)", line=1)

# Plot the new probabilities
lines(num.range, average.uncertainty.probs, type="l", col=myheatcolors[9], lwd=1.5)
points(fit_q_year[year100prob]/100, 1-q[year100prob], pch=21, bg=myheatcolors[3])
points(fit_q_year[year100prob]/100, new.ssurge.prob, pch=21, bg=myheatcolors[3])
points(mean.stormSLR.2100HET[year100prob], 1-q[year100prob], pch=21, bg=myheatcolors[5])
points(mean.stormSLR.2100HET[year100prob], new.meanMSLR.prob, pch=21, bg=myheatcolors[5])
points(heberger09.2100[year100prob], 1-q[year100prob], pch=21, bg=myheatcolors[7])
points(heberger09.2100[year100prob], new.average.prob, pch=21, bg=myheatcolors[7])
points(num.range[max.returnL+1], average.uncertainty.probs[max.returnL+1], pch=21, bg=myheatcolors[9])

legend("topright", c("Observations", "Baseline storm surge +\npotential SLR", "Baseline storm surge",
                     "Flood height accounting\nfor mean SLR", "Flood height accounting\nfor Heberger et al. (2009) SLR", 
                     "Flood height accounting\nfor SLR uncertainty"),
       col = c("black", "snow2", myheatcolors[3], myheatcolors[5], myheatcolors[7], myheatcolors[9]),
       pch = c(21, NA, 21, 21, 21, 21),
       pt.bg = c("black", NA, myheatcolors[3], myheatcolors[5], myheatcolors[7], myheatcolors[9]),
       lty = c(NA, 1, 1, 1, 1, 1), y.intersp=c(0.9,0.9,0.9,0.9,1,1.1),
       lwd = c(NA, 1.5, 1.5, 1.5, 1.5, 1.5),
       cex=0.8)

put.fig.letter("b.",font=2, cex=1)

dev.off()

############################# Figure #5 ########################################

# Figure 5 is created in ArcGIS ArcMap. Refer to the tutorial on how this is done.

############################# Supp. Figure #1 ##################################

# Supp. Figure 1 is a photograph.

############################# Supp. Figure #3 ##################################
#----------------------- Code Checking Plots  ----------------------------------
# This plot displays why the probability occuring increases when accounting for uncertainty
# General the sea-level rise with +/- 1 sigma
sigma1 = sd(prob_proj2100/100) #estimate the standard deviation of SLR
plus1sig = mean(prob_proj2100/100) + sigma1
minus1sig = mean(prob_proj2100/100) - sigma1

#--- GEVfit to 10^5
Psig1.ss.slr <- fit_q_year/100 + plus1sig
Msig1.ss.slr <- fit_q_year/100 + minus1sig

new.prob.plussig = inv.sf(Psig1.ss.slr, mean.stormSLR.2100HET[year100prob])
new.prob.minussig = inv.sf(Msig1.ss.slr, mean.stormSLR.2100HET[year100prob])
new.prob.minussig <- 0 #if NA is returned

average.sig.prob = mean(c(new.prob.plussig, 1-q[year100prob], new.prob.minussig))

#--- Generate GEVfit to 10^6
SSL <- 10^6 # desired storm surge level
quant = seq(0,1,length.out= SSL +1)  # quantile array

# Find closed-form solution of GEV fit
SSfit_q_year = qgev(quant, year.res.max.fit2@fit$par.ests[1], year.res.max.fit2@fit$par.ests[2],
                    year.res.max.fit2@fit$par.ests[3])
SSfit_q_year = SSfit_q_year[SSfit_q_year< max(SSfit_q_year)]

# Find which q is the 100-yr value 
my.num = which(quant <= 0.99)
my.num.max = which.max(my.num)
msig.year100prob <- my.num.max
# check to make sure the msig.year100prob value is the 100-yr value (10^-2)
print(1-quant[msig.year100prob])

Msig1.ss.slr.mil <- SSfit_q_year/100 + minus1sig
Psig1.ss.slr.mil <- SSfit_q_year/100 + plus1sig
mean.ss.slr.mil <- SSfit_q_year/100 + mean(prob_proj2100/100)

prob.mil.minussig = inv.sf(Msig1.ss.slr.mil, mean.stormSLR.2100HET[year100prob])

real.average.sig.prob = mean(c(new.prob.plussig, 1-q[year100prob], prob.mil.minussig))
#-------------------------------------------

png(file="Figures/S3_Fig.tif", family="Times", units="in", width=text_column_width, height=full_page_height, pointsize=13, res=300)
par(mfrow=c(3,1), mgp=c(1.5,.5,0),mar=c(4, 3, 1, 2.5))

plot(pdf2100, col="black", main="", lwd=1.5, xlab="Sea-level anomaly in 2100 (m)",
sub="With respect to the year 2000",ylim=c(0, 2), xlim=c(-0.1, 3.1), ylab="Probability density",
     yaxt="n") # SLR (MCMC heteroskedastic calib.)

lines(c(mean(prob_proj2100/100),mean(prob_proj2100/100)), c(-1,1.73), col=mygreen[4], lty=3,lwd=1.5)
lines(c(minus1sig,minus1sig), c(-1,1.28), col=mygreen[2], lty=3,lwd=1.5)
lines(c(plus1sig,plus1sig), c(-1,0.8), col=mygreen[6], lty=3,lwd=1.5)

points(minus1sig, 1.28, pch=21, bg=mygreen[2])
points(mean(prob_proj2100/100), 1.73, pch=21, bg=mygreen[4])
points(plus1sig, .8, pch=21, bg=mygreen[6])

text(mean(prob_proj2100/100), 1.9, "mean", col=mygreen[4], cex=1)
text((minus1sig - .3), 1.25, "-1 standard
deviation", col=mygreen[2], cex=1)
text((plus1sig + .3), 0.8, "+1 standard
deviation", col=mygreen[6], cex=1)
text(2.3, 0.2, "full probability density
using the simple model", col="black", cex=1)

put.fig.letter("a.",font=2, cex=1)

#-- Return level/ survival function plot ------------------------
par(mgp=c(1.5,.5,0),mar=c(4, 3, 1, 2.5))
plot.sf(coredata(year.res.max)/100, pch = 21, bg = "black",
        ylab = "Survival function [1 - cumulative frequency]",
        xlab = "Return level in 2100 (m)",sub="San Francisco Bay",
        yaxt = "n", yaxs = 'i',
        ylim = c(10^-6, 10^0+1), 
        xlim = c(min(mean.stormSLR.2100HET), 4.5))

abline(v=mean.stormSLR.2100HET[year100prob], lty=2, col=mygreen[4])
abline(v=Psig1.ss.slr[year100prob], lty=2, col=mygreen[6])
abline(v=Msig1.ss.slr.mil[msig.year100prob], lty=2, col=mygreen[2])
axis(2, at=10^(-6:-3), label=parse(text=paste("10^", -6:-3, sep="")))

axis(4, at=c(10^0, 10^-1, 10^-2, 10^-3, 10^-4, 10^-5), label=c("1", "10", "100", "1000", "10000", "100000"), cex=0.9)
mtext(4, text="Return period (years)", line=1, cex=0.65)

lines(SSfit_q_year/100, 1-quant[1:SSL], type="l", col = "skyblue", lwd=1.5)
lines(mean.ss.slr.mil, 1-quant[1:SSL], type="l",col=mygreen[4], lwd=1.5)
lines(Psig1.ss.slr.mil , 1-quant[1:SSL], type="l",col=mygreen[6], lwd=1.5)
lines(Msig1.ss.slr.mil, 1-quant[1:SSL], type="l",col=mygreen[2], lwd=1.5)
abline(h=0.01, lty=2) # add in the 1 in 100
text(4, 0.014, "100-yr return period", cex=0.85)

points(mean.stormSLR.2100HET[year100prob], new.prob.plussig, pch=21, bg=mygreen[6])
points(mean.stormSLR.2100HET[year100prob], prob.mil.minussig, pch=21, bg=mygreen[2])
points(mean.stormSLR.2100HET[year100prob], 1-q[year100prob], pch=21, bg=mygreen[4])
points(mean.stormSLR.2100HET[year100prob], real.average.sig.prob , pch=15, col="black", cex=1.5)
points(mean.stormSLR.2100HET[year100prob], average.sig.prob , pch=8, col="pink")
put.fig.letter("b.",font=2, cex=1)

legend("bottomright", ncol=1,
       c("Observations", "Baseline storm surge", "Baseline storm surge +\nmean - 1 standard deviation SLR",
         "Baseline storm surge + mean SLR", "Baseline storm surge +\nmean + 1 standard deviation SLR", 
         "Averaged 100-yr return period", "Approximated 100-yr return period"),
       col = c("black", "skyblue", mygreen[2], mygreen[4], mygreen[6], "black", "pink"),
       pch = c(21, NA, 21, 21, 21, 15, 8),
       pt.bg = c("black", NA, mygreen[2], mygreen[4], mygreen[6], NA, NA),
       lty = c(NA, 1, 1, 1, 1, NA, NA),
       #lwd = c(NA, 1.5, 1.5, 1.5, 1.5, NA, NA),
       y.intersp=c(0.9,1,1,1,1,1,0.9),
       lwd = 1.5,
       bty = 'n' , cex=0.85)

par(mgp=c(1.5,.5,0),mar=c(4, 3, 1, 2.5))
return.period = c(round(1-q[year100prob],2), round(average.sig.prob,6), round(real.average.sig.prob,6))
barplot(return.period, xlab="100-yr return period", ylab="Return period at 2.2m (decimals)",
        names.arg=c("Mean 100-yr return period", "Approximated", "Averaged 100-yr return period"), col=c(mygreen[4], "pink", "black"))
put.fig.letter("c.",font=2, cex=1)

dev.off()

############################# Supp. Figure #2 ##################################
# Add the 90% quantile to the figure #3

#-- Probability density plot ------------------------
png(file="Figures/S2_Fig.tif", family="Times", units="in", width=text_column_width, height=single_panel_height*2, pointsize=11, res=300)
par(mfrow=c(2,1), mgp=c(1.5,.5,0),mar=c(4, 3, 1, 2.5))

#pdf2100DamsRes <- density(proj2100DamsRes)
pdf2100 <- density(prob_proj2100/100)

plot(pdf2100, col=myheatcolors[9], main="", xlab="Sea-level anomaly in 2100 (m)",
sub="With respect to the year 2000",ylim=c(0, 2), xlim=c(0,3.1), lwd=1.5, ylab="Probability density",
yaxt="n") # SLR (MCMC heteroskedastic calib.)

lines(c(mean(prob_proj2100/100),mean(prob_proj2100/100)), c(-1,1.73), col=myheatcolors[5], lty=3,lwd=1.5)
points(mean(prob_proj2100/100), 1.73, pch=21, bg=myheatcolors[5])

lines(c(a2new[45],a2new[45]), c(-1,0.95), col=myheatcolors[7], lty=3,lwd=1.5)
points(a2new[45], .94, pch=21, bg=myheatcolors[7])

lines(c(quantile(prob_proj2100/100, 0.90), quantile(prob_proj2100/100, 0.90)), c(-1,0.6), col=myheatcolors[8], lty=3,lwd=1.5)
points(quantile(prob_proj2100/100, 0.90), .6, pch=21, bg=myheatcolors[8])

text(mean(prob_proj2100/100), 1.9, "mean SLR", col=myheatcolors[5], cex=0.9)
text((0.8+0.5), 0.9, "Heberger et al.\n(2009) SLR", col=myheatcolors[7], cex=0.9)
text(2.3, 0.2, "full probability density\nusing the simple model", col=myheatcolors[9], cex=0.9)
text((quantile(prob_proj2100/100, 0.90)+0.5), 0.6, "90% quantile SLR", col=myheatcolors[8], cex=0.9)

put.fig.letter("a.",font=2, cex=1)

#--Return level/ survival function plot ------------------------
par(mgp=c(1.5,.5,0),mar=c(4, 3, 1, 2.5))
plot.sf(coredata(year.res.max)/100, pch = 21, bg = "black",
ylab = "Survival function [1 - cumulative frequency]",
xlab = "Return level in 2100 (m)",sub="San Francisco Bay",
yaxt = "n", yaxs = 'i',
ylim = c(10^-2.5, 10^0+0.25),
xlim = c((fit_q_year[2]/100), 5))

for(i in 1:length(prob_proj2100)){
    lines(slr.storm.DR2100[,i], 1-q[1:storm_surgeL], col=mygray, lwd=1)
}

lines(heberger09.2100, 1-q[1:(storm_surgeL)], type="l",col = myheatcolors[7],lwd=1.5)
lines(mean.stormSLR.2100HET, 1-q[1:(storm_surgeL)], type="l",col = myheatcolors[5],lwd=1.5)
lines(quant90.stormSLR.2100HET, 1-q[1:(storm_surgeL)], type="l",col = myheatcolors[8],lwd=1.5)
lines(fit_q_year/100, 1-q[1:storm_surgeL], type="l",lwd=1.5, col = myheatcolors[3])

abline(h=0.01, lty=2) # add in the 1 in 100
axis(2, at=10^(-4:-2), label=parse(text=paste("10^", -4:-2, sep="")))
text(4, 0.0115, "100-yr return period", cex=0.85)

axis(4, at=c(10^0, 10^-1, 10^-2), label=c("1", "10", "100"), cex=0.9)
mtext(4, text="Return period (years)", line=1)

# Plot the new probabilities
lines(num.range, average.uncertainty.probs, type="l", col=myheatcolors[9], lwd=1.5)
points(fit_q_year[year100prob]/100, 1-q[year100prob], pch=21, bg=myheatcolors[3])
points(fit_q_year[year100prob]/100, new.ssurge.prob, pch=21, bg=myheatcolors[3])
points(mean.stormSLR.2100HET[year100prob], 1-q[year100prob], pch=21, bg=myheatcolors[5])
points(mean.stormSLR.2100HET[year100prob], new.meanMSLR.prob, pch=21, bg=myheatcolors[5])
points(quant90.stormSLR.2100HET[year100prob], 1-q[year100prob], pch=21, bg=myheatcolors[8])
points(quant90.stormSLR.2100HET[year100prob], new.quant90MSLR.prob, pch=21, bg=myheatcolors[8])
points(heberger09.2100[year100prob], 1-q[year100prob], pch=21, bg=myheatcolors[7])
points(heberger09.2100[year100prob], new.average.prob, pch=21, bg=myheatcolors[7])
points(num.range[max.returnL+1], average.uncertainty.probs[max.returnL+1], pch=21, bg=myheatcolors[9])

legend("topright", c("Observations", "Baseline storm surge +\npotential SLR", "Baseline storm surge",
                     "Flood height accounting\nfor mean SLR", "Flood height accounting\nfor Heberger et al. (2009) SLR", 
                     "Flood height accounting\nfor 90% quantile SLR", "Flood height accounting\nfor SLR uncertainty"),
       col = c("black", "snow2", myheatcolors[3], myheatcolors[5], myheatcolors[7], myheatcolors[8], myheatcolors[9]),
       pch = c(21, NA, 21, 21, 21, 21, 21),
       pt.bg = c("black", NA, myheatcolors[3], myheatcolors[5], myheatcolors[7], myheatcolors[8], myheatcolors[9]),
       lty = c(NA, 1, 1, 1, 1, 1, 1), y.intersp=c(0.9,0.9,0.9,0.9,1,1.1,1.2),
       lwd = c(NA, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5),
       cex=0.8)

put.fig.letter("b.",font=2, cex=1)

dev.off()

################################ PRINT NUMBERS #################################
print(paste("Print the 90% confidence interval, median, and mean for the year 2100"))
# (without consideration of Dams and Res:
print(paste("Median SLR = ",round(median.slr.proj$sle[221]/100,2))) #median based on median parameters
print(paste("Mean SLR = ",round(mean(prob_proj2100/100),2))) #mean
print(paste("95% quartile SLR = ",round(het_95[221],2))) #95% quartile
print(paste("5% quartile SLR = ",round(het_5[221],2))) #5% quartile

# Print the value to account for changes in Dams and Reserviors in California:
#round(DamsReservoirs,2)

print(paste("Current 100-yr storm surge values:"))
print(paste("Current storm surge for San Fran. = ",round(fit_q_year[year100prob]/100,2)))
print(paste("Current storm surge + Mean SLR = ",round(mean.stormSLR.2100HET[year100prob],2)))
print(paste("Heberger et al 09 storm surge w/o changes in land water storage = ",round(heberger09.2100[year100prob],2)))
print(paste("Future storm surge accounting for uncertain SLR = ",round(num.range[max.returnL+1],2)))

# using the 90% quantile
print(paste("Current storm surge + 90% quantile SLR = ",round(quant90.stormSLR.2100HET[year100prob],2)))

print(paste("The range of return level values used for accounting for SLR uncertainty:"))
print(paste("From = ", round(num.range[1], 2), " m to = ", round(num.range[length(num.range)],2), " m"))

print(paste("Factor of increase:"))
print(paste("From Heberger et al 2009 to new estimate:"))
print(paste("0.01 to: ", round(new.average.prob, 2)))

print(paste("From mean slr estimate to new estimate:"))
print(paste("0.01 to: ", round(new.meanMSLR.prob, 2)))

print(paste("From current storm surge (w/o slr) estimate to new estimate:"))
print(paste("0.01 to: ", round(new.ssurge.prob, 2)))

# using the 90% quantile
print(paste("From the 90% quantile slr estimate to new estimate:"))
print(paste("0.01 to: ", round(new.quant90MSLR.prob, 2)))

print(paste("100-yr Storm surge values and probabilities for the mean +/- sigma test"))
print(paste("Current storm surge + Mean SLR = ",round(mean.stormSLR.2100HET[year100prob],2)))
print(paste("Current storm surge + mean - 1 sigma = ",round(Msig1.ss.slr.mil[msig.year100prob],2)))
print(paste("Current storm surge + mean + 1 sigma = ",round(Psig1.ss.slr[year100prob],2)))

print(paste("Probability of occurence at the mean 100-yr storm surge return level"))
print(paste("mean = ", round(1-q[year100prob], 2), " ; should be 0.01"))
print(paste("mean - 1 sigma = ",round(prob.mil.minussig, 6))) # minus 1 sigma
print(paste("mean + 1 sigma = ",round(new.prob.plussig, 2))) # plus 1 sigma

print(paste("Factor of increase by accounting for uncertainty in simple example:"))
print(paste("From mean slr estimate to new estimate:"))
print(paste("0.01 to: ", round(real.average.sig.prob, 2)))

#################################### END #######################################



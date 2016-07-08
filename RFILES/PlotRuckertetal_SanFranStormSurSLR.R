#########################################################################
# file: PlotRuckertetal_SanFranStormSurSLR.R
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
#------------------------------------------------------------------------
# Plot figures for Ruckert et al. (in prep)
# Figure 1-3 and Supp. Fig. 2 & 6
# 
#########################################################################

rm(list =ls()) #Clear global environment
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
source("Scripts/plot_rangefn.R")
source("Scripts/plot_sf.r")
source("Scripts/put_fig_letter.r")
source("Scripts/return_level_plot.R")

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
# 100-yr storm surge data
load("Workspace/SanFranSLR_StormSurge_Analysis.RData")

############################# Figure #1 ##################################
# SLR and Tide guage plot 
# The data is the mean for each month from 1784-2013
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

#-------------------------- Save pdf plot ------------------------------
#pdf(file="Figures/Fig_1.pdf", family="Helvetica", width=5.2, height=4, pointsize=11)
png(file="Figures/Fig_1.tif", family="Helvetica", units="in", width=5.2, height=4, pointsize=11, res=300)
par(mfrow=c(1,1), mgp=c(1.5,.5,0),mar=c(4, 3, 3, 3.5))
plot(datasanfrantide$Year[55], datasanfrantide$Monthly_MSL[55], typ="l", 
     ylab="Sea-Level anomalies (m)", xlab="Year", xlim=c(1860, 2100),
     ylim=c(-0.32, 1.25))

polygon(het_y_90, het_x_90, col=mygray, border=NA)
lines(datasanfrantide[55:1991,9], datasanfrantide$Monthly_MSL[55:1991])
lines(years.mod[122:421], het_mean[122:421], col=mygreen[4], lwd=2)
lines(year, slr/100, col=mygreen[4], lwd=2)
# lines(years.mod, median.slr.proj$sle/100, col="blue")

abline(v=2002, lty=2, col=mygreen[4])
text(1980, 1.2, "Observed", cex=0.7)
text(2050, 1.2, cex=0.7, "Projection mean and 90% 
credible interval")
axis(side=4, labels=FALSE)
text(1950, 0.3, cex=0.7, "Church and White (2006)
Global mean SLR", col=mygreen[4])
text(1880, 0.3, cex=0.7, "San Francisco Bay
Tide gauge
9414290", col="black")
dev.off()

############################# Figure #2 ##################################
### Block Maxima plots
png(file="Figures/Fig_2_col.tif", family="Helvetica", units="in", width=5, height=8.6, pointsize=12, res=300)
#pdf(file="Figures/Fig_2_col.pdf", family="Helvetica", width=5, height=8.6, pointsize=12)
par(mfrow=c(2,1), mgp=c(1.5,.5,0),mar=c(4, 3, 3.5, 2))
plot(index(year.res.max), coredata(year.res.max)/100, type = 'h', col = "skyblue", xaxs = 'i',
     ylab = "Annual block maxima (m)", xlab = "Year", lwd=1.5)
box(lwd = 1)
put.fig.letter("a.",font=2, cex=1)


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
    q <- gevq(a, 1 - f[which((ceiling(-1/log(f)) - max_return_period <= 0))])
    d <- t(gev.rl.gradient(a = a, p = 1 - f))
    v <- apply(d, 1, q.form, m = mat)
    plot(-1/log(f[1:length(q)]), q, log = "x", type = "n", xlim = c(0.1, 10^3),
    ylim = c(min(dat, q), 1.85),#max(dat, q + 1.96 * sqrt(v[1:length(q)]))),
    xaxt = 'n', cex=1, xlab = "Return Period (years)", ylab = "Return Level (m)")
    axis(1, lwd = 1, at=10^(seq(-1,log10(10^3), by = 1)), label=c(0.1, 1, 10, 100, 1000))
    axis(2, lwd = 1)
    lines(-1/log(f[1:length(q)]), q, lty = 2, lwd = 1.5)
    lines(-1/log(f[1:length(q)]), q + 1.96 * sqrt(v[1:length(q)]), col = "skyblue", lwd = 1.5)
    lines(-1/log(f[1:length(q)]), q - 1.96 * sqrt(v[1:length(q)]), col = "skyblue", lwd = 1.5)
    points(-1/log((1:length(dat))/(length(dat) + 1)), sort(dat), lwd = 1, cex = 0.75, pch = 21, bg = "white")
    box(lwd = 1)
    
    if(legend == TRUE | legend == T)
    {
        legend("topleft",
        c("Annual block maxima observations", "95% confidence interval", "Best estimate"),
        col = c("black", "skyblue", "black"),
        pt.bg = c("white", NA, NA),
        pch = c(21, NA, NA),
        lty = c(NA, 1, 2), cex=0.8,
        lwd = c(1.5, 1.5, 1.5),
        bty = 'n',
        inset = c(0.01, -0.01))
    }
    
}

par(mgp=c(1.5,.5,0),mar=c(4, 3, 3.5, 2))
return_level_plot(coredata(year.res.max)/100, 10^4, legend = T)
put.fig.letter("b.",font=2, cex=1)
dev.off()

############################# Figure #3 ##################################
#Plot pdf of sea-level plus accounting for runoff changes versus 
# the current SLR accounting for runoff in Heberger 09
# source("Scripts/plot_rangefn.R")

#-- Probability density plot ------------------------
png(file="Figures/Fig_3_col.tif", family="Helvetica", units="in", width=5, height=8.6, pointsize=12, res=300)
par(mfrow=c(2,1), mgp=c(1.5,.5,0),mar=c(4, 3, 3.5, 2))

pdf2100DamsRes <- density(proj2100DamsRes)

plot(pdf2100DamsRes, col=myheatcolors[9], main="", xlab="Sea-level anomaly in 2100 (m)",
sub="With respect to the year 2000",ylim=c(0, 2.25), xlim=c(0,4.5), lwd=1.5, ylab="Probability density",
     yaxt="n") # SLR (MCMC heteroskedastic calib.) accounting for runoff changes

lines(c(mean(proj2100DamsRes),mean(proj2100DamsRes)), c(-1,1.73), col=myheatcolors[5], lty=3,lwd=1.5)
points(mean(proj2100DamsRes), 1.73, pch=21, bg=myheatcolors[5])

lines(c(1.38,1.38), c(-1,0.95), col=myheatcolors[7], lty=3,lwd=1.5)
points(1.38, .94, pch=21, bg=myheatcolors[7])

quant = c(0.05, 0.95)
add.hor.box(proj2100DamsRes, quant, width.size = 0.2, where.at = 2.1, tick.length = 0.05, line.width = 1.5, color = myheatcolors[9])
put.fig.letter("a.",font=2, cex=1)

legend("topright", ncol=1,
c("Observations", "Current storm surge", "Storm surge + potential SLR",
"Storm surge + mean SLR", "Storm surge + Heberger
et al. (2009) SLR", "Storm surge + SLR uncertainty"),
col = c("black", myheatcolors[3], mygray, myheatcolors[5], myheatcolors[7], myheatcolors[9]),
pch = c(21, 21, NA, 21, 21, 21),
pt.bg = c("black", myheatcolors[3], NA, myheatcolors[5], myheatcolors[7], myheatcolors[9]),
lty = c(NA, 1, 1, 1, 1, 1), y.intersp=c(0.9,0.9,0.9,0.9,1,1.1),
lwd = c(NA, 1.5, 1.5, 1.5, 1.5, 1.5),
bty = 'n' , cex=0.8)

#--Return level/ survival function plot ------------------------
par(mgp=c(1.5,.5,0),mar=c(4, 3, 3.5, 2))
plot.sf(coredata(year.res.max)/100, pch = 21, bg = "black",
        ylab = "Survival function [1 - cdf]",
        xlab = "Return Level (m)",sub="San Francisco Bay",
        yaxt = "n", yaxs = 'i',
        ylim = c(10^-2.5, 10^0+0.25), 
        xlim = c((fit_q_year[2]/100), 5.0))

for(i in 1:length(proj2100DamsRes)){
  lines(slr.storm.DR2100[,i], 1-q[1:storm_surgeL], col=mygray, lwd=1)
}

lines(heberger09.2100, 1-q[1:(storm_surgeL)], type="l",col = myheatcolors[7],lwd=1.5)
lines(mean.stormSLR.2100HET, 1-q[1:(storm_surgeL)], type="l",col = myheatcolors[5],lwd=1.5)
lines(fit_q_year/100, 1-q[1:storm_surgeL], type="l",lwd=1.5, col = myheatcolors[3])

abline(h=0.01, lty=2) # add in the 1 in 100
axis(2, at=10^(-4:-2), label=parse(text=paste("10^", -4:-2, sep="")))
text(4.5, 0.015, "1:100 level", cex=1)

# Plot the new probabilities
lines(num.range, average.uncertainty.probs, type="l", col=myheatcolors[9], lwd=1.5)
points(fit_q_year[year100prob]/100, 1-q[year100prob], pch=21, bg=myheatcolors[3])
points(fit_q_year[year100prob]/100, new.ssurge.prob, pch=21, bg=myheatcolors[3])
points(mean.stormSLR.2100HET[year100prob], 1-q[year100prob], pch=21, bg=myheatcolors[5])
points(mean.stormSLR.2100HET[year100prob], new.meanMSLR.prob, pch=21, bg=myheatcolors[5])
points(heberger09.2100[year100prob], 1-q[year100prob], pch=21, bg=myheatcolors[7])
points(heberger09.2100[year100prob], new.average.prob, pch=21, bg=myheatcolors[7])
points(num.range[max.returnL], average.uncertainty.probs[max.returnL], pch=21, bg=myheatcolors[9])
put.fig.letter("b.",font=2, cex=1)

dev.off()

############################# Figure #4 ########################################

# Figure 4 is created in ArcGIS ArcMap. Refer to the tutorial on how this is done.

############################# Supp. Figure #1 ##################################

# Supp. Figure 1 is a photograph.

############################# Supp. Figure #2 ##################################
#----------------------- Code Checking Plots  ----------------------------------
# This plot displays why the probability occuring increases when accounting for uncertainty
# General the sea-level rise with +/- 1 sigma
sigma1 = sd(proj2100DamsRes) #estimate the standard deviation of SLR + runoff changes
plus1sig = mean(proj2100DamsRes) + sigma1
minus1sig = mean(proj2100DamsRes) - sigma1

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
mean.ss.slr.mil <- SSfit_q_year/100 + mean(proj2100DamsRes)

prob.mil.minussig = inv.sf(Msig1.ss.slr.mil, mean.stormSLR.2100HET[year100prob])

real.average.sig.prob = mean(c(new.prob.plussig, 1-q[year100prob], prob.mil.minussig))
#-------------------------------------------

png(file="Figures/SuppFigures/S2_Fig.tif", family="Helvetica", units="in", width=5, height=8.6, pointsize=12, res=300)
par(mfrow=c(2,1), mgp=c(1.5,.5,0),mar=c(4, 3, 3.5, 2))

plot(pdf2100DamsRes, col="black", main="", lwd=1.5, xlab="Sea-level anomaly in 2100 (m)", 
sub="With respect to the year 2000",ylim=c(0, 2.5), xlim=c(-0.1, 3.6), ylab="Probability density",
     yaxt="n") # SLR (MCMC heteroskedastic calib.) accounting for runoff changes

lines(c(mean(proj2100DamsRes),mean(proj2100DamsRes)), c(-1,1.73), col=mygreen[4], lty=3,lwd=1.5)
lines(c(minus1sig,minus1sig), c(-1,1.28), col=mygreen[2], lty=3,lwd=1.5)
lines(c(plus1sig,plus1sig), c(-1,0.8), col=mygreen[6], lty=3,lwd=1.5)

points(minus1sig, 1.28, pch=21, bg=mygreen[2])
points(mean(proj2100DamsRes), 1.73, pch=21, bg=mygreen[4])
points(plus1sig, .8, pch=21, bg=mygreen[6])

text(mean(proj2100DamsRes), 1.9, "mean", col=mygreen[4], cex=1)
text((minus1sig - .6), 1.25, "-1 standard
deviation", col=mygreen[2], cex=1)
text((plus1sig + .5), 0.8, "+1 standard
deviation", col=mygreen[6], cex=1)

quantiles = c(0.05, 0.95)
add.hor.box(proj2100DamsRes, quantiles, width.size = 0.2, where.at = 2.1, tick.length = 0.05,
            line.width = 1, color = mygreen[4])
put.fig.letter("a.",font=2, cex=1)

legend("topright", ncol=1,c("Observations", "Current storm surge", "New 1:100-yr level"),
col = c("black", "skyblue","navy"),
pch = c(21, NA, 8),pt.bg = c("black", NA, NA),
lty = c(NA, 1, NA),lwd = c(NA, 1.5, NA),bty = 'n' , cex=1)

#-- Return level/ survival function plot ------------------------
par(mgp=c(1.5,.5,0),mar=c(4, 3, 3.5, 2))
plot.sf(coredata(year.res.max)/100, pch = 21, bg = "black",
        ylab = "Survival function [1 - cdf]",
        xlab = "Return Level (m)",sub="San Francisco Bay",
        yaxt = "n", yaxs = 'i',
        ylim = c(10^-6, 10^0+1), 
        xlim = c(min(mean.stormSLR.2100HET), 4))

abline(v=mean.stormSLR.2100HET[year100prob], lty=2, col=mygreen[4])
abline(v=Psig1.ss.slr[year100prob], lty=2, col=mygreen[6])
abline(v=Msig1.ss.slr.mil[msig.year100prob], lty=2, col=mygreen[2])
axis(2, at=10^(-6:-3), label=parse(text=paste("10^", -6:-3, sep="")))

lines(SSfit_q_year/100, 1-quant[1:SSL], type="l", col = "skyblue", lwd=1.5)
lines(mean.ss.slr.mil, 1-quant[1:SSL], type="l",col=mygreen[4], lwd=1.5)
lines(Psig1.ss.slr.mil , 1-quant[1:SSL], type="l",col=mygreen[6], lwd=1.5)
lines(Msig1.ss.slr.mil, 1-quant[1:SSL], type="l",col=mygreen[2], lwd=1.5)
abline(h=0.01, lty=2) # add in the 1 in 100
text(3.5, 0.015, "1:100 level", cex=1)

points(mean.stormSLR.2100HET[year100prob], new.prob.plussig, pch=21, bg=mygreen[6])
points(mean.stormSLR.2100HET[year100prob], prob.mil.minussig, pch=21, bg=mygreen[2])
points(mean.stormSLR.2100HET[year100prob], 1-q[year100prob], pch=21, bg=mygreen[4])
points(mean.stormSLR.2100HET[year100prob], real.average.sig.prob , pch=8, col="navy")
put.fig.letter("b.",font=2, cex=1)
dev.off()

############################# Supp. Figure #3 ##################################

# Supp. Figure 3 is created in ArcGIS ArcMap. Refer to the tutorial on how this is done.

############################# Supp. Figure #4 ##################################

# Supp. Figure 4 is created in ArcGIS ArcMap. Refer to the tutorial on how this is done.

############################# Supp. Figure #5 ##################################

# Supp. Figure 5 is created in ArcGIS ArcMap. Refer to the tutorial on how this is done.

############################# Supp. Figure #6 ##################################
png(file="Figures/SuppFigures/S6_Fig.tif", family="Helvetica", units="in", width=5, height=8.6, pointsize=12, res=300)
par(mfrow=c(2,1), mgp=c(1.5,.5,0),mar=c(4, 3, 3.5, 2))

plot.sf(coredata(year.res.max)/100, pch = 21, bg = "black",
        ylab = "Survival function [1 - cdf]",
        xlab = "Return Level (m)",sub="San Francisco Bay",
        yaxt = "n", yaxs = 'i',
        ylim = c(10^-6, 10^0+1), 
        xlim = c(-0.5, 3.5))

abline(v=mean.stormSLR.2100HET[year100prob], lty=2, col=mygreen[4])
axis(2, at=10^(-6:-3), label=parse(text=paste("10^", -6:-3, sep="")))

lines(SSfit_q_year/100, 1-quant[1:SSL], type="l", col = "skyblue", lwd=1.5)
lines(mean.ss.slr.mil, 1-quant[1:SSL], type="l",col=mygreen[4], lwd=1.5)
lines(Psig1.ss.slr.mil , 1-quant[1:SSL], type="l",col=mygreen[6], lwd=1.5)
lines(Msig1.ss.slr.mil, 1-quant[1:SSL], type="l",col=mygreen[2], lwd=1.5)
abline(h=0.01, lty=2) # add in the 1 in 100
text(0.5, 0.02, "1:100 level", cex=1)

points(mean.stormSLR.2100HET[year100prob], new.prob.plussig, pch=21, bg=mygreen[6])
points(mean.stormSLR.2100HET[year100prob], prob.mil.minussig, pch=21, bg=mygreen[2])
points(mean.stormSLR.2100HET[year100prob], 1-q[year100prob], pch=21, bg=mygreen[4])
points(mean.stormSLR.2100HET[year100prob], real.average.sig.prob , pch=8, col="navy")
put.fig.letter("a.",font=2, cex=1)

legend("bottomleft", ncol=1,
c("Observations", "Current storm surge", "-1 standard deviation",
"Storm surge + mean SLR", "+1 standard deviation", "New 100-yr level", "Approximate 100-yr level"),
col = c("black", "skyblue", mygreen[2], mygreen[4], mygreen[6], "black", "blue"),
pch = c(21, NA, 21, 21, 21, 8, 8),
pt.bg = c("black", NA, mygreen[2], mygreen[4], mygreen[6], NA, NA),
lty = c(NA, 1, 1, 1, 1, NA, NA),
lwd = c(NA, 1.5, 1.5, 1.5, 1.5, NA, NA),
bty = 'n' , cex=0.85)

par(mgp=c(1.5,.5,0),mar=c(4, 3, 3.5, 2))
return.period = c(round(1-q[year100prob],2), round(average.sig.prob,6), round(real.average.sig.prob,6))
barplot(return.period, xlab="100yr level", ylab="Return period at 2.76m (decimals)",
names.arg=c("Mean 100yr level", "Approximate", "New 100yr level"), col=c(mygreen[4], "blue", "black"))
put.fig.letter("b.",font=2, cex=1)

dev.off()

################################ PRINT NUMBERS #################################
# #Print the 90% confidence interval, median, and mean for the year 2100
# (without consideration of Dams and Res:
print(round(median.slr.proj$sle[221]/100,2)) #median based on median parameters
print(round(mean(prob_proj2100/100),2)) #mean
print(round(het_95[221],2)) #95% quartile
print(round(het_5[221],2)) #5% quartile

# Print the value to account for changes in Dams and Reserviors in California:
round(DamsReservoirs,2)

# Current 100-yr storm surge values:
round(fit_q_year[year100prob]/100,2) #Current storm surge for San Fran.
round(mean.stormSLR.2100HET[year100prob],2) #Current storm surge + Mean SLR for San Fran.
round(heberger09.2100[year100prob],2) #Current storm surge + H09 SLR for San Fran.
round(num.range[max.returnL],2)  # New storm surge + SLR for San Fran. [Accounting for uncertainty]

#Print the range of return level values used for accounting for SLR uncertainty
round(num.range[1], 2) ; round(num.range[length(num.range)],2)

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

# Storm surge values and probabilities for the mean +/- sigma test
#100-yr storm surge
round(mean.stormSLR.2100HET[year100prob],2) #Current storm surge + Mean SLR for San Fran.
round(Msig1.ss.slr.mil[msig.year100prob],2) # mean minus 1 sigma storm surge
round(Psig1.ss.slr[year100prob],2) # mean plus 1 sigma storm surge

# Probability of occurence at the mean 100-yr storm surge return level
round(1-q[year100prob], 2) #mean; should be 0.01
round(prob.mil.minussig, 6) # minus 1 sigma
round(new.prob.plussig, 2) # plus 1 sigma

#Factor of increase by accounting for uncertainty
# From mean slr estimate to new estimate:
# 0.01 to:
round(real.average.sig.prob, 2)

#################################### END #######################################



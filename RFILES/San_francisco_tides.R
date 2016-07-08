###################################
# file: San_francisco_tides.R
###################################
# Author and copyright: Perry Oddo
# Pennsylvania State University
# poddo@psu.edu
# Edits by Kelsey Ruckert
# klr324@psu.edu
# Code written Nov. 2015
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOR IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
###################################
# GEV Analysis of tide gauge data from
# San Francisco, CA
#################################### 

# Clear environment and graphics
rm(list = ls())
graphics.off()

# Install and load required libraries
# install.packages("extRemes")
# install.packages("fExtremes")
# install.packages("ismev")
# install.packages("lubridate")
# install.packages("zoo")
library(extRemes)
library(fExtremes)
library(ismev)
library(lubridate)
library(zoo)

# Read in tide gauge data
data = read.table("Data/SFB_TG_1914-2014.txt", header = TRUE, sep = '\t')
data$sl <- data$sl * 100
data[,4] <- NULL
data[,4] <- NULL

# Define function to aggregate data by years
as.year <- function(x) {floor(as.numeric(as.yearmon(x)))}

# Format data - create Datetime column and reference columns for each year in dataset
# date.time column used to create unique index for each entry 
data$date2 <- as.Date(as.character(data$date), format="%Y%m%d")
data$date.time <- as.POSIXct(paste(data$date, data$time), format = "%Y%m%d %H:%M", "GMT")
data$year.id <- as.numeric(as.factor(with(data, as.year(paste(date2)))))

# Create zoo object for tide gauge data for entire time series, ordered by date.time
sl <- zoo(data$sl, order.by=data$date.time)

# Aggregate annual block maximas and mean values
year.max <- aggregate(sl, as.year, max)
year.mean <- aggregate(sl, as.year, mean)

# Create data frame for block maxima, means, and index years
annual <- data.frame(index(year.max), coredata(year.mean), coredata(year.max))
annual$year.id <- index(annual[,1])

# Find residuals of annual means
# Create data frame to index and match annual means based on year.id
annual.res <- merge(data[, c("year.id", "sl")], annual[, c("year.id", "coredata.year.mean.")])
annual.res$date <- as.year(data$date2)
annual.res$residual <- annual.res$sl - annual.res$coredata.year.mean.   # Residuals of Annual Means
annual.res$date.time <- data$date.time   # Unique index for residuals

# Create zoo obhect for residuals and aggregate block maximas
# year.res.zoo are 'detrended' sea levels
year.res.zoo <- zoo(annual.res$residual, order.by = annual.res$date.time)
year.res.max <- aggregate(year.res.zoo, as.year, max)

# Save object for MCMC analysis
MCMC_coredata <- coredata(year.res.max)
save(MCMC_coredata, file = "Workspace/coredata.RData")

### Fit GEV of residuals ###
year.res.max.fit <- fevd(coredata(year.res.max))   # extRemes package
year.res.max.fit2 <- gevFit(coredata(year.res.max))   # fExtremes package
year.res.max.fit3 <- gev.fit(coredata(year.res.max), show = FALSE)   # ismev package

# Print GEV estimates
print(year.res.max.fit2@fit$par.ests)
#xi         mu       beta 
#0.22964812 1.30177669 0.06732941  

# Determine return levels using maximum likelihood estimate (95% confidence interval):
year.return <- return.level(year.res.max.fit, return.period = 100000, do.ci = TRUE)

# Print estimates of 2:100000-year flood, with confidence intervals
print(year.return)

save.image(file = "Workspace/SFB_stormSurge.RData")
################################## END #############################################

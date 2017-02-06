#######################################################################
#
# plot_rangefn.R    September 2015
#
# Author: Kelsey Ruckert (klr324@psu.edu)
#
# Function that plots a range or a box and whisker bar on an existing plot
#
# To use this function, simply source this file:
#   source("plot_rangefn.R")
#
# Note: I wrote this code because I've been asked fairly often for
# code that does this type of plot.  The original code I had included
# a lot of things the regular user wouldn't need, so I wrote a simple
# plot() wrapper to share with people.  It's not fancy, but it does the
# job.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOR IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL, 
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# Function Name: plotrange
#                - plots a range
#
# Function Name: add.hor.box
#                - plots a blox plot
#
#######################################################################
# Function to plot a range on an existing plot
plotrange <- function(low.number, median, high.number, year=T, height=F, color){
  
  range <- c(low.number, median, high.number)
  
  if(year){
    plot.position <- c(year, year, year)
    
    points(plot.position[1:2], c(range[1], range[3]), col=color, cex=1, pch="-")
    points(plot.position[1], range[2], col=color, cex=1, pch=8)
    segments(plot.position[1], range[1], plot.position[3], range[3], lwd=2, col=color)
  }
  
  if(height){
    plot.position <- c(height, height, height)
    
    points(c(range[1], range[3]), plot.position[1:2], col=color, cex=1, pch="|")
    points(range[2], plot.position[1], col=color, cex=1, pch=8)
    segments(range[1], plot.position[1], range[3], plot.position[3], lwd=2, col=color)
  }
}

#Function to plot a blox plot on an existing plot with adding lines showing select probabilities
add.hor.box <- function(data.numbers, probabilities, width.size, where.at, tick.length, line.width, color){
  quants <- quantile(data.numbers, probabilities)
  boxplot(data.numbers, add=TRUE, horizontal=TRUE, axes=FALSE, outline=FALSE, col=color, boxwex=width.size, at=where.at)
  segments(x0 = quants, y0 = rep(where.at - tick.length, 2),
           x1 = quants, y1 = rep(where.at + tick.length, 2), col = c("blue", "blue"), lwd = line.width)
}
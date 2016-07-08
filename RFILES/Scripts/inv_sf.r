#######################################################################
#
# inv_sf.r    15 Nov 2015
#
# Author: Gregory Garner (ggg121@psu.edu)
#
# Function that estimates the inverse of a survival function of a given vector of data
# and a value. It will interpolate the probability of a value.
#
# To use this function, simply source this file:
#   source("inv_sf.r")
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# Function Name: inv.sf
#
# Note that the value must be within the range of x values otherwise NA is returned.
#
#######################################################################

inv.sf <- function(x, val) {
  
  # Is val outside the range?
  if(max(x) < val | min(x) > val){
    return(NA)
  }
  
  # Sort the vector and get its length
  sort.x <- sort(x)
  num.x <- length(x)
  sf <- seq(1,1/num.x, by=-1/num.x)
  
  # Determine the bounding values
  lb.i <- findInterval(val, sort.x)
  ub.i <- lb.i + 1
  
  # Note, if the lower bound is the last
  # element in the vector, return the 
  # SF value for the last item
  # (i.e. x[num.x] == val)
  if (lb.i == num.x) {
    return(1/num.x)
  }
  
  # Set the lower and upper bounds for interpolation
  lb.x <- sort.x[lb.i]
  ub.x <- sort.x[ub.i]
  lb.sf <- sf[lb.i]
  ub.sf <- sf[ub.i]
  
  # Fit the interpolation
  ((val - lb.x)/(ub.x - lb.x)) * (ub.sf - lb.sf) + lb.sf
  
}
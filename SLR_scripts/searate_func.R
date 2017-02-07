###################################################################################
#
#  -file = "searate_func.R"   Code written March 2014
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This function sources the Rahmstorf(2007) model to find the rate of change in
#       sea level per year.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOR IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# To use this function, simply source this file:
#   source("searate_func.R")
#
# INPUTS:
#   vector of parameter values:
#       parameters[1]: alpha; sensitivity of sea-level to temperature changes
#       parameters[2]: T_0; temperature when sea-level anomaly is zero
#
#   vector of temperatures:
#
###################################################################################

thermal = function(parameters, Temp){ #inputs are parameters and temperature

    # Extract parameter values
    a = parameters[1] # sensitivity of sea level to temperature changes
    Ti = parameters[2] # temperature when sea-level anomaly is zero

    # Estimate the rate of sea-level change each year, Equation (1)
    therm_rate=a*(Temp-Ti)

    # Return the rates of sea-level change each year.
    list(therm_rate)
}

#################################### END #########################################
#***************************************************************
# Performs a search for an appropriate dementia cutpoint 
# Reference data found in dem_calcs(Gross, et.al).R
#***************************************************************

#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "magrittr")

#Suppress warnings
options(warn = -1)

#---- Specify source files ----
source("sex-dementia_sim_parA.R")
source("dementia_incidence2000-2013.R")







#---- Search for dementia cut-point ----

#Function we are trying to optimize
dem_1000py <- function(Cslopeb, Cslopec, cutpoint, obs){
  C_gen <- function(Cslopeb, Cslopec, obs){
    knots = c(0, 20, 35)
    Cij <- vector(length = length(visit_times))
    for(j in 1:length(Cij)){
      t = visit_times[j]
      test_num = j - 1
      eps <- paste("eps", test_num, sep = "")
      Cij[j] = b00 + obs[, "z0i"] + b01*obs[, "sex"] + b02*obs[, "age0_c50"] + 
        b03*obs[, "U"] + obs[, eps] + 
        (b10a - Cslopeb)*knots[2]*(t >= knots[2]) + 
        (Cslopeb - Cslopec)*knots[3]*(t >= knots[3]) + 
        (obs[, "z1i"] + b11*obs[, "sex"] + b12*obs[, "age0_c50"] + 
           b13*obs[, "U"] + b10a*(t >= knots[1] & t< knots[2]) + 
           Cslopeb*(t >= knots[2] & t< knots[3]) +
           Cslopec*(t >= knots[3]))*t
    }
    return(Cij)
  }
  
}

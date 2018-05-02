#***************************************************************
# Performs a parameter search using a specified function and 
# specified file containing reference data
#***************************************************************

#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse")

#---- Specify source files ----
source("sex-demensia_sim_parA.R")
source("sex-demensia_sim_script.R")
source("life_table2014.R")

#---- Search Function ----


#***************************************************************
# Performs a parameter search using a specified function and 
# specified file containing reference data
#***************************************************************

#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse")

#---- Specify source files ----
source("sex-dementia_sim_parA.R")
source("sex-dementia_sim_script.R")
source("life_table2014.R")

#---- Search Function ----
survival <- function(obs){
  #Calculate survival times for each interval
  Sij <- vector(length = (length(visit_times) - 1))
  for(j in 1:length(Sij)){
    test_num = j - 1
    US <- paste("U", test_num, test_num + 1, sep = "")
    agec <- paste("age", test_num, "_c50", sep = "")
    slope <- paste("slope", test_num, test_num + 1, sep = "")
    C <- paste("Ci", test_num, sep = "")
    Sij[j] = -log(obs[US])/
      (lambda*exp(g1*obs["sex"] + g2*obs[agec] + g3*obs["U"] + 
                    g4*obs["sex"]*obs[agec] + g5*obs[slope] + 
                    g6*obs[C]))
  }
  return("Sij" = Sij)
}

test <- 

#---- Quantiles of Cij Distribution ----
#Looking for a reasonable dementia cut point
#Use 4 simulated datasets and find quantiles of baseline Cij
dem_cut <- replicate(4, sex_dem_sim()) %>% map("Ci0") %>% unlist() %>% 
  quantile(0.05)
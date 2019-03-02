#---- Package Loading, Options, and Seed ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("future.apply", "parallel")

set.seed(20190213)

#---- Source Files ----
source("RScripts/sex-dementia_sim_parA.R")    #The parameter file
source("RScripts/variable_names.R")           #Creates all the variable names
source("RScripts/sex-dementia_sim_script.R")  #The simulation script
source("RScripts/misc_custom_functions.R")    #Other functions needed

#---- Running the simulation in parallel----
runs = 1000
plan(multiprocess, workers = (detectCores() - 2))
sim_results <- future_replicate(runs, sex_dem_sim())

#---- Converting results to usable format ----
results_mat <- matrix(unlist(sim_results), 
                         nrow = runs, ncol = ncol(sim_results), byrow = TRUE)
  
results_mat %<>% as.data.frame() %>% 
  set_colnames(dimnames(sim_results)[[2]]) %>%
  saveRDS("Data/logIRRs_A_20190223")








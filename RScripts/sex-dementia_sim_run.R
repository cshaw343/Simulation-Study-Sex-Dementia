#---- Package Loading, Options, and Seed ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("future.apply", "parallel")

set.seed(20190213)

#---- Source Files ----
source("RScripts/sex-dementia_sim_parA.R")          #The parameter file
source("RScripts/sex-dementia_sim_data_gen.R")      #The data generation script
source("RScripts/sex-dementia_sim_data_analysis.R") #The data analysis script
source("RScripts/misc_custom_functions.R")          #Other functions needed

#---- Generating one cohort ----
data_gen() %>% saveRDS("Data/test_sim_results_A_20190310")

#---- Running the simulation in parallel----
start_time <- Sys.time()
runs = 100
plan(multiprocess, workers = (detectCores() - 4))
sim_IRRs <- future_replicate(runs, sex_dem_sim())

#Saving results
results_mat <- matrix(unlist(sim_IRRs), 
                      nrow = runs, ncol = ncol(sim_IRRs), byrow = TRUE)

results_mat %>% as.data.frame() %>% 
  set_colnames(dimnames(sim_IRRs)[[2]]) %>%
  saveRDS(here("Data", "logIRRs_A_20190304"))

Sys.time() - start_time









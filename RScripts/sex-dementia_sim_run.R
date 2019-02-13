#---- Package Loading, Options, and Seed ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("parallel")

set.seed(20190213)

#---- Source Files ----
source("RScripts/sex-dementia_sim_parA.R")    #The parameter file
source("RScripts/variable_names.R")           #Creates all the variable names
source("RScripts/sex-dementia_sim_script.R")  #The simulation script
source("RScripts/misc_custom_functions.R")    #Other functions needed

# #---- Running the simulation in parallel----
# use_cores <- detectCores() - 1                 #Use one less than number of available cores
# cl <- makeCluster(use_cores, type = "FORK")    #Make a cluster from your cores
#                                                #"FORKING" environment only works on Mac
# runs = 100
# #Storing the results of the simulation
# start <- proc.time()  #This times the code
# sim_results <- parSapply(cl, 1:runs, sex_dem_sim)
# #stop the cluster
# stopCluster(cl)
# proc.time() - start

#---- Running simulation on one core ----
runs = 100
sim_results <- replicate(runs, sex_dem_sim())

#---- Converting results to usable format ----
results_mat <- matrix(nrow = runs*num_obs, ncol = nrow(sim_results))
for(r in 1:nrow(sim_results)){
  results_mat[, r] <- unlist(sim_results[r, ])
}

results_mat %<>% as.data.frame() %>% 
  set_colnames(dimnames(sim_results)[[1]]) %>%
  saveRDS("Data/test_sim_results_A_20190213")






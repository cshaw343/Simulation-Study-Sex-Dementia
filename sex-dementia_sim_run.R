#---- Specify the simulation script ----
source("sex-dementia_sim_script.R")
source("sex-dementia_sim_script_testSlope.R")

#---- Specify the parameter file ----
source("sex-dementia_sim_parA.R")

#---- Running the simulation----
#Storing the results of the simulation
set.seed(0)
sim_results <- replicate(1, sex_dem_sim())
set.seed(0)
sim_results_compare <- replicate(1, sex_dem_sim_test())

#---- Looking at simulation results ----
#Finding the mean Cij by sex across all simulations
mean_Cij_sim <- sim_results[[2]]
mean_Cij_sim_compare <- sim_results_compare[[2]]

#---- Specify the parameter file ----
source("RScripts/sex-dementia_sim_parA.R")

#---- Create all variable names ----
source("RScripts/variable_names.R")

#---- Specify the simulation script ----
source("RScripts/sex-dementia_sim_script.R")

#---- Running the simulation----
#Storing the results of the simulation
sim_results <- replicate(1, sex_dem_sim())

#---- Looking at simulation results ----
sim_obs <- sim_results["obs", ]
#Finding the mean Cij by sex across all simulations
mean_Cij_sim <- sim_results["mean_Cij", ]


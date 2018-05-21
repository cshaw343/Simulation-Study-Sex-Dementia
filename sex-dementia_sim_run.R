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
mean_Cij_sim <- as_tibble(do.call(rbind, sim_results$mean_Cij)) 

mean_Cij_sim_compare <- as_tibble(do.call(rbind, sim_results_compare)) %>% 
  group_by(sex) %>% summarise_all(mean)

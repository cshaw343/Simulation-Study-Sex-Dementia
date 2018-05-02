#---- Specify the source file ----
source("sex-demensia_sim_script.R")

#---- Specify the parameter file ----
par_file <- "sex-demensia_sim_parA.R"
source(par_file)

#---- Running the simulation----
#Storing the results of the simulation
sim_results <- replicate(1, sex_dem_sim())

#---- Looking at simulation results ----
#Finding the mean Cij by sex across all simulations
mean_Cij_sim <- as_tibble(do.call(rbind, sim_results)) %>% group_by(sex) %>% 
  summarise_all(mean)